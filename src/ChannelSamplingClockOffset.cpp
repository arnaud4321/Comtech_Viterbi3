/**
 * @file ChannelSamplingClockOffset.cpp
 * @brief Async resampling thread: float IQ in → Lagrange ring at @f$F_s(1+\mathrm{ppm}\cdot10^{-6})@f$ → shorts to @c OutputBuffer.
 *
 * @details **EnqueueNoisyInterleaved** — After NCO rotation: pushes @c TxOutputBatchSize floats into the input
 * queue. **ThreadMain** — Ring + @c tScoAbs with @c step = @f$1/(1+\varepsilon)@f$, Lagrange-4 outputs,
 * latency warmup, safe-band checks, then @c FloatBatchToShortsInterleaved → @c BufferShort.
 * **Configure/UpdateTotalPpm** — Set @f$\varepsilon@f$ and the atomic @c step_.
 *
 * Simulations may run for weeks or months; @ref ScoRingBuffer::RewrapForLongRun keeps absolute indices and
 * the Lagrange time cursor in a moderate range (see header).
 */

#include "ChannelSamplingClockOffset.h"
#include "ConsoleAlert.h"
#include "Lagrange4Simd.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstring>
#include <iostream>
extern mutex mtxfilethr;

void ChannelSamplingClockOffset::ScoRingBuffer::PushBlock(const float* inI, const float* inQ, int n)
{
    for (int i = 0; i < n; ++i)
    {
        const long long a = absWrite + static_cast<long long>(i);
        const size_t idx = static_cast<size_t>(a) & static_cast<size_t>(kMask);
        ringI[idx] = inI[i];
        ringQ[idx] = inQ[i];
    }
    absWrite += static_cast<long long>(n);
}

long long ChannelSamplingClockOffset::ScoRingBuffer::OldestAbs() const
{
    return absWrite - static_cast<long long>(kSize);
}

float ChannelSamplingClockOffset::ScoRingBuffer::GetI(long long absIdx) const
{
    return ringI[static_cast<size_t>(absIdx) & static_cast<size_t>(kMask)];
}

float ChannelSamplingClockOffset::ScoRingBuffer::GetQ(long long absIdx) const
{
    return ringQ[static_cast<size_t>(absIdx) & static_cast<size_t>(kMask)];
}

void ChannelSamplingClockOffset::ScoRingBuffer::RewrapForLongRun(double* tCursor)
{
    // Keep absWrite and the Lagrange time base small: double ULP on fractional mu, safe indexing via mask.
    // Subtracting one full ring length preserves (index mod kSize) for live samples (kSize is power of 2).
    while (absWrite >= 2LL * static_cast<long long>(kSize))
    {
        absWrite -= static_cast<long long>(kSize);
        if (tCursor != nullptr)
            *tCursor -= static_cast<double>(kSize);
    }
}

float ChannelSamplingClockOffset::ScoRingBuffer::InterpLagrange4I(double t) const
{
    const long long k = static_cast<long long>(std::floor(t));
    const double mu = t - static_cast<double>(k);
    return Lagrange4Simd::EvalSet_ps(GetI(k - 1), GetI(k + 0), GetI(k + 1), GetI(k + 2), Lagrange4Simd::Coeffs_ps(mu));
}

float ChannelSamplingClockOffset::ScoRingBuffer::InterpLagrange4Q(double t) const
{
    const long long k = static_cast<long long>(std::floor(t));
    const double mu = t - static_cast<double>(k);
    return Lagrange4Simd::EvalSet_ps(GetQ(k - 1), GetQ(k + 0), GetQ(k + 1), GetQ(k + 2), Lagrange4Simd::Coeffs_ps(mu));
}

void ChannelSamplingClockOffset::FloatBatchToShortsInterleaved(const float* iPtr, const float* qPtr,
                                                                 int nComplex, short* dst)
{
    for (int i = 0; i < nComplex; ++i)
    {
        int vi = static_cast<int>(std::lrintf(iPtr[i]));
        int vq = static_cast<int>(std::lrintf(qPtr[i]));
        vi = std::max(-32768, std::min(32767, vi));
        vq = std::max(-32768, std::min(32767, vq));
        dst[2 * i] = static_cast<short>(vi);
        dst[2 * i + 1] = static_cast<short>(vq);
    }
}

ChannelSamplingClockOffset::ChannelSamplingClockOffset() = default;

ChannelSamplingClockOffset::~ChannelSamplingClockOffset()
{
    StopJoin();
}

void ChannelSamplingClockOffset::Configure(double totalPpm)
{
    totalPpm_ = totalPpm;
    UpdateTotalPpm(totalPpm);
}

void ChannelSamplingClockOffset::UpdateTotalPpm(double totalPpm)
{
    totalPpm_ = totalPpm;
    const double eps = totalPpm_ * 1.0e-6;
    step_.store(1.0 / (1.0 + eps), std::memory_order_relaxed);
}

void ChannelSamplingClockOffset::Start(BufferShort* outputBuffer, std::mutex* outputMutex,
                                       std::condition_variable* cvSpace,
                                       std::condition_variable* cvData, std::atomic<bool>* stopAll)
{
    StopJoin();
    pOutBuf_ = outputBuffer;
    pOutMtx_ = outputMutex;
    pCvSpace_ = cvSpace;
    pCvData_ = cvData;
    pStopAll_ = stopAll;
    threadRunning_ = true;
    thread_ = std::thread(&ChannelSamplingClockOffset::ThreadMain, this);
}

void ChannelSamplingClockOffset::StopJoin()
{
    threadRunning_ = false;
    if (thread_.joinable())
    {
        cvInData_.notify_all();
        cvInSpace_.notify_all();
        thread_.join();
    }
    pOutBuf_ = nullptr;
    pOutMtx_ = nullptr;
    pCvSpace_ = nullptr;
    pCvData_ = nullptr;
    pStopAll_ = nullptr;
    std::lock_guard<std::mutex> lk(mtxIn_);
    inQueue_.clear();
}

void ChannelSamplingClockOffset::EnqueueNoisyInterleaved(const float* interleaved, int nFloats)
{
    if (!pStopAll_ || !threadRunning_)
        return;
    std::unique_lock<std::mutex> lk(mtxIn_);
    static auto lastScoInFifoSatLog = std::chrono::steady_clock::now();
    if (inQueue_.size() >= kInputQueueDepth)
    {
        const auto nowSat = std::chrono::steady_clock::now();
        if (nowSat - lastScoInFifoSatLog >= std::chrono::seconds(1))
        {
            lastScoInFifoSatLog = nowSat;
            CONSOLE_ALERT_STMT(std::cout << ConsoleAlert::kRedOpen << "[SCO] input FIFO saturated; depth="
                                         << inQueue_.size() << "/" << kInputQueueDepth
                                         << " batches (producer blocked)" << ConsoleAlert::kReset << std::endl;);
        }
    }
    cvInSpace_.wait(lk, [&] {
        return (pStopAll_ && *pStopAll_) || inQueue_.size() < kInputQueueDepth;
    });
    if (pStopAll_ && *pStopAll_)
        return;
    inQueue_.emplace_back(interleaved, interleaved + nFloats);
    cvInData_.notify_one();
}

/**
 * @brief SCO worker: dequeue float batches, resample ring, write shorts to @c OutputBuffer.
 */
void ChannelSamplingClockOffset::ThreadMain()
{

    #ifdef WRITE_LOG_THR
	mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"ChannelSamplingClockOffset  Thread %d\n", gettid());
    fclose(fidthr);
    mtxfilethr.unlock();
    #endif


    alignas(32) float deintI[kComplexPerBatch];
    alignas(32) float deintQ[kComplexPerBatch];
    alignas(32) float scoChunkI[kScoOutMax];
    alignas(32) float scoChunkQ[kScoOutMax];

    ScoRingBuffer scoRing;
    bool scoInit = false;
    double tScoAbs = 0.0;

    const long long latencySamples = static_cast<long long>(kComplexPerBatch) * 64LL;
    assert(latencySamples + 16 < ScoRingBuffer::kSize);

    auto lastShortOutSatLog = std::chrono::steady_clock::now();

    for (;;)
    {
        std::vector<float> batch;
        {
            std::unique_lock<std::mutex> lk(mtxIn_);
            cvInData_.wait(lk, [&] {
                return !inQueue_.empty() || !threadRunning_ ||
                       (pStopAll_ != nullptr && *pStopAll_);
            });
            if (inQueue_.empty())
            {
                if (!threadRunning_ || (pStopAll_ != nullptr && *pStopAll_))
                    break;
                continue;
            }
            batch = std::move(inQueue_.front());
            inQueue_.pop_front();
            lk.unlock();
            cvInSpace_.notify_one();
        }

        if (static_cast<int>(batch.size()) != TxOutputBatchSize)
            continue;

        const int nC = TxOutputBatchSize / 2;
        for (int i = 0; i < nC; ++i)
        {
            deintI[i] = batch[static_cast<size_t>(2 * i)];
            deintQ[i] = batch[static_cast<size_t>(2 * i + 1)];
        }

        scoRing.PushBlock(deintI, deintQ, nC);

        const double step = step_.load(std::memory_order_relaxed);

        int scoOutLen = 0;
        if (!scoInit)
        {
            if (scoRing.absWrite > latencySamples + 4)
            {
                tScoAbs = static_cast<double>(scoRing.absWrite - latencySamples);
                scoInit = true;
            }
        }
        if (scoInit)
        {
            const double newestSafe = static_cast<double>(scoRing.absWrite - 4);
            const double oldestSafe = static_cast<double>(scoRing.OldestAbs() + 2);

            if ((tScoAbs < oldestSafe) || (tScoAbs > newestSafe))
            {
                scoInit = false;
                scoOutLen = 0;
            }
            else
            {
                const double targetLag = static_cast<double>(latencySamples);
                const double desiredTAfter = static_cast<double>(scoRing.absWrite) - targetLag;
                const double nOutReal = (desiredTAfter - tScoAbs) / step;
                int nOut = static_cast<int>(std::llround(nOutReal));
                if (nOut < 0)
                    nOut = 0;

                int nOutByNewest = static_cast<int>(std::floor((newestSafe - tScoAbs) / step)) + 1;
                if (nOutByNewest < 0)
                    nOutByNewest = 0;
                nOut = std::min(nOut, nOutByNewest);
                nOut = std::min(nOut, kScoOutMax);
                scoOutLen = nOut;

                for (int n = 0; n < scoOutLen; ++n)
                {
                    const double t = tScoAbs;
                    const long long k = static_cast<long long>(std::floor(t));
                    const double mu = t - static_cast<double>(k);
                    Lagrange4Simd::EvalIQ_ps(scoRing.GetI(k - 1), scoRing.GetI(k + 0), scoRing.GetI(k + 1),
                                             scoRing.GetI(k + 2), scoRing.GetQ(k - 1), scoRing.GetQ(k + 0),
                                             scoRing.GetQ(k + 1), scoRing.GetQ(k + 2), mu, &scoChunkI[n],
                                             &scoChunkQ[n]);
                    tScoAbs += step;
                }
            }
        }

        // Multi-week / multi-month runs: bound absWrite and tScoAbs (when active) without changing ring physics.
        scoRing.RewrapForLongRun(scoInit ? &tScoAbs : nullptr);

        if (scoOutLen <= 0 || pOutBuf_ == nullptr || pOutMtx_ == nullptr)
            continue;

        const int nShorts = 2 * scoOutLen;
        std::unique_lock<std::mutex> lkOut(*pOutMtx_);
        while (pStopAll_ != nullptr && !*pStopAll_ && pOutBuf_->AlmostFull())
        {
            const auto nowSat = std::chrono::steady_clock::now();
            if (nowSat - lastShortOutSatLog >= std::chrono::seconds(1))
            {
                lastShortOutSatLog = nowSat;
                CONSOLE_ALERT_STMT(std::cout << ConsoleAlert::kRedOpen
                                             << "[SCO] AWGN output short ring almost full; fill="
                                             << pOutBuf_->GetSizeInBuffer() << "/" << pOutBuf_->GetBufferSize() - 1
                                             << ConsoleAlert::kReset << std::endl;);
            }
            pCvSpace_->wait(lkOut);
        }
        if (pStopAll_ != nullptr && *pStopAll_)
            continue;

        short* out = pOutBuf_->GetWriteBuffer(nShorts);
        FloatBatchToShortsInterleaved(scoChunkI, scoChunkQ, scoOutLen, out);
        pOutBuf_->AdvancePtrWr(nShorts);
        lkOut.unlock();
        pCvData_->notify_one();
    }
}
