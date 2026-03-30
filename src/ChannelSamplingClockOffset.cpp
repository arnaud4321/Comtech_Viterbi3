#include "ChannelSamplingClockOffset.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>

void ChannelSamplingClockOffset::ScoRingBuffer::PushBlock(const float* inI, const float* inQ, int n)
{
    for (int i = 0; i < n; ++i)
    {
        const long long a = absWrite + static_cast<long long>(i);
        ringI[static_cast<int>(a) & kMask] = inI[i];
        ringQ[static_cast<int>(a) & kMask] = inQ[i];
    }
    absWrite += static_cast<long long>(n);
}

long long ChannelSamplingClockOffset::ScoRingBuffer::OldestAbs() const
{
    return absWrite - static_cast<long long>(kSize);
}

float ChannelSamplingClockOffset::ScoRingBuffer::GetI(long long absIdx) const
{
    return ringI[static_cast<int>(absIdx) & kMask];
}

float ChannelSamplingClockOffset::ScoRingBuffer::GetQ(long long absIdx) const
{
    return ringQ[static_cast<int>(absIdx) & kMask];
}

float ChannelSamplingClockOffset::ScoRingBuffer::InterpLagrange4I(double t) const
{
    const long long k = static_cast<long long>(std::floor(t));
    const double mu = t - static_cast<double>(k);
    const float x0 = GetI(k - 1);
    const float x1 = GetI(k + 0);
    const float x2 = GetI(k + 1);
    const float x3 = GetI(k + 2);

    const double c0 = -mu * (mu - 1.0) * (mu - 2.0) / 6.0;
    const double c1 = (mu + 1.0) * (mu - 1.0) * (mu - 2.0) / 2.0;
    const double c2 = -(mu + 1.0) * mu * (mu - 2.0) / 2.0;
    const double c3 = (mu + 1.0) * mu * (mu - 1.0) / 6.0;

    return static_cast<float>(c0 * x0 + c1 * x1 + c2 * x2 + c3 * x3);
}

float ChannelSamplingClockOffset::ScoRingBuffer::InterpLagrange4Q(double t) const
{
    const long long k = static_cast<long long>(std::floor(t));
    const double mu = t - static_cast<double>(k);
    const float x0 = GetQ(k - 1);
    const float x1 = GetQ(k + 0);
    const float x2 = GetQ(k + 1);
    const float x3 = GetQ(k + 2);

    const double c0 = -mu * (mu - 1.0) * (mu - 2.0) / 6.0;
    const double c1 = (mu + 1.0) * (mu - 1.0) * (mu - 2.0) / 2.0;
    const double c2 = -(mu + 1.0) * mu * (mu - 2.0) / 2.0;
    const double c3 = (mu + 1.0) * mu * (mu - 1.0) / 6.0;

    return static_cast<float>(c0 * x0 + c1 * x1 + c2 * x2 + c3 * x3);
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
    const double eps = totalPpm_ * 1.0e-6;
    step_ = 1.0 / (1.0 + eps);
}

void ChannelSamplingClockOffset::Start(BufferShort* outputBuffer, std::mutex* outputMutex,
                                       std::condition_variable* cvSpace,
                                       std::condition_variable* cvData, bool* stopAll)
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
    cvInSpace_.wait(lk, [&] {
        return (pStopAll_ && *pStopAll_) || inQueue_.size() < kInputQueueDepth;
    });
    if (pStopAll_ && *pStopAll_)
        return;
    inQueue_.emplace_back(interleaved, interleaved + nFloats);
    cvInData_.notify_one();
}

void ChannelSamplingClockOffset::ThreadMain()
{
    alignas(32) float deintI[kComplexPerBatch];
    alignas(32) float deintQ[kComplexPerBatch];
    alignas(32) float scoChunkI[kScoOutMax];
    alignas(32) float scoChunkQ[kScoOutMax];

    ScoRingBuffer scoRing;
    bool scoInit = false;
    double tScoAbs = 0.0;

    const long long latencySamples = static_cast<long long>(kComplexPerBatch) * 64LL;
    assert(latencySamples + 16 < ScoRingBuffer::kSize);

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

        const double step = step_;

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
                    scoChunkI[n] = scoRing.InterpLagrange4I(t);
                    scoChunkQ[n] = scoRing.InterpLagrange4Q(t);
                    tScoAbs += step;
                }
            }
        }

        if (scoOutLen <= 0 || pOutBuf_ == nullptr || pOutMtx_ == nullptr)
            continue;

        const int nShorts = 2 * scoOutLen;
        std::unique_lock<std::mutex> lkOut(*pOutMtx_);
        while (pStopAll_ != nullptr && !*pStopAll_ && pOutBuf_->AlmostFull())
            pCvSpace_->wait(lkOut);
        if (pStopAll_ != nullptr && *pStopAll_)
            continue;

        short* out = pOutBuf_->GetWriteBuffer(nShorts);
        FloatBatchToShortsInterleaved(scoChunkI, scoChunkQ, scoOutLen, out);
        pOutBuf_->AdvancePtrWr(nShorts);
        lkOut.unlock();
        pCvData_->notify_one();
    }
}
