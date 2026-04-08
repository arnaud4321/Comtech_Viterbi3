/**
 * @file ReceiverResampler.cpp
 * @brief Worker thread: drain @c oBufferFilter → internal FIFO → @ref Resampler::CreateOutputs → downstream ring.
 *
 * @details **ThreadMain** — While running: refill FIFO from input ring up to @c kFifoMax (else block on
 * @c inCvData_ — backpressure). Calls @c resampler_ with atomic @c advance_ from
 * @ref ReceiverResampler::UpdateFromSymbolRateHz. Pushes output blocks to @c outRing_, notifies @c outCvData_.
 * @c UpdateFromSymbolRateHz sets @f$a \approx F_s^{\mathrm{in}}/(2R_s)@f$ clamped to @f$[0.5,2]@f$.
 */

#include "ReceiverResampler.h"
#include "ConsoleAlert.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
extern mutex mtxfilethr;

ReceiverResampler::ReceiverResampler() = default;

ReceiverResampler::~ReceiverResampler()
{
    StopJoin();
}

void ReceiverResampler::Start(BufferFloat* inRing, std::mutex* inMtx, std::condition_variable* inCvData,
                              BufferFloat* outRing, std::mutex* outMtx, std::condition_variable* outCvData,
                              std::condition_variable* outCvSpace, std::atomic<bool>* stopAll)
{
    StopJoin();
    inRing_ = inRing;
    inMtx_ = inMtx;
    inCvData_ = inCvData;
    outRing_ = outRing;
    outMtx_ = outMtx;
    outCvData_ = outCvData;
    outCvSpace_ = outCvSpace;
    stopAll_ = stopAll;
    frac_ = 0.0;
    advance_.store(1.0, std::memory_order_relaxed);
    fifoI_.clear();
    fifoQ_.clear();
    fifoI_.reserve(static_cast<size_t>(std::min(kFifoMax, SPB * 64)));
    fifoQ_.reserve(static_cast<size_t>(std::min(kFifoMax, SPB * 64)));
    fifoRd_ = 0;
    running_ = true;
    thread_ = std::thread(&ReceiverResampler::ThreadMain, this);
}

void ReceiverResampler::StopJoin()
{
    running_ = false;
    if (thread_.joinable())
    {
        if (inCvData_)
            inCvData_->notify_all();
        if (outCvData_)
            outCvData_->notify_all();
        if (outCvSpace_)
            outCvSpace_->notify_all();
        thread_.join();
    }
    inRing_ = nullptr;
    inMtx_ = nullptr;
    inCvData_ = nullptr;
    outRing_ = nullptr;
    outMtx_ = nullptr;
    outCvData_ = nullptr;
    outCvSpace_ = nullptr;
    stopAll_ = nullptr;}

void ReceiverResampler::UpdateFromSymbolRateHz(double symbolRateHz)
{
    if (symbolRateHz <= 0.0)
        return;
    const double FsIn = SamplingFrequency;
    const double FsOut = 2.0 * symbolRateHz;
    double a = FsIn / FsOut;
    // Keep within a safe range.
    a = std::max(0.5, std::min(2.0, a));
    advance_.store(a, std::memory_order_relaxed);
}

/**
 * @brief Resampler loop: bounded FIFO, Lagrange outputs, write to @c outRing_.
 */
void ReceiverResampler::ThreadMain()
{

    #ifdef WRITE_LOG_THR
	mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"Receiver Resampler  Thread %d\n", gettid());
    fclose(fidthr);
    mtxfilethr.unlock();
    #endif

    constexpr int kReadChunk = SPB;
    constexpr int kOutBlockMax = SPB * 2; // worst-case upsampling (Advance=0.5)
    auto lastFifoFullLog = std::chrono::steady_clock::now();
    auto lastOutRingFullLog = std::chrono::steady_clock::now();
    auto lastOutCommitBlockLog = std::chrono::steady_clock::now();

    alignas(32) float outI[kOutBlockMax];
    alignas(32) float outQ[kOutBlockMax];

    while (running_ && stopAll_ && !*stopAll_)
    {
        // 1) Refill internal FIFO from input ring if needed (bounded by kFifoMax → backpressure on inRing_).
        {
            std::unique_lock<std::mutex> lk(*inMtx_);
            while (running_ && !*stopAll_ &&
                   (static_cast<int>(fifoI_.size() - fifoRd_) < (kReadChunk + kOverlap)) &&
                   (inRing_->GetSizeInBuffer() < kReadChunk))
            {
                inCvData_->wait_for(lk, std::chrono::milliseconds(1));
            }
            if (!running_ || *stopAll_)
                break;

            while (inRing_->GetSizeInBuffer() >= kReadChunk)
            {
                const size_t fifoUsed = fifoI_.size() - fifoRd_;
                if (fifoUsed + static_cast<size_t>(kReadChunk) > static_cast<size_t>(kFifoMax))
                {
                    // Backpressure: do not AdvancePtrRd; input ring stays full → filter thread waits on AlmostFull.
                    const auto now = std::chrono::steady_clock::now();
                    if (now - lastFifoFullLog >= std::chrono::seconds(1))
                    {
                        CONSOLE_ALERT_STMT(
                            std::cout << ConsoleAlert::kRedOpen << "[ReceiverResampler] internal FIFO at cap "
                                      << fifoUsed << "/" << kFifoMax
                                      << " samples; pausing drain of matched-filter ring (pending "
                                      << inRing_->GetSizeInBuffer() << " samples)" << ConsoleAlert::kReset
                                      << std::endl;);
                        lastFifoFullLog = now;
                    }
                    break;
                }

                float* inI = nullptr;
                float* inQ = nullptr;
                inRing_->GetReadBuffer(inI, inQ, kReadChunk);
                if (!inI)
                    break;
                fifoI_.insert(fifoI_.end(), inI, inI + kReadChunk);
                fifoQ_.insert(fifoQ_.end(), inQ, inQ + kReadChunk);
                inRing_->AdvancePtrRd(kReadChunk);
                inCvData_->notify_one();
            }
        }

        const size_t avail = (fifoI_.size() > fifoRd_) ? (fifoI_.size() - fifoRd_) : 0u;
        if (avail < 4u)
            continue;

        // 2) Determine output capacity.
        int outMax = 0;
        {
            std::unique_lock<std::mutex> lk(*outMtx_);
            const int fill = outRing_->GetSizeInBuffer();
            const int cap = outRing_->GetBufferSize() - 1;
            const int free = std::max(0, cap - fill);
            outMax = std::min(kOutBlockMax, free);
            while (outMax <= 0 && running_ && !*stopAll_)
            {
                const auto nowOr = std::chrono::steady_clock::now();
                if (nowOr - lastOutRingFullLog >= std::chrono::seconds(1))
                {
                    lastOutRingFullLog = nowOr;
                    CONSOLE_ALERT_STMT(std::cout << ConsoleAlert::kRedOpen
                                                 << "[ReceiverResampler] output float ring saturated (no free "
                                                    "slots); fill="
                                                 << outRing_->GetSizeInBuffer() << "/" << cap << ConsoleAlert::kReset
                                                 << std::endl;);
                }
                outCvSpace_->wait_for(lk, std::chrono::milliseconds(1));
                const int fill2 = outRing_->GetSizeInBuffer();
                const int free2 = std::max(0, cap - fill2);
                outMax = std::min(kOutBlockMax, free2);
            }
            if (!running_ || *stopAll_)
                break;
        }
        if (outMax <= 0)
            continue;

        // 3) Resample from FIFO to output block.
        const double a = advance_.load(std::memory_order_relaxed);
        unsigned int start = 0;
        const unsigned int produced = resampler_.CreateOutputs(
            fifoI_.data() + fifoRd_, fifoQ_.data() + fifoRd_,
            outI, outQ,
            a, start, frac_,
            static_cast<unsigned int>(avail),
            static_cast<unsigned int>(outMax));

        // Consume input samples from FIFO, keep kOverlap lookahead implicitly (resampler stops at Length-3).
        const size_t consumed = static_cast<size_t>(start);
        fifoRd_ += consumed;
        if (fifoRd_ > static_cast<size_t>(SPB * 8))
        {
            fifoI_.erase(fifoI_.begin(), fifoI_.begin() + static_cast<long>(fifoRd_));
            fifoQ_.erase(fifoQ_.begin(), fifoQ_.begin() + static_cast<long>(fifoRd_));
            fifoRd_ = 0;
        }

        if (produced == 0)
            continue;

        // 4) Commit outputs to ring.
        {
            std::unique_lock<std::mutex> lk(*outMtx_);
            while (!*stopAll_ && !outRing_->CanCommitWrite(static_cast<int>(produced)))
            {
                const auto nowOc = std::chrono::steady_clock::now();
                if (nowOc - lastOutCommitBlockLog >= std::chrono::seconds(1))
                {
                    lastOutCommitBlockLog = nowOc;
                    CONSOLE_ALERT_STMT(std::cout << ConsoleAlert::kRedOpen
                                                 << "[ReceiverResampler] output float ring cannot commit " << produced
                                                 << " samples (backpressure)" << ConsoleAlert::kReset << std::endl;);
                }
                outCvSpace_->wait(lk);
            }
            if (*stopAll_)
                break;
            float* wI = nullptr;
            float* wQ = nullptr;
            outRing_->GetWriteBuffer(wI, wQ, static_cast<int>(produced));
            std::copy(outI, outI + produced, wI);
            std::copy(outQ, outQ + produced, wQ);
            outRing_->AdvancePtrWr(static_cast<int>(produced));
        }
        outCvData_->notify_one();
    }
}

