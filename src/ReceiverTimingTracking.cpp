/**
 * @file ReceiverTimingTracking.cpp
 * @brief Thread: read 2 sps ring → @ref GardnerTiming::ProcessBlock → queue 1 sps frames for phase DD / Viterbi.
 *
 * @details **ThreadMain** — Under mutex, waits for @c SPB complex pairs on @c pFilterBuf_, appends to
 * @c pendingI_/Q_, runs Gardner until a full @c kSymFrame can be pushed to @c framePool_ queue
 * (@c pushOneFrameToQueue). @c WaitPopSymbolFrame blocks consumers on @c cvFrameReady_. Optional
 * @c DEBUG_STRESS_BACKPRESSURE adds sleep to stress upstream FIFOs.
 */

// Uncomment to stress backpressure (slow Gardner consumer; expect [ReceiverResampler] FIFO-at-cap logs).
// #define DEBUG_STRESS_BACKPRESSURE

#include "ReceiverTimingTracking.h"
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <immintrin.h>
#include <sys/stat.h>
#include <thread>
extern mutex mtxfilethr;

namespace
{
void TakeEvenDebug(float* x, float* output, int length)
{
    for (int i = 0; i < length; i += 16)
    {
        __m256 a0 = _mm256_loadu_ps(x + i);
        __m256 a1 = _mm256_loadu_ps(x + i + 8);
        __m256 e0 = _mm256_shuffle_ps(a0, a1, _MM_SHUFFLE(2, 0, 2, 0));
        __m256d e0d = _mm256_castps_pd(e0);
        e0d = _mm256_permute4x64_pd(e0d, 0xD8);
        e0 = _mm256_castpd_ps(e0d);
        _mm256_storeu_ps(output + i / 2, e0);
    }
}
} // namespace

ReceiverTimingTracking::ReceiverTimingTracking()
{
    const int outMax = 2 * ReceiverInputBatchIQSymbols;
    oneSpsI_ = static_cast<float*>(_mm_malloc(sizeof(float) * static_cast<size_t>(outMax), 32));
    oneSpsQ_ = static_cast<float*>(_mm_malloc(sizeof(float) * static_cast<size_t>(outMax), 32));
}

ReceiverTimingTracking::~ReceiverTimingTracking()
{
    StopJoin();
    _mm_free(oneSpsI_);
    _mm_free(oneSpsQ_);
}

void ReceiverTimingTracking::ResetGardner(double nominalRatio, double kp, double ki, int lockAvg)
{
    gardner_.Reset(nominalRatio, kp, ki, lockAvg);
}

void ReceiverTimingTracking::Start(BufferFloat* filterRing, std::mutex* mtxFilterRing,
                                   std::condition_variable* cvFilterData, std::atomic<bool>* stopAll,
                                   double displayPeriodSec)
{
    StopJoin();
    pFilterBuf_ = filterRing;
    pMtxFilter_ = mtxFilterRing;
    pCvFilterData_ = cvFilterData;
    pStopAll_ = stopAll;
    displayPeriodSec_ = displayPeriodSec;
    pendingCount_ = 0;
    frameQHead_ = 0;
    frameQTail_ = 0;
    frameQCount_ = 0;
    threadRunning_ = true;
    thread_ = std::thread(&ReceiverTimingTracking::ThreadMain, this);
}

void ReceiverTimingTracking::StopJoin()
{
    threadRunning_ = false;
    if (thread_.joinable())
    {
        cvFrameReady_.notify_all();
        cvFrameSpace_.notify_all();
        if (pCvFilterData_)
            pCvFilterData_->notify_all();
        thread_.join();
    }
    pFilterBuf_ = nullptr;
    pMtxFilter_ = nullptr;
    pCvFilterData_ = nullptr;
    pStopAll_ = nullptr;
    {
        std::lock_guard<std::mutex> lk(mtxFrameQ_);
        frameQHead_ = frameQTail_ = frameQCount_ = 0;
    }
}

bool ReceiverTimingTracking::WaitPopSymbolFrame(float* dstI, float* dstQ, int nSym)
{
    assert(nSym == kSymFrame);
    std::unique_lock<std::mutex> lk(mtxFrameQ_);
    cvFrameReady_.wait(lk, [&] {
        return (pStopAll_ != nullptr && *pStopAll_) || frameQCount_ > 0 || !threadRunning_;
    });
    if (frameQCount_ == 0)
        return false;
    const FrameSlot& f = framePool_[frameQHead_];
    std::memcpy(dstI, f.i, sizeof(float) * static_cast<size_t>(kSymFrame));
    std::memcpy(dstQ, f.q, sizeof(float) * static_cast<size_t>(kSymFrame));
    frameQHead_ = (frameQHead_ + 1) % kFrameQueueDepth;
    frameQCount_--;
    lk.unlock();
    cvFrameSpace_.notify_one();
    return true;
}

bool ReceiverTimingTracking::pushOneFrameToQueue(const float* i, const float* q)
{
    std::unique_lock<std::mutex> lk(mtxFrameQ_);
    cvFrameSpace_.wait(lk, [&] {
        return (pStopAll_ != nullptr && *pStopAll_) || !threadRunning_ ||
               frameQCount_ < kFrameQueueDepth;
    });
    if ((pStopAll_ != nullptr && *pStopAll_) || !threadRunning_)
        return false;
    FrameSlot& slot = framePool_[frameQTail_];
    std::memcpy(slot.i, i, sizeof(float) * static_cast<size_t>(kSymFrame));
    std::memcpy(slot.q, q, sizeof(float) * static_cast<size_t>(kSymFrame));
    frameQTail_ = (frameQTail_ + 1) % kFrameQueueDepth;
    frameQCount_++;
    lk.unlock();
    cvFrameReady_.notify_one();
    return true;
}

/**
 * @brief Gardner input pump and framed 1 sps output queue.
 */
void ReceiverTimingTracking::ThreadMain()
{

    #ifdef WRITE_LOG_THR
	mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"Receiver Timing Tracking  Thread %d\n", gettid());
    fclose(fidthr);
    mtxfilethr.unlock();
    #endif

    alignas(32) float filterChunkI[ReceiverInputBatchIQSamples];
    alignas(32) float filterChunkQ[ReceiverInputBatchIQSamples];
    uint64_t samples_total = 0;
    const auto wall_start = std::chrono::steady_clock::now();
    bool prevLocked = false;
    auto lastStatusDisplay = wall_start;

    while (threadRunning_ && (pStopAll_ == nullptr || !*pStopAll_))
    {
        float* filterOutI = nullptr;
        float* filterOutQ = nullptr;
        {
            std::unique_lock<std::mutex> lk(*pMtxFilter_);
            pFilterBuf_->GetReadBuffer(filterOutI, filterOutQ, ReceiverInputBatchIQSamples);
            while ((filterOutI == nullptr) && (pStopAll_ == nullptr || !*pStopAll_) && threadRunning_)
            {
                pCvFilterData_->wait_for(lk, std::chrono::duration<double>(1e-3));
                pFilterBuf_->GetReadBuffer(filterOutI, filterOutQ, ReceiverInputBatchIQSamples);
            }
            if (pStopAll_ != nullptr && *pStopAll_)
                break;
            if (!threadRunning_)
                break;
            if (filterOutI == nullptr)
                break;
            assert(filterOutI != nullptr && filterOutQ != nullptr);
            std::memcpy(filterChunkI, filterOutI,
                        sizeof(float) * static_cast<size_t>(ReceiverInputBatchIQSamples));
            std::memcpy(filterChunkQ, filterOutQ,
                        sizeof(float) * static_cast<size_t>(ReceiverInputBatchIQSamples));
        }

#ifdef DEBUG_STRESS_BACKPRESSURE
        // Slow consumer: resampler output ring fills → internal FIFO hits kFifoMax → filter blocks on AlmostFull.
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
#endif

#ifdef BYPASS_GARDNER_DUMP_VITERBI_INPUT
        TakeEvenDebug(filterChunkI, oneSpsI_, ReceiverInputBatchIQSamples);
        TakeEvenDebug(filterChunkQ, oneSpsQ_, ReceiverInputBatchIQSamples);
        {
            static FILE* gVitBypassDump = nullptr;
            if (!gVitBypassDump)
            {
                mkdir("../data", 0755);
                gVitBypassDump = std::fopen("../data/viterbi_input_bypass.bin", "wb");
            }
            if (gVitBypassDump)
            {
                for (int i = 0; i < kSymFrame; ++i)
                {
                    float iq[2] = {oneSpsI_[i], oneSpsQ_[i]};
                    std::fwrite(iq, sizeof(float), 2, gVitBypassDump);
                }
                std::fflush(gVitBypassDump);
            }
        }
        pendingCount_ = kSymFrame;
#else
        const float* gInI = filterChunkI;
        const float* gInQ = filterChunkQ;
        const int gInLen = ReceiverInputBatchIQSamples;
#ifdef DEBUG_GARDNER_OUTPUTS
        if (gInLen > 0)
        {
            static FILE* gGardnerInputDump = nullptr;
            if (!gGardnerInputDump)
            {
                mkdir("../data", 0755);
                gGardnerInputDump = std::fopen("../data/gardner_input_iq.bin", "wb");
            }
            if (gGardnerInputDump)
            {
                for (int i = 0; i < gInLen; ++i)
                {
                    float iq[2] = {gInI[i], gInQ[i]};
                    std::fwrite(iq, sizeof(float), 2, gGardnerInputDump);
                }
                std::fflush(gGardnerInputDump);
            }
        }
#endif
        int nSym = 0;
        if (gInLen > 0)
        {
            nSym = gardner_.ProcessBlock(gInI, gInQ, gInLen, oneSpsI_, oneSpsQ_,
                                         2 * ReceiverInputBatchIQSymbols);
        }
        samples_total += static_cast<uint64_t>(ReceiverInputBatchIQSamples);
        const bool nowLocked = gardner_.IsLocked();
        locked_.store(nowLocked, std::memory_order_relaxed);
        if (nowLocked != prevLocked)
        {
            const double t_rate = static_cast<double>(samples_total) / SamplingFrequency;
            const double t_sim =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
            
            std::cout << "[GardnerTiming] "
                      << (nowLocked ? "\033[32mLOCKED\033[0m" : "\033[31mUNLOCKED\033[0m")
                      << " omegaSpan=" << gardner_.GetLastOmegaLockSpan()
                      << " thresh=" << gardner_.GetOmegaLockSpanThreshold()
                      << " t_sim=" << t_sim << " s"
                      << " t_rate=" << t_rate << " s"
                      << std::endl;
            
            prevLocked = nowLocked;
        }

        // Periodic status line (once per second) to see lock-related values without waiting for transitions.
        {
            const auto now = std::chrono::steady_clock::now();
            if (displayPeriodSec_ > 0.0 &&
                now - lastStatusDisplay >= std::chrono::duration<double>(displayPeriodSec_))
            {
                const double t_rate = static_cast<double>(samples_total) / SamplingFrequency;
                const double t_sim = std::chrono::duration<double>(now - wall_start).count();
                std::cout << "[GardnerTiming] status"
                          << " locked=" << (nowLocked ? 1 : 0)
                          << " omegaSpan=" << gardner_.GetLastOmegaLockSpan()
                          << " thresh=" << gardner_.GetOmegaLockSpanThreshold()
                          << " t_sim=" << t_sim << " s"
                          << " t_rate=" << t_rate << " s"
                          << std::endl;
                lastStatusDisplay = now;
            }
        }
        // Always produce symbols/frames, even when Gardner is not locked.
        // Downstream blocks can use IsLocked() to decide how to interpret status, but they should not stall.
        if (nSym > 0)
        {
            const int copyCount = std::min(nSym, kPendingCap - pendingCount_);
            std::copy(oneSpsI_, oneSpsI_ + copyCount, pendingI_ + pendingCount_);
            std::copy(oneSpsQ_, oneSpsQ_ + copyCount, pendingQ_ + pendingCount_);
            pendingCount_ += copyCount;
        }
#endif

        while (pendingCount_ >= kSymFrame)
        {
            if (pStopAll_ != nullptr && *pStopAll_)
                break;
            if (!pushOneFrameToQueue(pendingI_, pendingQ_))
                break;
            pendingCount_ -= kSymFrame;
            if (pendingCount_ > 0)
            {
                std::copy(pendingI_ + kSymFrame, pendingI_ + kSymFrame + pendingCount_, pendingI_);
                std::copy(pendingQ_ + kSymFrame, pendingQ_ + kSymFrame + pendingCount_, pendingQ_);
            }
        }

        {
            std::unique_lock<std::mutex> lk(*pMtxFilter_);
            pFilterBuf_->AdvancePtrRd(ReceiverInputBatchIQSamples);
        }
        if (pCvFilterData_ != nullptr)
            pCvFilterData_->notify_one();
    }
}
