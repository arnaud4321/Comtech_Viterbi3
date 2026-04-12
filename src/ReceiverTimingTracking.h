/**
 * @file ReceiverTimingTracking.h
 * @brief Thread: Gardner timing recovery from 2 sps stream → framed 1 sps symbols for downstream PLL/Viterbi.
 */
#pragma once

#include "BufferFloat.h"
#include "definitions.h"
#include "GardnerTiming.h"
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

/**
 * @brief Gardner PI loop and update cadence (JSON section @c TimingTracking).
 *
 * @c UpdatePeriodSymbols is the number of output symbols between each @f$\omega@f$ update
 * (@ref GardnerTiming::Reset @a updatePeriod).
 */
struct TimingTrackingConfig
{
    double NominalOmega = 2.0;
    double Kp = 1.0e-3;
    double Ki = 1.0e-5;
    int UpdatePeriodSymbols = 64;
};

/**
 * @brief Thread: 2 sps stream → @ref GardnerTiming → framed 1 sps symbols for PLL and @ref Viterbi.
 *
 * @details Reads @ref BufferFloat from the stage upstream of Gardner (after optional NCO/AGC), pushes
 * samples through @ref GardnerTiming::ProcessBlock, and queues fixed-length frames of
 * @c ReceiverInputBatchIQSymbols complex symbols. The Viterbi manager and @ref ReceiverPhaseTrackingDD
 * block on @ref WaitPopSymbolFrame.
 *
 * **Debug:** uncomment @c DEBUG_STRESS_BACKPRESSURE in @c ReceiverTimingTracking.cpp to inject an
 * artificial per-batch delay (stress backpressure into the resampler / filter ring).
 */
class ReceiverTimingTracking
{
public:
    ReceiverTimingTracking();
    ~ReceiverTimingTracking();

    void ResetGardner(double nominalRatio, double kp, double ki, int lockAvg);

    void Start(BufferFloat* filterRing, std::mutex* mtxFilterRing,
               std::condition_variable* cvFilterData, std::atomic<bool>* stopAll,
               double displayPeriodSec = 1.0);

    void StopJoin();

    /// Blocks until a full frame is available or shutdown. Copies nSym symbols into dstI/dstQ; false if stopped with no frame.
    bool WaitPopSymbolFrame(float* dstI, float* dstQ, int nSym);
    bool IsLocked() const { return locked_.load(std::memory_order_relaxed); }
    void ForceUnlock() { 
        locked_.store(false, std::memory_order_relaxed);
        gardner_.ForceUnlock();
    }
    double GetGardnerOmega() const { return gardner_.GetOmega(); }
    double GetGardnerOmegaNom() const { return gardner_.GetOmegaNom(); }

    void ThreadMain();

private:
    static constexpr int kSymFrame = ReceiverInputBatchIQSymbols;
    static constexpr int kPendingCap = 2 * kSymFrame;
    static constexpr int kFrameQueueDepth = 32;

    struct FrameSlot
    {
        alignas(32) float i[kSymFrame];
        alignas(32) float q[kSymFrame];
    };

    GardnerTiming gardner_;
    std::atomic<bool> locked_{false};
    float* oneSpsI_ = nullptr;
    float* oneSpsQ_ = nullptr;

    BufferFloat* pFilterBuf_ = nullptr;
    std::mutex* pMtxFilter_ = nullptr;
    std::condition_variable* pCvFilterData_ = nullptr;
    std::atomic<bool>* pStopAll_ = nullptr;

    alignas(32) float pendingI_[kPendingCap];
    alignas(32) float pendingQ_[kPendingCap];
    int pendingCount_ = 0;

    FrameSlot framePool_[kFrameQueueDepth];
    int frameQHead_ = 0;
    int frameQTail_ = 0;
    int frameQCount_ = 0;
    std::mutex mtxFrameQ_;
    std::condition_variable cvFrameReady_;
    std::condition_variable cvFrameSpace_;

    std::thread thread_;
    bool threadRunning_ = false;
    double displayPeriodSec_ = 1.0;

    bool pushOneFrameToQueue(const float* i, const float* q);
};
