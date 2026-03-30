#pragma once

#include "BufferFloat.h"
#include "definitions.h"
#include "GardnerTiming.h"
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

/// Dedicated thread: reads the filtered 2 sps stream, runs timing recovery (Gardner),
/// and outputs 1 sps symbol frames (kSymFrame = ReceiverInputBatchIQSymbols) to the Viterbi manager.
class ReceiverTimingTracking
{
public:
    ReceiverTimingTracking();
    ~ReceiverTimingTracking();

    void ResetGardner(double nominalRatio, double kp, double ki, int lockAvg);

    void Start(BufferFloat* filterRing, std::mutex* mtxFilterRing,
               std::condition_variable* cvFilterData, bool* stopAll);

    void StopJoin();

    /// Blocks until a full frame is available or shutdown. Copies nSym symbols into dstI/dstQ; false if stopped with no frame.
    bool WaitPopSymbolFrame(float* dstI, float* dstQ, int nSym);
    bool IsLocked() const { return locked_.load(std::memory_order_relaxed); }

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
    bool* pStopAll_ = nullptr;

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

    bool pushOneFrameToQueue(const float* i, const float* q);
};
