#pragma once

#include "BufferFloat.h"
#include "Resampler.h"
#include "definitions.h"
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

/// Dedicated thread that resamples the matched-filter output stream (2 sps nominal)
/// to correct symbol-rate mismatch. Uses ring buffers to handle variable I/O rates.
class ReceiverResampler
{
public:
    ReceiverResampler();
    ~ReceiverResampler();

    void Start(BufferFloat* inRing, std::mutex* inMtx, std::condition_variable* inCvData,
               BufferFloat* outRing, std::mutex* outMtx, std::condition_variable* outCvData,
               std::condition_variable* outCvSpace, bool* stopAll);

    void StopJoin();

    /// Update resampling ratio from a new symbol-rate estimate (Hz).
    /// Internally converts to Advance = Fs_in / Fs_out, with Fs_out = 2 * Rs_est.
    void UpdateFromSymbolRateHz(double symbolRateHz);

    void ThreadMain();

private:
    static constexpr int kOverlap = 3; // cubic interpolation needs 4 points => keep last 3
    BufferFloat* inRing_ = nullptr;
    std::mutex* inMtx_ = nullptr;
    std::condition_variable* inCvData_ = nullptr;

    BufferFloat* outRing_ = nullptr;
    std::mutex* outMtx_ = nullptr;
    std::condition_variable* outCvData_ = nullptr;
    std::condition_variable* outCvSpace_ = nullptr;

    bool* stopAll_ = nullptr;

    std::atomic<double> advance_{1.0};
    Resampler resampler_;
    double frac_ = 0.0;
    std::vector<float> fifoI_;
    std::vector<float> fifoQ_;
    size_t fifoRd_ = 0;

    std::thread thread_;
    bool running_ = false;
};

