/**
 * @file ReceiverResampler.h
 * @brief Dedicated thread: rate conversion between matched filter output and ~2 sps symbol stream.
 */
#pragma once

#include "BufferFloat.h"
#include "Resampler.h"
#include "definitions.h"
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

/**
 * @brief Resamples matched-filter output toward nominal 2 samples/symbol using @ref Resampler (Lagrange).
 *
 * @details The matched filter runs at fixed simulation @c SamplingFrequency; true @f$R_s@f$ may differ. This
 * block does **closed-loop rate conversion**: @ref UpdateFromSymbolRateHz sets @c advance_ to
 * @f$\approx F_s/(2R_s)@f$, clamped to @f$[0.5,2]@f$ (see @ref ReceiverResampler::UpdateFromSymbolRateHz).
 *
 * **ThreadMain:** drains @c inRing_ in chunks of @c kReadChunk = @c SPB when at least @c SPB samples are
 * available and the internal FIFO has room below @c kFifoMax (else backpressure — filter ring can fill).
 * @ref Resampler::CreateOutputs uses atomic @c advance_, maintaining @c frac_ and consuming FIFO samples while
 * keeping @c kOverlap = 3 samples of lookahead; @c kOutBlockMax = @c 2·SPB caps one write burst. Output commits
 * to @c outRing_ when space exists (waits on @c outCvSpace_).
 */
class ReceiverResampler
{
public:
    ReceiverResampler();
    ~ReceiverResampler();

    /**
     * @brief Start worker thread; wires matched-filter output ring to the next RX stage (e.g. frequency corrector).
     * @param inRing Matched-filter @c BufferFloat (producer: filter thread).
     * @param inMtx Mutex protecting @a inRing (same as filter thread).
     * @param inCvData Signaled when @a inRing gains data (filter commits); worker waits/refills FIFO.
     * @param outRing Downstream @c BufferFloat (e.g. NCO/AGC input).
     * @param outMtx Mutex protecting @a outRing.
     * @param outCvData Signaled after writing to @a outRing.
     * @param outCvSpace Signaled when @a outRing has space (consumer advances read pointer).
     * @param stopAll Shared shutdown flag (worker exits when set).
     */
    void Start(BufferFloat* inRing, std::mutex* inMtx, std::condition_variable* inCvData,
               BufferFloat* outRing, std::mutex* outMtx, std::condition_variable* outCvData,
               std::condition_variable* outCvSpace, bool* stopAll);

    void StopJoin();

    /**
     * @brief Update interpolation ratio from a symbol-rate estimate (Hz).
     * @details Sets Advance \f$\approx F_s^{\mathrm{in}} / (2 R_s)\f$ (clamped in implementation).
     */
    void UpdateFromSymbolRateHz(double symbolRateHz);

    /** @brief Thread body: drain input ring → internal FIFO → @ref Resampler::CreateOutputs → output ring. */
    void ThreadMain();

private:
    static constexpr int kOverlap = 3; // cubic interpolation needs 4 points => keep last 3
    /// Max samples queued in the internal FIFO (fifoI_/fifoQ_ after fifoRd_). Backpressure: stop draining inRing_ above this.
    static constexpr int kFifoMax = SPB * 512;
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

