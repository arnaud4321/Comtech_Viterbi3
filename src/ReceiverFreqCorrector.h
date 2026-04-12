/**
 * @file ReceiverFreqCorrector.h
 * @brief Thread: AGC + NCO driven by optional @ref CentralFreqEstimatorViva (VIVA) between resampler and Gardner.
 */
#pragma once

#include "BufferFloat.h"
#include "CentralFreqEstimatorViva.h"
#include "FrequencyOffset.h"
#include "ReceiverPhaseTrackingDD.h"
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <thread>

/** @brief Configuration for NCO, estimator cadence, and AGC. */
struct ReceiverFreqCorrectorConfig
{
    bool Enable = true;
    int EstimationBlockSamples = 65536; // complex samples at 2 sps
    CentralFreqVivaConfig Estimator{};
    double EstimatePeriodSec = 1.0;     // re-estimate periodically (s)
    double NcoHzEmaAlpha = 0.2;         // smooth NCO updates: ncoHz = (1-a)*ncoHz + a*rawHz
    // Slew-rate limiter in Hz/s (limits NCO changes smoothly vs dynamics).
    // Set <=0 to disable limiting.
    double MaxHzSlewRate = 8000.0;
    double PowerEmaAlpha = 0.01;
    // If > 0, normalize to this absolute average power (instead of "first observed").
    double TargetAvgPower = -1.0;
};

/**
 * @brief Dedicated thread pulling 2 sps samples, applying NCO + AGC, forwarding to timing recovery.
 *
 * @details **Coarse frequency:** when enabled, accumulates complex samples until
 * @c EstimationBlockSamples, runs @ref CentralFreqEstimatorViva::Run, then sets @c targetHz_ to the
 * raw offset; the NCO then tracks @c targetHz_ with @c NcoHzEmaAlpha and
 * @c MaxHzSlewRate before calling @ref FrequencyOffset::SetFrequency on
 * the NCO. Until a first valid estimate is available, forwarding may be gated (@c freqReady_).
 *
 * **AGC:** each chunk uses the **sample mean** of @f$I^2+Q^2@f$ (mean instantaneous power over complex
 * samples). That value is low-pass filtered with @c PowerEmaAlpha into @c pwrEma_. The reference power
 * @c pwrRef_ is @c TargetAvgPower when @f$>0@f$; otherwise it is seeded from the first chunk and, if no
 * absolute target is set, slowly tracks @c CentralFreqVivaResult::AvgPower from each VIVA run. I and Q are scaled by the **linear** gain
 * @f$g=\sqrt{P_{\mathrm{ref}}/P_{\mathrm{ema}}}@f$; @c gainDb_ is @f$20\log_{10} g@f$ for logging only.
 *
 * **Re-acquisition:** if @ref ReceiverPhaseTrackingDD loses lock, estimation can be re-enabled to steer
 * the NCO after large slips.
 */
class ReceiverFreqCorrector
{
public:
    ReceiverFreqCorrector();
    ~ReceiverFreqCorrector();

    void Configure(const ReceiverFreqCorrectorConfig& cfg, double displayPeriodSec);

    void Start(BufferFloat* inRing, std::mutex* inMtx, std::condition_variable* inCvData,
               BufferFloat* outRing, std::mutex* outMtx, std::condition_variable* outCvData,
               ReceiverPhaseTrackingDD* phaseDD, std::atomic<bool>* stopAll);

    void StopJoin();

    void ThreadMain();

    bool IsFreqReady() const { return freqReady_.load(std::memory_order_relaxed); }
    double GetNcoHz() const { return ncoHz_.load(std::memory_order_relaxed); }
    double GetTargetHz() const { return targetHz_.load(std::memory_order_relaxed); }
    double GetGainDb() const { return gainDb_.load(std::memory_order_relaxed); }

    // Force estimation to restart, used when Viterbi fails to lock for too long
    void ForceEstimationRestart() { forceRestart_.store(true, std::memory_order_relaxed); }

private:
    ReceiverFreqCorrectorConfig cfg_{};
    double displayPeriodSec_ = 1.0;
    
    std::atomic<bool> forceRestart_{false};

    BufferFloat* inRing_ = nullptr;
    std::mutex* inMtx_ = nullptr;
    std::condition_variable* inCvData_ = nullptr;

    BufferFloat* outRing_ = nullptr;
    std::mutex* outMtx_ = nullptr;
    std::condition_variable* outCvData_ = nullptr;

    ReceiverPhaseTrackingDD* phaseDD_ = nullptr;
    std::atomic<bool>* stopAll_ = nullptr;

    std::thread thread_;
    bool running_ = false;

    // NCO
    FrequencyOffset nco_;
    std::atomic<double> ncoHz_{0.0};
    std::atomic<double> targetHz_{0.0};
    bool hasTargetHz_ = false;

    // Gain normalization (keep power close to reference)
    double pwrEma_ = 0.0;
    double pwrRef_ = 0.0;
    std::atomic<double> gainDb_{0.0};

    // Estimator
    CentralFreqEstimatorViva est_;
    std::vector<float> accI_;
    std::vector<float> accQ_;
    int accCount_ = 0;
    bool allowEstimation_ = true;
    bool prevPhaseLocked_ = false;
    std::atomic<bool> freqReady_{false};

    uint64_t samplesTotal_ = 0;   // complex samples at 2 sps domain
    uint64_t droppedChunks_ = 0;  // chunks not forwarded because !freqReady_
};

