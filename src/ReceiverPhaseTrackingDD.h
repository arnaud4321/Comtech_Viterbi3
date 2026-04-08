/**
 * @file ReceiverPhaseTrackingDD.h
 * @brief Decision-directed carrier phase tracking (2nd-order PLL) after Gardner.
 */
#pragma once

#include "ReceiverTimingTracking.h"
#include "definitions.h"
#include <atomic>
#include <cstdio>
#include <condition_variable>
#include <mutex>
#include <thread>

// Uncomment to enable debug dumps (input/output constellation + phase + phase error)
// Files are written to ../data/:
// - phase_dd_in.bin      (float32 interleaved I,Q per symbol, input symbols to DD loop)
// - phase_dd_out.bin     (float32 interleaved I,Q per symbol, output after phase correction)
// - phase_dd_phase.bin   (float32 phase estimate in rad, per symbol)
// - phase_dd_err.bin     (float32 DD phase error in rad, per symbol)
// - phase_dd_err_filt.bin (float32 filtered error magnitude used for lock: EMA(|err|), per symbol)
// - phase_dd_thresh.bin   (float32 lock threshold in rad, per symbol)
// #define DEBUG_PHASE_DD

/**
 * @brief Decision-directed carrier phase tracking (QPSK) after @ref ReceiverTimingTracking.
 *
 * @details Per symbol: derotate @f$(I,Q)@f$ with LUT → @f$(z_{rI},z_{rQ})@f$; QPSK slice @f$\pm1@f$.
 * @f$e=\mathrm{atan2}(z_{rQ}d_i-z_{rI}d_q,\, z_{rI}d_i+z_{rQ}d_q)@f$ (= @c dd_phase_error). Then @c freqRadPerSym_ plus @c ki_ times @f$e@f$;
 * @c phaseRad_ plus @c freqRadPerSym_ plus @c kp_ times @f$e@f$, wrapped to @f$(-\pi,\pi]@f$; @c outI/outQ store current derotated @f$z_r@f$.
 * Lock when EMA of @f$|e|@f$ stays below @c lockThresholdRad_ for @c lockCountRequired_ symbols. Delock can re-enable VIVA in @ref ReceiverFreqCorrector.
 */
class ReceiverPhaseTrackingDD
{
public:
    ReceiverPhaseTrackingDD();
    ~ReceiverPhaseTrackingDD();

    /// Start the phase-tracking thread.
    /// Input is pulled from timingTracking (WaitPopSymbolFrame).
    void Start(ReceiverTimingTracking* timingTracking, std::atomic<bool>* stopAll,
               double displayPeriodSec = 1.0);

    void StopJoin();

    /// Blocks until a full corrected frame is available or shutdown.
    /// Copies nSym symbols into dstI/dstQ; false if stopped with no frame.
    bool WaitPopSymbolFrame(float* dstI, float* dstQ, int nSym);

    bool IsLocked() const { return locked_.load(std::memory_order_relaxed); }
    double GetLastFreqEstHz() const { return lastFreqEstHz_.load(std::memory_order_relaxed); }
    /**
     * @brief Get the last computed Error Vector Magnitude (EVM) RMS.
     * 
     * Computed independently of the AGC scaling using the ratio of signal power and
     * the mean symbol magnitude, ensuring an accurate SNR estimation even at high SNR.
     * 
     * @return double EVM RMS value
     */
    double GetLastEvmRms() const { return lastEvmRms_.load(std::memory_order_relaxed); }

    void ThreadMain();

private:
    static constexpr int kSymFrame = ReceiverInputBatchIQSymbols;
    static constexpr int kFrameQueueDepth = 32;

    struct FrameSlot
    {
        alignas(32) float i[kSymFrame];
        alignas(32) float q[kSymFrame];
    };

    ReceiverTimingTracking* timingTracking_ = nullptr;
    std::atomic<bool>* stopAll_ = nullptr;

    // 2nd-order PLL state (radians, radians/sample)
    float phaseRad_ = 0.0f;
    float freqRadPerSym_ = 0.0f;

    // Loop gains (tune empirically)
    float kp_ = 0.02f;
    float ki_ = 1.0e-5f;

    // Lock detection (based on smoothed |error|)
    float errEma_ = 0.0f;
    float errEmaAlpha_ = 0.001f;
    float lockThresholdRad_ = 0.35f;
    int lockCount_ = 0;
    int lockCountRequired_ = 10 * kSymFrame; // ~10 frames worth of symbols
    std::atomic<bool> locked_{false};
    std::atomic<double> lastFreqEstHz_{0.0};
    std::atomic<double> lastEvmRms_{0.0};

    // Phase rotator LUT (cos/sin table)
    static constexpr int kLutSize = 16384;
    float cosLut_[kLutSize];
    float sinLut_[kLutSize];

    // Frame queue
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

    FILE* dbgIn_ = nullptr;
    FILE* dbgOut_ = nullptr;
    FILE* dbgPhase_ = nullptr;
    FILE* dbgErr_ = nullptr;
    FILE* dbgErrFilt_ = nullptr;
    FILE* dbgThresh_ = nullptr;
    void openDebugFilesIfNeeded();
    void closeDebugFiles();
};

