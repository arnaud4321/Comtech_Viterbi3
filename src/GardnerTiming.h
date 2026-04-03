/**
 * @file GardnerTiming.h
 * @brief Gardner TED timing recovery with Lagrange interpolation at 2 samples/symbol.
 */
#pragma once
// Uncomment to dump early/on-time/late (3 *.bin files): block-based writes, not per-symbol.
//#define DEBUG_GARDNER_OUTPUTS
// #define DEBUG_GARDNER_DEBIT

#include <vector>

/**
 * @brief Interpolating Gardner detector producing one sample per symbol (QPSK-oriented).
 *
 * Maintains a large ring buffer, fractional-time Lagrange-4 interpolation (@ref Lagrange4Simd, AVX for I/Q pairs), and PI loop on \f$\omega\f$
 * (samples per symbol in the 2 sps domain). Lock is inferred from stability of \f$\omega\f$.
 *
 * **Update cadence for @f$\omega@f$:** @f$\omega@f$ is **constant** between updates. It is recomputed **once every
 * @c updatePeriod_ output symbols** (default @c 64 from @c Reset): i.e. after @c updatePeriod_ TED accumulations in
 * @c errorAcc_. Approximate time between two @f$\omega@f$ updates: @f$\mathrm{updatePeriod\_}\cdot\omega/F_s@f$ seconds
 * (@f$F_s@f$ = sample rate of the 2 sps branch, @f$\omega@f$ = samples per symbol at update time).
 *
 * @details **Per output symbol:** Lagrange-4 samples at @f$t-\omega/2@f$, @f$t@f$, @f$t+\omega/2@f$; accumulate
 * @f$(I_{\mathrm{late}}-I_{\mathrm{early}})I_{\mathrm{mid}}+(Q_{\mathrm{late}}-Q_{\mathrm{early}})Q_{\mathrm{mid}}@f$
 * into @c errorAcc_. On every @c updatePeriod_-th output symbol: mean = @c errorAcc_/updatePeriod_; @c integ_ += @c ki_ * mean;
 * @f$\omega@f$ := @c omegaNom_ + @c integ_ + @c kp_ * mean, then clamp @f$\omega\in[1.7,2.3]@f$, @c errorAcc_=0, @c PushOmegaHistoryOnUpdate.
 * Between PI ticks @f$\omega@f$ stays constant.
 * Advance @f$t@f$ by @f$\omega@f$ after each symbol output.
 */
class GardnerTiming 
{
public:
    GardnerTiming();
    ~GardnerTiming();

    /**
     * @brief Reset loop filter, ring, and lock state.
     * @param omegaNom Nominal @f$\omega@f$ (samples/symbol at 2 sps), typically 2.
     * @param kp Proportional gain on mean TED.
     * @param ki Integral gain on mean TED.
     * @param updatePeriod Number of **output symbols** between each @f$\omega@f$ / PI update (@f$\omega@f$ fixed in between).
     */
    void Reset(double omegaNom = 2.0, double kp = 1e-3, double ki = 1e-5, int updatePeriod = 64);

    /**
     * @brief Consume @a inLen input samples; write up to @a outMax one-sps outputs.
     * @return Number of complex symbols written to @a outI / @a outQ.
     */
    int ProcessBlock(const float* inI, const float* inQ, int inLen,
                     float* outI, float* outQ, int outMax);
    bool IsLocked() const { return locked_; }
    /// Last omega span used for lock decision: |max(historyOmega)-min(historyOmega)|.
    double GetLastOmegaLockSpan() const { return lastOmegaLockSpan_; }
    /// Omega span threshold for lock decision.
    double GetOmegaLockSpanThreshold() const { return kOmegaLockSpanThreshold; }
    /// Current omega (samples per symbol in the 2sps domain). Nominal is omegaNom.
    double GetOmega() const { return omega_; }
    /// Nominal omega (samples per symbol).
    double GetOmegaNom() const { return omegaNom_; }

private:
    void PushToRing(const float* inI, const float* inQ, int inLen);
    long long OldestAbs() const;
    double NewestSafe() const;
    bool CanInterp(double t) const;
    float GetRingI(long long absIdx) const;
    float GetRingQ(long long absIdx) const;
    float InterpLagrange4I(double t) const;
    float InterpLagrange4Q(double t) const;
    void InterpLagrange4IQ(double t, float* outI, float* outQ) const;
    void PushOmegaHistoryOnUpdate();

    double omegaNom_ = 2.0;
    double omega_ = 0.0; // set in Reset()
    double kp_ = 0.0; // set in Reset()
    double ki_ = 0.0; // set in Reset()
    double integ_ = 0.0;
    int updatePeriod_ = 64;
    int updateCounter_ = 0;
    // Input ring buffer (absolute time cursor, no clamping at block boundaries).
    static constexpr int kRingSize = 1 << 18; // power of 2
    static constexpr int kRingMask = kRingSize - 1;
    std::vector<float> ringI_;
    std::vector<float> ringQ_;
    long long absWrite_ = 0; // absolute sample index of next write
    double tAbs_ = 0.0;      // absolute time cursor in input-sample units
    bool tInit_ = false;
    double errorAcc_ = 0.0;
    // Lock detection: omega stability over a time-based window, using decimated omega updates.
    // Target: keep ~kOmegaLockMeasures omega points spanning ~kOmegaLockWindowSeconds.
    static constexpr double kOmegaLockWindowSeconds = 0.1;
    static constexpr int kOmegaLockMeasures = 100;
    double kOmegaLockSpanThreshold ; // lock if (max-min) <= threshold
    // value depends on ki and kp !!
    std::vector<double> historyOmega_;
    int historyOmegaWr_ = 0;
    int historyOmegaFill_ = 0;
    int omegaLockStrideUpdates_ = 1;  // push 1 point every N omega recomputations
    int omegaLockStrideCounter_ = 0;
    bool locked_ = false;
    double lastOmegaLockSpan_ = 0.0;

#ifdef DEBUG_GARDNER_OUTPUTS
    void OpenDebugFilesIfNeeded();
    void CloseDebugFiles();
    void* fidEarly_ = nullptr;
    void* fidOnTime_ = nullptr;
    void* fidLate_ = nullptr;
    void* fidErr_ = nullptr;
    void* fidOmega_ = nullptr;
    void* fidPhi_ = nullptr;
    // Interleaved float32 I,Q: 8 complex samples per symbol (evenly spaced over omega_) for eye plots.
    void* fidOut8Sps_ = nullptr;
#endif
};

