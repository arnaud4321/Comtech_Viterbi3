#pragma once
// Uncomment to dump early/on-time/late (3 *.bin files): block-based writes, not per-symbol.
//#define DEBUG_GARDNER_OUTPUTS
// #define DEBUG_GARDNER_DEBIT

#include <vector>

class GardnerTiming 
{
public:
    GardnerTiming();
    ~GardnerTiming();

    void Reset(double omegaNom = 2.0, double kp = 1e-3, double ki = 1e-5, int updatePeriod = 64);
    int ProcessBlock(const float* inI, const float* inQ, int inLen,
                     float* outI, float* outQ, int outMax);
    bool IsLocked() const { return locked_; }

private:
    void PushToRing(const float* inI, const float* inQ, int inLen);
    long long OldestAbs() const;
    double NewestSafe() const;
    bool CanInterp(double t) const;
    float GetRingI(long long absIdx) const;
    float GetRingQ(long long absIdx) const;
    float InterpLagrange4I(double t) const;
    float InterpLagrange4Q(double t) const;
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

