#pragma once
// Uncomment to dump early/on-time/late (3 *.bin files): block-based writes, not per-symbol.
#define DEBUG_GARDNER_OUTPUTS
// #define DEBUG_GARDNER_DEBIT

class GardnerTiming
{
public:
    GardnerTiming();
    ~GardnerTiming();

    void Reset(double omegaNom = 2.0, double kp = 1e-3, double ki = 1e-5, int updatePeriod = 64);
    int ProcessBlock(const float* inI, const float* inQ, int inLen,
                     float* outI, float* outQ, int outMax);

private:
    static inline float InterpLagrange4(const float* x, int len, double t);
    static inline int ClampIndex(int idx, int lo, int hi);

    double omegaNom_ = 2.0;
    double omega_ = 2.0;
    double kp_ = 1e-3;
    double ki_ = 1e-6;
    double integ_ = 0.0;
    int updatePeriod_ = 64;
    int updateCounter_ = 0;
    double tCursor_ = 0.0; // consistent with Reset() and TakeEven grid (0,2,4,…)
    static constexpr int kOverlap = 4;
    float overlapI_[kOverlap] = {0.0f, 0.0f, 0.0f, 0.0f};
    float overlapQ_[kOverlap] = {0.0f, 0.0f, 0.0f, 0.0f};
    bool hasOverlap_ = false;

#ifdef DEBUG_GARDNER_OUTPUTS
    void OpenDebugFilesIfNeeded();
    void CloseDebugFiles();
    void* fidEarly_ = nullptr;
    void* fidOnTime_ = nullptr;
    void* fidLate_ = nullptr;
    void* fidErr_ = nullptr;
    void* fidOmega_ = nullptr;
#endif
};

