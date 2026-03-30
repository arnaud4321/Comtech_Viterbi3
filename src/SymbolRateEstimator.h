#pragma once

#include <vector>

//#define DEBUG_SYMBOL_RATE_ESTIMATOR_DUMP

struct SymbolRateEstimatorConfig
{
    int FftSize = 65536;
    int MaxOffsetHz = 200000;
    double EstimatePeriodSec = 1.0;
    double PeakToMedianThreshold = 8.0;
    double MaxRelativeJump = 0.02;
};

struct SymbolRateEstimateResult
{
    bool Detected = false;
    double SymbolRateHz = 0.0;
    double PeakToMedian = 0.0;
    double FftResolutionHz = 0.0;
    int SearchKMaxBins = 0;        // search range is k in [-SearchKMaxBins, +SearchKMaxBins]
    double SearchMaxOffsetHz = 0.0; // effective max |fOffset| searched around DC
    const char* FailReason = "";    // non-empty only when Detected == false
};

/*
 * Jablon-like symbol-rate estimator for 2x sampled matched-filter output.
 * It forms a non-linear sequence from |I+jQ|^2, shifts Fs/2 to DC via (-1)^n,
 * then finds the dominant FFT peak around DC (parabolic interpolation).
 */
class SymbolRateEstimator
{
public:
    SymbolRateEstimator();
    ~SymbolRateEstimator();

    void Reset(double samplingFreqHz, int fftSize = 65536, int maxOffsetHz = 200000);
    int GetRequiredBatches(int samplesPerBatch) const;
    bool PushBatch(const float* I, const float* Q, int length);
    SymbolRateEstimateResult RunEstimation();

private:
    double FsHz = 21.42e6;
    int Nfft = 65536;
    int MaxOffsetHz = 200000;
    int Fill = 0;
    std::vector<float> BufferI;
    std::vector<float> BufferQ;
    void* FftIn = nullptr;
    void* FftOut = nullptr;
    void* FftPlan = nullptr;

    // Jablon input oversampling via Lagrange interpolation (time domain).
    double OsFactor = 2.0;
    int NfftOs = 0;
    std::vector<float> BufferOsI;
    std::vector<float> BufferOsQ;
    void* FftInOs = nullptr;
    void* FftOutOs = nullptr;
    void* FftPlanOs = nullptr;

    static bool IsPowerOfTwo(int n);
    static int ClampPow2(int n);
};

