/**
 * @file SymbolRateEstimator.h
 * @brief Symbol-rate estimate from oversampled matched-filter IQ: |r[n]| → real FFT, line near Fs/2, parabolic peak.
 */
#pragma once

#include "IppComplexDft1d.h"

#include <vector>

//#define DEBUG_SYMBOL_RATE_ESTIMATOR_DUMP

/** @brief JSON-mapped parameters for @ref SymbolRateEstimator. */
struct SymbolRateEstimatorConfig
{
    int FftSize = 65536;
    int MaxOffsetHz = 200000;
    double EstimatePeriodSec = 1.0;
    double PeakToMedianThreshold = 8.0;
    double MaxRelativeJump = 0.02;
};

/** @brief Single estimation outcome (may be failure with @ref FailReason). */
struct SymbolRateEstimateResult
{
    bool Detected = false;
    double SymbolRateHz = 0.0;
    double PeakToMedian = 0.0;
    double FftResolutionHz = 0.0;
    int SearchKMaxBins = 0;        // half-width in bins: k in [k0-kMax, k0+kMax] (center k0 ≈ Fs/2 tone)
    double SearchMaxOffsetHz = 0.0; // kMax translated to Hz (half-width of search window around Fs/2)
    const char* FailReason = "";    // non-empty only when Detected == false
};

/**
 * @brief FFT-based symbol-rate estimator for oversampled matched-filter IQ (custom pipeline in this project).
 *
 * @details **Window** — @c PushBatch fills @c Nfft complex samples at base rate @f$F_s@f$.
 * **Optional upsampling** — If @c OsFactor @f$>1@f$ (default 2), Lagrange-4 interpolation in time yields @c NfftOs
 * samples at effective @f$f_{s,\mathrm{used}} = F_s \cdot OsFactor@f$ for the FFT stage.
 * **Nonlinearity** — Per sample @f$a_n = \sqrt{I_n^2+Q_n^2}@f$.
 * **FFT** — Forward **complex** DFT via Intel IPP (@c ippsDFTFwd_CToC_32fc); length @c Nfft or @c NfftOs (@c OsFactor path). In @c RunEstimation, length @f$N@f$
 * and rate @f$f_{s,\mathrm{used}}@f$ are the local @c nFftUsed / @c fsUsed. Input @f$z[n]=a_n+j0@f$. Output @f$X[k]@f$ is complex; @f$|X[k]|^2@f$ uses Re and Im.
 * **Reference bin @f$k_0@f$** — @f$k_0=\mathrm{round}\bigl((F_s/2)\cdot N/f_{s,\mathrm{used}}\bigr)@f$ with @f$F_s=@f$@c FsHz (member). So bin @f$k_0@f$ corresponds to **physical** @f$F_s/2@f$ Hz
 * on a grid sampled at @f$f_{s,\mathrm{used}}@f$ (either @f$F_s@f$ or @f$F_s\cdot@f$@c OsFactor).
 * **Search** — Half-width in bins @f$k_{\max}=\min(N/2-2,\,\mathrm{round}(\texttt{MaxOffsetHz}\cdot N/f_{s,\mathrm{used}}))@f$. For integers @f$m@f$ from @f$k_0-k_{\max}@f$
 * to @f$k_0+k_{\max}@f$, map @f$m@f$ to @f$[0,N)@f$ by wrap-around, score @f$|X|^2@f$, keep strongest bin whose wrapped index is **not** @f$k_0@f$.
 * **Refinement** — Parabolic interpolation on three @f$|X|^2@f$ samples; peak frequency → @ref SymbolRateEstimateResult::SymbolRateHz.
 * **Gating** — Peak/median ratio and relative change vs previous estimate are applied in @ref Receiver when consuming results.
 */
class SymbolRateEstimator
{
public:
    SymbolRateEstimator();
    ~SymbolRateEstimator() = default;

    /** @brief Reinitialize FFT sizes, thresholds, and internal buffers. */
    void Reset(double samplingFreqHz,
               int fftSize = 65536,
               int maxOffsetHz = 200000,
               double peakToMedianThreshold = 8.0);

    /** @brief Batches of @a samplesPerBatch complex samples needed before @ref RunEstimation. */
    int GetRequiredBatches(int samplesPerBatch) const;

    /** @brief Accumulate one processing batch; returns true when window is full. */
    bool PushBatch(const float* I, const float* Q, int length);

    /** @brief Execute FFT pipeline on the filled window (may set @c Detected false). */
    SymbolRateEstimateResult RunEstimation();

private:
    double FsHz = 21.42e6;
    int Nfft = 65536;
    int MaxOffsetHz = 200000;
    double PeakToMedianThreshold = 8.0;
    int Fill = 0;
    std::vector<float> BufferI;
    std::vector<float> BufferQ;
    IppComplexDft1d dft_;
    std::vector<Ipp32fc> fftIn_;
    std::vector<Ipp32fc> fftOut_;

    // Optional time-domain upsampling (Lagrange) before FFT when OsFactor > 1.
    double OsFactor = 2.0;
    int NfftOs = 0;
    std::vector<float> BufferOsI;
    std::vector<float> BufferOsQ;
    IppComplexDft1d dftOs_;
    std::vector<Ipp32fc> fftInOs_;
    std::vector<Ipp32fc> fftOutOs_;

    static bool IsPowerOfTwo(int n);
    static int ClampPow2(int n);
};

