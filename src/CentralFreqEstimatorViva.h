/**
 * @file CentralFreqEstimatorViva.h
 * @brief QPSK coarse carrier estimator (fourth-power / VIVA sequence + FFT peak + parabolic fit).
 */
#pragma once

#include <cstddef>

/** @brief Tunables for @ref CentralFreqEstimatorViva. */
struct CentralFreqVivaConfig
{
    double MaxOffsetHz = 200000.0;          ///< Search \f$|f| \le\f$ MaxOffsetHz
    double PeakToMedianThreshold = 8.0;     ///< Accept if peak/median magnitude ratio ≥ threshold
    int FftSizeMultiplier = 4;              ///< @c Nfft = mult * nIn (zero-pad)
};

/** @brief Output of one VIVA FFT estimation window. */
struct CentralFreqVivaResult
{
    bool Detected = false;
    double OffsetHz = 0.0;          ///< Estimated carrier offset (Hz) at RF/baseband convention used by caller
    double PeakToMedian = 0.0;
    double AvgPower = 0.0;          ///< Mean \f$|r|^2\f$ on the input block
    double FftResolutionHz = 0.0;
};

/**
 * @brief Coarse carrier offset (Hz) for QPSK via fourth-power line and FFT (VIVA-equivalent nonlinearity).
 *
 * @details **Nonlinearity** — For each sample @f$r=I+jQ@f$, FFT input is
 * @f$v = \frac{|r|^2}{\overline{|r|^2}+\varepsilon}\,e^{j4\angle(r)}@f$, where \f$\overline{|r|^2}\f$ is the mean
 * of \f$|r|^2\f$ over the input block (reported as @ref CentralFreqVivaResult::AvgPower). This matches the
 * common VIVA/fourth-power construction and provides simple power normalization. Pad with zeros to length
 * @c FftSizeMultiplier × @c nIn.
 * **FFT search** — Forward DFT; compare @f$|Y[k]|^2@f$ for @f$k=1\ldots k_{\max}@f$ and symmetric index @f$N_{\mathrm{fft}}-k@f$,
 * excluding @f$k=0@f$. Bin limit @f$k_{\max}@f$ maps @f$4@f$ times @c MaxOffsetHz (Hz on true carrier) into the DFT grid.
 * **Acceptance** — peak/median ≥ @c PeakToMedianThreshold. **Refine** — parabola on three @f$|Y|^2@f$ samples in signed-bin space.
 * **Output** — interpolated line frequency @f$f_4@f$ on the fourth-power spectrum; @c OffsetHz is @f$f_4/4@f$.
 * Drives @ref ReceiverFreqCorrector.
 */
class CentralFreqEstimatorViva
{
public:
    CentralFreqEstimatorViva();
    ~CentralFreqEstimatorViva();

    /** @brief Set sample rate and estimator thresholds / FFT zero-padding. */
    void Configure(double samplingFreqHz, const CentralFreqVivaConfig& cfg);

    /**
     * @brief Run one FFT-based estimate over @a nIn complex samples.
     * @param inI,inQ Interleaved float IQ (same length).
     * @param nIn Input length (typically thousands of samples; minimum enforced in implementation).
     */
    CentralFreqVivaResult Run(const float* inI, const float* inQ, int nIn);

private:
    double fsHz_ = 0.0;
    CentralFreqVivaConfig cfg_{};

    void* fftIn_ = nullptr;   // fftwf_complex*
    void* fftOut_ = nullptr;  // fftwf_complex*
    void* fftPlan_ = nullptr; // fftwf_plan
    int nFft_ = 0;

    void ensureFftSize(int nFft);
};

