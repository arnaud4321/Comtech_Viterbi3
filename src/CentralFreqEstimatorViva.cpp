/**
 * @file CentralFreqEstimatorViva.cpp
 * @brief VIVA nonlinearity @f$v = |r|^2 e^{j4\\angle(r)}@f$ (optionally normalized), zero-padded FFT (FFTW float),
 * peak + parabolic refine → Hz/4.
 *
 * @details **Run** — First pass computes @c AvgPower = mean @f$|r|^2@f$ over @c nIn. Then for
 * @f$n<n_{\mathrm{in}}@f$: @f$X[n] \\leftarrow \\frac{|r[n]|^2}{\\mathrm{AvgPower}+\\varepsilon}\\,e^{j4\\angle(r[n])}@f$;
 * else zero. @c fftwf_execute; @f$k_{\max}@f$ from @f$4\\cdot@f$@c MaxOffsetHz. Loop @f$k=1\\ldots k_{\max}@f$:
 * compare @f$|Y[k]|^2@f$ and @f$|Y[N-k]|^2@f$; median of collected mags; parabolic refine;
 * @c OffsetHz = @f$f_{\mathrm{interp}}/4@f$.
 */

#include "CentralFreqEstimatorViva.h"

#include <algorithm>
#include <cmath>
#include <fftw3.h>
#include <vector>

CentralFreqEstimatorViva::CentralFreqEstimatorViva() = default;

CentralFreqEstimatorViva::~CentralFreqEstimatorViva()
{
    if (fftPlan_)
        fftwf_destroy_plan(static_cast<fftwf_plan>(fftPlan_));
    if (fftIn_)
        fftwf_free(fftIn_);
    if (fftOut_)
        fftwf_free(fftOut_);
}

void CentralFreqEstimatorViva::Configure(double samplingFreqHz, const CentralFreqVivaConfig& cfg)
{
    fsHz_ = samplingFreqHz;
    cfg_ = cfg;
    if (cfg_.FftSizeMultiplier < 1)
        cfg_.FftSizeMultiplier = 1;
}

void CentralFreqEstimatorViva::ensureFftSize(int nFft)
{
    if (nFft_ == nFft && fftPlan_ && fftIn_ && fftOut_)
        return;

    if (fftPlan_)
    {
        fftwf_destroy_plan(static_cast<fftwf_plan>(fftPlan_));
        fftPlan_ = nullptr;
    }
    if (fftIn_)
    {
        fftwf_free(fftIn_);
        fftIn_ = nullptr;
    }
    if (fftOut_)
    {
        fftwf_free(fftOut_);
        fftOut_ = nullptr;
    }

    nFft_ = std::max(1024, nFft);
    fftIn_ = fftwf_alloc_complex(static_cast<size_t>(nFft_));
    fftOut_ = fftwf_alloc_complex(static_cast<size_t>(nFft_));
    fftPlan_ = fftwf_plan_dft_1d(
        nFft_,
        static_cast<fftwf_complex*>(fftIn_),
        static_cast<fftwf_complex*>(fftOut_),
        FFTW_FORWARD,
        FFTW_ESTIMATE);
}

/**
 * @brief One-shot coarse carrier estimate over @a nIn complex samples.
 */
CentralFreqVivaResult CentralFreqEstimatorViva::Run(const float* inI, const float* inQ, int nIn)
{
    CentralFreqVivaResult out;
    if (!inI || !inQ || nIn < 1024 || fsHz_ <= 0.0)
        return out;

    const int nFft = cfg_.FftSizeMultiplier * nIn;
    ensureFftSize(nFft);
    if (!fftPlan_ || !fftIn_ || !fftOut_)
        return out;

    // Build VIVA sequence: v = |r|^2 * exp(j*4*angle(r)), normalized by AvgPower for stability.
    const float eps = 1e-12f;
    double sumP = 0.0;
    for (int n = 0; n < nIn; ++n)
    {
        const float i = inI[n];
        const float q = inQ[n];
        const float p = i * i + q * q;
        sumP += static_cast<double>(p);
    }

    out.AvgPower = sumP / static_cast<double>(nIn);
    const float invAvgP = 1.0f / static_cast<float>(out.AvgPower + static_cast<double>(eps));

    fftwf_complex* X = static_cast<fftwf_complex*>(fftIn_);
    for (int n = 0; n < nFft_; ++n)
    {
        float re = 0.0f, im = 0.0f;
        if (n < nIn)
        {
            const float i = inI[n];
            const float q = inQ[n];
            const float p = i * i + q * q;
            const float theta = std::atan2(q, i);
            const float phase = 4.0f * theta;
            const float amp = p * invAvgP;
            re = amp * std::cos(phase);
            im = amp * std::sin(phase);
        }
        X[n][0] = re;
        X[n][1] = im;
    }

    out.FftResolutionHz = fsHz_ / static_cast<double>(nFft_);

    fftwf_execute(static_cast<fftwf_plan>(fftPlan_));

    const double maxHz4 = 4.0 * std::max(0.0, cfg_.MaxOffsetHz);
    const int kMax = std::min(nFft_ / 2 - 2,
                              static_cast<int>(std::llround(maxHz4 * static_cast<double>(nFft_) / fsHz_)));
    if (kMax < 2)
        return out;

    // Search peak around DC (exclude k=0) in +/-kMax.
    float bestM = -1.0f;
    int bestK = 1;
    std::vector<float> mags;
    mags.reserve(static_cast<size_t>(2 * kMax));

    fftwf_complex* Y = static_cast<fftwf_complex*>(fftOut_);
    for (int k = 1; k <= kMax; ++k)
    {
        const float reP = Y[k][0];
        const float imP = Y[k][1];
        const float mP = reP * reP + imP * imP;
        mags.push_back(mP);
        if (mP > bestM)
        {
            bestM = mP;
            bestK = k;
        }

        const int kn = nFft_ - k;
        const float reN = Y[kn][0];
        const float imN = Y[kn][1];
        const float mN = reN * reN + imN * imN;
        mags.push_back(mN);
        if (mN > bestM)
        {
            bestM = mN;
            bestK = -k;
        }
    }

    if (mags.empty())
        return out;

    // Median magnitude for peak/median ratio
    const size_t mid = mags.size() / 2;
    std::nth_element(mags.begin(), mags.begin() + static_cast<long>(mid), mags.end());
    const float med = std::max(1e-20f, mags[mid]);
    out.PeakToMedian = static_cast<double>(bestM / med);

    if (out.PeakToMedian < cfg_.PeakToMedianThreshold)
        return out;

    // Parabolic refinement on |FFT|^2 at neighbors in signed-bin space (same idea as SymbolRateEstimator).
    auto fftIndexForSignedBin = [this](int sk) -> int {
        if (sk > 0)
            return sk;
        if (sk < 0)
            return nFft_ + sk;
        return 0;
    };
    const int skL = bestK - 1;
    const int skR = bestK + 1;
    const int idxL = fftIndexForSignedBin(skL);
    const int idxC = fftIndexForSignedBin(bestK);
    const int idxR = fftIndexForSignedBin(skR);

    const float reL = Y[idxL][0];
    const float imL = Y[idxL][1];
    const float reC = Y[idxC][0];
    const float imC = Y[idxC][1];
    const float reR = Y[idxR][0];
    const float imR = Y[idxR][1];
    const float yL = reL * reL + imL * imL;
    const float yC = reC * reC + imC * imC;
    const float yR = reR * reR + imR * imR;

    double delta = 0.0;
    const double denom = static_cast<double>(yL) - 2.0 * static_cast<double>(yC) + static_cast<double>(yR);
    if (std::abs(denom) > 1e-20)
        delta = 0.5 * (static_cast<double>(yL) - static_cast<double>(yR)) / denom;
    delta = std::max(-0.5, std::min(0.5, delta));

    const double signedInterp = static_cast<double>(bestK) + delta;
    const double freqHzPow4 = signedInterp * fsHz_ / static_cast<double>(nFft_);
    // The estimator runs on (I+jQ)^4 for QPSK, so measured frequency is 4x the true offset.
    out.OffsetHz = freqHzPow4 / 4.0;
    out.Detected = true;
    return out;
}

