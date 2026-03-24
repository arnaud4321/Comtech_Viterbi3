#include "SymbolRateEstimator.h"

#include <algorithm>
#include <cmath>
#include <fftw3.h>

namespace
{
constexpr double kPi = 3.14159265358979323846;
} // namespace

SymbolRateEstimator::SymbolRateEstimator()
{
    Reset(FsHz, Nfft, MaxOffsetHz);
}

SymbolRateEstimator::~SymbolRateEstimator()
{
    if (FftPlan)
        fftwf_destroy_plan(static_cast<fftwf_plan>(FftPlan));
    if (FftIn)
        fftwf_free(FftIn);
    if (FftOut)
        fftwf_free(FftOut);
}

void SymbolRateEstimator::Reset(double samplingFreqHz, int fftSize, int maxOffsetHz)
{
    FsHz = samplingFreqHz;
    Nfft = ClampPow2(fftSize);
    MaxOffsetHz = std::max(1000, maxOffsetHz);
    Fill = 0;
    BufferI.assign(Nfft, 0.0f);
    BufferQ.assign(Nfft, 0.0f);
    if (FftPlan)
    {
        fftwf_destroy_plan(static_cast<fftwf_plan>(FftPlan));
        FftPlan = nullptr;
    }
    if (FftIn)
    {
        fftwf_free(FftIn);
        FftIn = nullptr;
    }
    if (FftOut)
    {
        fftwf_free(FftOut);
        FftOut = nullptr;
    }
    FftIn = fftwf_alloc_complex(static_cast<size_t>(Nfft));
    FftOut = fftwf_alloc_complex(static_cast<size_t>(Nfft));
    FftPlan = fftwf_plan_dft_1d(
        Nfft,
        static_cast<fftwf_complex*>(FftIn),
        static_cast<fftwf_complex*>(FftOut),
        FFTW_FORWARD,
        FFTW_ESTIMATE);
}

int SymbolRateEstimator::GetRequiredBatches(int samplesPerBatch) const
{
    if (samplesPerBatch <= 0)
        return 1;
    return (Nfft + samplesPerBatch - 1) / samplesPerBatch;
}

bool SymbolRateEstimator::PushBatch(const float* I, const float* Q, int length)
{
    if (!I || !Q || length <= 0)
        return false;

    int copied = 0;
    while (copied < length && Fill < Nfft)
    {
        const int chunk = std::min(length - copied, Nfft - Fill);
        std::copy(I + copied, I + copied + chunk, BufferI.begin() + Fill);
        std::copy(Q + copied, Q + copied + chunk, BufferQ.begin() + Fill);
        Fill += chunk;
        copied += chunk;
    }
    return Fill >= Nfft;
}

SymbolRateEstimateResult SymbolRateEstimator::RunEstimation()
{
    SymbolRateEstimateResult out;
    out.FftResolutionHz = FsHz / static_cast<double>(Nfft);
    if (Fill < Nfft)
        return out;

    if (!FftPlan || !FftIn || !FftOut)
        return out;

    double meanPow = 0.0;
    for (int n = 0; n < Nfft; ++n)
    {
        const float ii = BufferI[n];
        const float qq = BufferQ[n];
        meanPow += static_cast<double>(ii) * ii + static_cast<double>(qq) * qq;
    }
    meanPow /= static_cast<double>(Nfft);

    for (int n = 0; n < Nfft; ++n)
    {
        const float ii = BufferI[n];
        const float qq = BufferQ[n];
        const float p = static_cast<float>(ii * ii + qq * qq - meanPow);
        const float w = 0.5f - 0.5f * std::cos(static_cast<float>(2.0 * kPi * n / (Nfft - 1)));
        const float shifted = ((n & 1) == 0) ? (p * w) : (-p * w); // shift Fs/2 to DC
        static_cast<fftwf_complex*>(FftIn)[n][0] = shifted;
        static_cast<fftwf_complex*>(FftIn)[n][1] = 0.0f;
    }

    fftwf_execute(static_cast<fftwf_plan>(FftPlan));

    const int kMax = std::min(Nfft / 2 - 2, static_cast<int>(std::round(MaxOffsetHz * Nfft / FsHz)));
    if (kMax < 2)
    {
        Fill = 0;
        return out;
    }

    std::vector<float> mags;
    mags.reserve(2 * kMax + 1);
    int bestK = 0;
    float bestM = -1.0f;
    for (int k = -kMax; k <= kMax; ++k)
    {
        const int idx = (k >= 0) ? k : (Nfft + k);
        const float re = static_cast<fftwf_complex*>(FftOut)[idx][0];
        const float im = static_cast<fftwf_complex*>(FftOut)[idx][1];
        const float m = re * re + im * im;
        mags.push_back(m);
        if (std::abs(k) >= 1 && m > bestM) // ignore exact DC
        {
            bestM = m;
            bestK = k;
        }
    }

    std::nth_element(mags.begin(), mags.begin() + mags.size() / 2, mags.end());
    const float medianM = std::max(1e-20f, mags[mags.size() / 2]);
    out.PeakToMedian = static_cast<double>(bestM / medianM);

    if (out.PeakToMedian < 8.0)
    {
        Fill = 0;
        return out;
    }

    const int leftK = bestK - 1;
    const int rightK = bestK + 1;
    const int idxL = (leftK >= 0) ? leftK : (Nfft + leftK);
    const int idxC = (bestK >= 0) ? bestK : (Nfft + bestK);
    const int idxR = (rightK >= 0) ? rightK : (Nfft + rightK);
    const float reL = static_cast<fftwf_complex*>(FftOut)[idxL][0];
    const float imL = static_cast<fftwf_complex*>(FftOut)[idxL][1];
    const float reC = static_cast<fftwf_complex*>(FftOut)[idxC][0];
    const float imC = static_cast<fftwf_complex*>(FftOut)[idxC][1];
    const float reR = static_cast<fftwf_complex*>(FftOut)[idxR][0];
    const float imR = static_cast<fftwf_complex*>(FftOut)[idxR][1];
    const float yL = reL * reL + imL * imL;
    const float yC = reC * reC + imC * imC;
    const float yR = reR * reR + imR * imR;

    double delta = 0.0;
    const double denom = static_cast<double>(yL) - 2.0 * static_cast<double>(yC) + static_cast<double>(yR);
    if (std::abs(denom) > 1e-20)
        delta = 0.5 * (static_cast<double>(yL) - static_cast<double>(yR)) / denom;
    delta = std::max(-0.5, std::min(0.5, delta));

    const double kInterp = static_cast<double>(bestK) + delta;
    const double fOffset = kInterp * FsHz / static_cast<double>(Nfft);
    out.SymbolRateHz = FsHz * 0.5 + fOffset;
    out.Detected = true;

    Fill = 0; // restart accumulation for next estimate window
    return out;
}

bool SymbolRateEstimator::IsPowerOfTwo(int n)
{
    return (n > 0) && ((n & (n - 1)) == 0);
}

int SymbolRateEstimator::ClampPow2(int n)
{
    if (n < 1024)
        return 1024;
    if (IsPowerOfTwo(n))
        return n;
    int p = 1;
    while (p < n)
        p <<= 1;
    return p;
}

