/**
 * @file SymbolRateEstimator.cpp
 * @brief Implementation: I/Q window → optional Lagrange upsample → real-valued @f$|r[n]|@f$ FFT → peak band @f$\pm k_{\max}@f$ around @f$F_s/2@f$.
 *
 * @details **PushBatch** — Append batches until @c Fill == @c Nfft.
 * **RunEstimation** — Optionally interpolate to @c NfftOs (@c InterpLagrange4, @c OsFactor). Build FFT input
 * @f$\Re = |r[n]|@f$, @f$\Im = 0@f$. Intel IPP @c ippsDFTFwd_CToC_32fc; collect @f$|X[k]|^2@f$ for @f$k@f$ near @f$k_0@f$
 * (code: @c f0Hz = 0.5*FsHz mapped through @c fsUsed / @c nFftUsed). Peak bin ties toward @f$k_0@f$ (nominal line near @f$F_s/2@f$);
 * parabola on @f$y_L,y_C,y_R@f$ → @c SymbolRateHz. Reset @c Fill when done. Acceptance logic lives in @ref Receiver.
 */

#include "SymbolRateEstimator.h"
#include "Lagrange4Simd.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <sys/stat.h>

namespace
{
constexpr double kPi = 3.14159265358979323846;

inline int ClampIndex(int x, int lo, int hi)
{
    return std::max(lo, std::min(hi, x));
}

} // namespace

SymbolRateEstimator::SymbolRateEstimator()
{
    Reset(FsHz, Nfft, MaxOffsetHz);
}

void SymbolRateEstimator::Reset(double samplingFreqHz, int fftSize, int maxOffsetHz, double peakToMedianThreshold)
{
    FsHz = samplingFreqHz;
    Nfft = ClampPow2(fftSize);
    MaxOffsetHz = std::max(1000, maxOffsetHz);
    PeakToMedianThreshold = std::max(1.0, peakToMedianThreshold);
    Fill = 0;
    BufferI.assign(Nfft, 0.0f);
    BufferQ.assign(Nfft, 0.0f);

    dft_.ensureLength(Nfft);
    fftIn_.resize(static_cast<size_t>(Nfft));
    fftOut_.resize(static_cast<size_t>(Nfft));

    // Oversampled plan/buffers
    if (OsFactor > 1.0)
    {
        const int want = static_cast<int>(std::llround(OsFactor * static_cast<double>(Nfft)));
        // IPP DFT supports arbitrary lengths; keep exact oversampling ratio.
        NfftOs = std::max(1024, want);
        BufferOsI.assign(static_cast<size_t>(NfftOs), 0.0f);
        BufferOsQ.assign(static_cast<size_t>(NfftOs), 0.0f);
    }
    else
    {
        NfftOs = 0;
        BufferOsI.clear();
        BufferOsQ.clear();
    }

    if (NfftOs > 0)
    {
        dftOs_.ensureLength(NfftOs);
        fftInOs_.resize(static_cast<size_t>(NfftOs));
        fftOutOs_.resize(static_cast<size_t>(NfftOs));
    }
    else
    {
        dftOs_.ensureLength(0);
        fftInOs_.clear();
        fftOutOs_.clear();
    }
}

int SymbolRateEstimator::GetRequiredBatches(int samplesPerBatch) const
{
    if (samplesPerBatch <= 0)
        return 1;
    return (Nfft + samplesPerBatch - 1) / samplesPerBatch;
}

/**
 * @brief Append one batch toward the @c Nfft-sample FFT window; true when full.
 */
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

/**
 * @brief Execute the FFT-based symbol-rate estimate on the accumulated @c Nfft-sample window.
 */
SymbolRateEstimateResult SymbolRateEstimator::RunEstimation()
{
    SymbolRateEstimateResult out;
    const bool useOs = (OsFactor > 1.0 && NfftOs > 0 && dftOs_.valid());
    const int nFftUsed = useOs ? NfftOs : Nfft;
    const double fsUsed = useOs ? (FsHz * OsFactor) : FsHz;
    out.FftResolutionHz = fsUsed / static_cast<double>(nFftUsed);
    if (Fill < Nfft)
    {
        out.FailReason = "not enough samples (Fill < Nfft)";
        return out;
    }

    if (!dft_.valid() || static_cast<int>(fftIn_.size()) != Nfft || static_cast<int>(fftOut_.size()) != Nfft)
    {
        out.FailReason = "FFT plan/buffers not initialized";
        return out;
    }
    if (useOs && (static_cast<int>(fftInOs_.size()) != NfftOs || static_cast<int>(fftOutOs_.size()) != NfftOs))
    {
        out.FailReason = "FFT oversampled buffers not initialized";
        return out;
    }

    const float* inI = BufferI.data();
    const float* inQ = BufferQ.data();
    int inLen = Nfft;
    if (useOs)
    {
        const double scale = (static_cast<double>(Nfft) - 1.0) / (static_cast<double>(NfftOs) - 1.0);
        for (int n = 0; n < NfftOs; ++n)
        {
            const double t = static_cast<double>(n) * scale;
            int k = static_cast<int>(std::floor(t));
            k = ClampIndex(k, 1, Nfft - 3);
            const double mu = t - static_cast<double>(k);
            Lagrange4Simd::EvalIQContiguous_ps(BufferI.data() + k - 1, BufferQ.data() + k - 1, mu, &BufferOsI[n],
                                               &BufferOsQ[n]);
        }
        inI = BufferOsI.data();
        inQ = BufferOsQ.data();
        inLen = NfftOs;
    }

    Ipp32fc* fftInPtr = useOs ? fftInOs_.data() : fftIn_.data();
    Ipp32fc* fftOutPtr = useOs ? fftOutOs_.data() : fftOut_.data();
    for (int n = 0; n < inLen; ++n)
    {
        const float ii = inI[n];
        const float qq = inQ[n];
        const float a = std::sqrt(ii * ii + qq * qq);
        fftInPtr[n].re = a;
        fftInPtr[n].im = 0.0f;
    }

    if (useOs)
        dftOs_.forward(fftInOs_.data(), fftOutOs_.data());
    else
        dft_.forward(fftIn_.data(), fftOut_.data());

    const double kMaxReal =
        std::round(static_cast<double>(MaxOffsetHz) * static_cast<double>(nFftUsed) / fsUsed);
    const int kMax = std::min(nFftUsed / 2 - 2, static_cast<int>(kMaxReal));
    if (kMax < 2)
    {
        out.SearchKMaxBins = kMax;
        out.SearchMaxOffsetHz = static_cast<double>(kMax) * FsHz / static_cast<double>(Nfft);
        out.FailReason = "kMax < 2 (MaxOffsetHz too small for this Fs/Nfft)";
        Fill = 0;
        return out;
    }
    out.SearchKMaxBins = kMax;
    out.SearchMaxOffsetHz = static_cast<double>(kMax) * fsUsed / static_cast<double>(nFftUsed);

    std::vector<float> mags;
    mags.reserve(2 * kMax + 1);
    std::vector<int> idxs;
    idxs.reserve(2 * kMax + 1);
    const double f0Hz = 0.5 * FsHz;
    const int k0 = std::max(0, std::min(nFftUsed - 1, static_cast<int>(std::llround(f0Hz * nFftUsed / fsUsed))));
    for (int k = (k0 - kMax); k <= (k0 + kMax); ++k)
    {
        const int idx = (k >= 0 && k < nFftUsed) ? k : ((k % nFftUsed + nFftUsed) % nFftUsed);
        const float re = fftOutPtr[idx].re;
        const float im = fftOutPtr[idx].im;
        const float m = re * re + im * im;
        mags.push_back(m);
        idxs.push_back(idx);
    }

    float maxM = 0.0f;
    for (float m : mags)
        maxM = std::max(maxM, m);

    auto circBinDist = [nFftUsed](int a, int b) {
        const int d = std::abs(a - b);
        return std::min(d, nFftUsed - d);
    };

    // Among bins within kMagTieRatio of the band maximum, prefer the one closest to k0 (nominal baud line near Fs/2
    // at 2 Sps), then the strongest — reduces wrong Rs when a sidelobe is slightly higher than the true line.
    constexpr float kMagTieRatio = 0.92f;
    const float mFloor = maxM * kMagTieRatio;
    int bestK = k0;
    float bestM = -1.0f;
    int bestD = nFftUsed + 1;
    for (size_t i = 0; i < mags.size(); ++i)
    {
        const float mi = mags[i];
        const int idx = idxs[i];
        if (mi < mFloor)
            continue;
        const int d = circBinDist(idx, k0);
        if (d < bestD || (d == bestD && mi > bestM))
        {
            bestD = d;
            bestM = mi;
            bestK = idx;
        }
    }
    if (bestM < 0.0f)
    {
        for (size_t i = 0; i < mags.size(); ++i)
        {
            if (mags[i] > bestM)
            {
                bestM = mags[i];
                bestK = idxs[i];
            }
        }
    }

#ifdef DEBUG_SYMBOL_RATE_ESTIMATOR_DUMP
    {
        mkdir("../data", 0755);
        const int32_t len = static_cast<int32_t>(mags.size());
        const int32_t kMax32 = static_cast<int32_t>(kMax);
        const int32_t bestK32 = static_cast<int32_t>(bestK);
        const int32_t nfftUsed32 = static_cast<int32_t>(nFftUsed);
        const double fsBase = FsHz;
        const int32_t nfftBase32 = static_cast<int32_t>(Nfft);
        const double osFactor = useOs ? OsFactor : 1.0;
        const int32_t nTime32 = static_cast<int32_t>(inLen);
        FILE* f = std::fopen("../data/symbol_rate_mags.bin", "wb");
        if (f)
        {
            std::fwrite(&len, sizeof(len), 1, f);
            std::fwrite(&kMax32, sizeof(kMax32), 1, f);
            std::fwrite(&bestK32, sizeof(bestK32), 1, f);
            std::fwrite(&fsUsed, sizeof(fsUsed), 1, f);
            std::fwrite(&nfftUsed32, sizeof(nfftUsed32), 1, f);
            std::fwrite(&fsBase, sizeof(fsBase), 1, f);
            std::fwrite(&nfftBase32, sizeof(nfftBase32), 1, f);
            std::fwrite(&osFactor, sizeof(osFactor), 1, f);
            std::fwrite(&nTime32, sizeof(nTime32), 1, f);
            if (nTime32 > 0)
            {
                std::fwrite(inI, sizeof(float), static_cast<size_t>(nTime32), f);
                std::fwrite(inQ, sizeof(float), static_cast<size_t>(nTime32), f);
            }
            if (len > 0)
                std::fwrite(mags.data(), sizeof(float), static_cast<size_t>(len), f);

            std::vector<float> magsFull;
            magsFull.resize(static_cast<size_t>(nFftUsed));
            for (int k = 0; k < nFftUsed; ++k)
            {
                const float re = fftOutPtr[k].re;
                const float im = fftOutPtr[k].im;
                magsFull[static_cast<size_t>(k)] = re * re + im * im;
            }
            std::fwrite(magsFull.data(), sizeof(float), static_cast<size_t>(nFftUsed), f);
            std::fclose(f);
        }
    }
#endif

    std::nth_element(mags.begin(), mags.begin() + mags.size() / 2, mags.end());
    const float medianM = std::max(1e-20f, mags[mags.size() / 2]);
    out.PeakToMedian = static_cast<double>(bestM / medianM);

    if (out.PeakToMedian < PeakToMedianThreshold)
    {
        static thread_local std::string fail;
        fail = "peak/median below threshold (peak/median=" + std::to_string(out.PeakToMedian)
             + ", threshold=" + std::to_string(PeakToMedianThreshold) + ")";
        out.FailReason = fail.c_str();
        Fill = 0;
        return out;
    }

    const int leftK = bestK - 1;
    const int rightK = bestK + 1;
    const int idxL = (leftK >= 0) ? leftK : (nFftUsed + leftK);
    const int idxC = bestK;
    const int idxR = (rightK < nFftUsed) ? rightK : (rightK - nFftUsed);
    const float reL = fftOutPtr[idxL].re;
    const float imL = fftOutPtr[idxL].im;
    const float reC = fftOutPtr[idxC].re;
    const float imC = fftOutPtr[idxC].im;
    const float reR = fftOutPtr[idxR].re;
    const float imR = fftOutPtr[idxR].im;
    const float yL = reL * reL + imL * imL;
    const float yC = reC * reC + imC * imC;
    const float yR = reR * reR + imR * imR;

    double delta = 0.0;
    const double denom = static_cast<double>(yL) - 2.0 * static_cast<double>(yC) + static_cast<double>(yR);
    if (std::abs(denom) > 1e-20)
        delta = 0.5 * (static_cast<double>(yL) - static_cast<double>(yR)) / denom;
    delta = std::max(-0.5, std::min(0.5, delta));

    const double kInterp = static_cast<double>(bestK) + delta;
    const double fPeak = kInterp * fsUsed / static_cast<double>(nFftUsed);
    out.SymbolRateHz = fPeak;
    out.Detected = true;

    Fill = 0;
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
