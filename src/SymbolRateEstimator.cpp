#include "SymbolRateEstimator.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fftw3.h>
#include <sys/stat.h>

namespace
{
constexpr double kPi = 3.14159265358979323846;

inline int ClampIndex(int x, int lo, int hi)
{
    return std::max(lo, std::min(hi, x));
}

inline float InterpLagrange4(const float* x, int len, double t)
{
    int k = static_cast<int>(std::floor(t));
    k = ClampIndex(k, 1, len - 3);
    const double mu = t - static_cast<double>(k);
    const float x0 = x[k - 1];
    const float x1 = x[k + 0];
    const float x2 = x[k + 1];
    const float x3 = x[k + 2];

    const double c0 = -mu * (mu - 1.0) * (mu - 2.0) / 6.0;
    const double c1 = (mu + 1.0) * (mu - 1.0) * (mu - 2.0) / 2.0;
    const double c2 = -(mu + 1.0) * mu * (mu - 2.0) / 2.0;
    const double c3 = (mu + 1.0) * mu * (mu - 1.0) / 6.0;
    return static_cast<float>(c0 * x0 + c1 * x1 + c2 * x2 + c3 * x3);
}
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
    if (FftPlanOs)
        fftwf_destroy_plan(static_cast<fftwf_plan>(FftPlanOs));
    if (FftInOs)
        fftwf_free(FftInOs);
    if (FftOutOs)
        fftwf_free(FftOutOs);
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

    // Oversampled plan/buffers
    if (OsFactor > 1.0)
    {
        const int want = static_cast<int>(std::llround(OsFactor * static_cast<double>(Nfft)));
        // FFTW supports non power-of-two sizes; keep exact oversampling ratio.
        NfftOs = std::max(1024, want);
        BufferOsI.assign(NfftOs, 0.0f);
        BufferOsQ.assign(NfftOs, 0.0f);
    }
    else
    {
        NfftOs = 0;
        BufferOsI.clear();
        BufferOsQ.clear();
    }
    if (FftPlanOs)
    {
        fftwf_destroy_plan(static_cast<fftwf_plan>(FftPlanOs));
        FftPlanOs = nullptr;
    }
    if (FftInOs)
    {
        fftwf_free(FftInOs);
        FftInOs = nullptr;
    }
    if (FftOutOs)
    {
        fftwf_free(FftOutOs);
        FftOutOs = nullptr;
    }
    if (NfftOs > 0)
    {
        FftInOs = fftwf_alloc_complex(static_cast<size_t>(NfftOs));
        FftOutOs = fftwf_alloc_complex(static_cast<size_t>(NfftOs));
        FftPlanOs = fftwf_plan_dft_1d(
            NfftOs,
            static_cast<fftwf_complex*>(FftInOs),
            static_cast<fftwf_complex*>(FftOutOs),
            FFTW_FORWARD,
            FFTW_ESTIMATE);
    }
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
    const bool useOs = (OsFactor > 1.0 && NfftOs > 0 && FftPlanOs && FftInOs && FftOutOs);
    const int nFftUsed = useOs ? NfftOs : Nfft;
    const double fsUsed = useOs ? (FsHz * OsFactor) : FsHz;
    out.FftResolutionHz = fsUsed / static_cast<double>(nFftUsed);
    if (Fill < Nfft)
    {
        out.FailReason = "not enough samples (Fill < Nfft)";
        return out;
    }

    if (!FftPlan || !FftIn || !FftOut)
    {
        out.FailReason = "FFT plan/buffers not initialized";
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
            BufferOsI[n] = InterpLagrange4(BufferI.data(), Nfft, t);
            BufferOsQ[n] = InterpLagrange4(BufferQ.data(), Nfft, t);
        }
        inI = BufferOsI.data();
        inQ = BufferOsQ.data();
        inLen = NfftOs;
    }

    void* fftInPtr = useOs ? FftInOs : FftIn;
    for (int n = 0; n < inLen; ++n)
    {
        const float ii = inI[n];
        const float qq = inQ[n];
        const float a = std::sqrt(ii * ii + qq * qq);
        // No frequency shift: keep tone around Fs/2 in-place (we will search around Fs/2).
        static_cast<fftwf_complex*>(fftInPtr)[n][0] = a;
        static_cast<fftwf_complex*>(fftInPtr)[n][1] = 0.0f;
    }

    fftwf_execute(static_cast<fftwf_plan>(useOs ? FftPlanOs : FftPlan));

    // Cast to double before multiply/divide to avoid int32 overflow when MaxOffsetHz*Nfft is large.
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
    // We search around the *physical* center FsBase/2 (not fsUsed/2).
    // After oversampling, that tone moves to FsBase/2 Hz which corresponds to k0 = f0 * N / fsUsed.
    const double f0Hz = 0.5 * FsHz;
    const int k0 = std::max(0, std::min(nFftUsed - 1, static_cast<int>(std::llround(f0Hz * nFftUsed / fsUsed))));
    int bestK = k0;
    float bestM = -1.0f;
    void* fftOutPtr = useOs ? FftOutOs : FftOut;
    for (int k = (k0 - kMax); k <= (k0 + kMax); ++k)
    {
        const int idx = (k >= 0 && k < nFftUsed) ? k : ((k % nFftUsed + nFftUsed) % nFftUsed);
        const float re = static_cast<fftwf_complex*>(fftOutPtr)[idx][0];
        const float im = static_cast<fftwf_complex*>(fftOutPtr)[idx][1];
        const float m = re * re + im * im;
        mags.push_back(m);
        if (idx != k0 && m > bestM) // ignore exact center bin
        {
            bestM = m;
            bestK = idx;
        }
    }

#ifdef DEBUG_SYMBOL_RATE_ESTIMATOR_DUMP
    {
        mkdir("../data", 0755);
        // Dump mags BEFORE nth_element mutates it.
        // File format:
        //   [int32 len][int32 kMax][int32 bestK]
        //   [double FsHzUsed][int32 NfftUsed][double FsHzBase][int32 NfftBase][double OsFactor]
        //   [int32 nTime][float I_time[nTime]][float Q_time[nTime]]
        //   [float mags_zoom[len]]  (k = (center bin - kMax) : (center bin + kMax))
        //   [float mags_full[NfftUsed]] (k = 0:NfftUsed-1)
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
                const float re = static_cast<fftwf_complex*>(fftOutPtr)[k][0];
                const float im = static_cast<fftwf_complex*>(fftOutPtr)[k][1];
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

    if (out.PeakToMedian < 50.0)
    {
        out.FailReason = "peak/median below threshold";
        Fill = 0;
        return out;
    }

    const int leftK = bestK - 1;
    const int rightK = bestK + 1;
    const int idxL = (leftK >= 0) ? leftK : (nFftUsed + leftK);
    const int idxC = bestK;
    const int idxR = (rightK < nFftUsed) ? rightK : (rightK - nFftUsed);
    const float reL = static_cast<fftwf_complex*>(fftOutPtr)[idxL][0];
    const float imL = static_cast<fftwf_complex*>(fftOutPtr)[idxL][1];
    const float reC = static_cast<fftwf_complex*>(fftOutPtr)[idxC][0];
    const float imC = static_cast<fftwf_complex*>(fftOutPtr)[idxC][1];
    const float reR = static_cast<fftwf_complex*>(fftOutPtr)[idxR][0];
    const float imR = static_cast<fftwf_complex*>(fftOutPtr)[idxR][1];
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

