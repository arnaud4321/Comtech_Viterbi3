#include "CentralFreqEstimator.h"
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <sys/stat.h>
#include <fftw3.h>

/* Path from run cwd (typically src/): project root data/outfft.bin */
static const char* OutFftPath = "../data/outfft.bin";

static const double SearchRangeHz = 1e6;  /* peak search ±1 MHz */

CentralFreqEstimator::CentralFreqEstimator() = default;

CentralFreqEstimator::~CentralFreqEstimator()
{
    if (fftwPlan)  fftwf_destroy_plan(static_cast<fftwf_plan>(fftwPlan));
    if (fftwIn)    fftwf_free(fftwIn);
    if (fftwOut)   fftwf_free(fftwOut);
}

void CentralFreqEstimator::RunEstimation(const float* I, const float* Q, int Length)
{
    if (!I || !Q || Length < 2 || SamplingFrequency <= 0)
        return;

    int N = Length;
    if (fftwN != N)
    {
        if (fftwPlan)  fftwf_destroy_plan(static_cast<fftwf_plan>(fftwPlan));
        if (fftwIn)    fftwf_free(fftwIn);
        if (fftwOut)   fftwf_free(fftwOut);
        fftwN = N;
        fftwIn  = fftwf_alloc_complex(static_cast<size_t>(N));
        fftwOut = fftwf_alloc_complex(static_cast<size_t>(N));
        fftwPlan = fftwf_plan_dft_1d(N,
            static_cast<fftwf_complex*>(fftwIn),
            static_cast<fftwf_complex*>(fftwOut),
            FFTW_FORWARD, FFTW_ESTIMATE);
    }

    for (int i = 0; i < N; i++)
    {
        float ii = I[i], qq = Q[i];
        float i2 = ii * ii, q2 = qq * qq, iq = ii * qq;
        float re_z4 = (i2 - q2) * (i2 - q2) - 4.f * i2 * q2;
        float im_z4 = 4.f * iq * (i2 - q2);
        static_cast<fftwf_complex*>(fftwIn)[i][0] = re_z4;
        static_cast<fftwf_complex*>(fftwIn)[i][1] = im_z4;
    }

    fftwf_execute(static_cast<fftwf_plan>(fftwPlan));

#ifdef DEBUG_ESTIMATOR
    if (!DebugFirstEstimationWritten)
    {
        mkdir("../data", 0755);  /* ignore error if exists */
        FILE* f = fopen(OutFftPath, "wb");
        if (f)
        {
            /* FFT output: N complex values as I,Q float pairs (real, imag per bin) */
            fwrite(fftwOut, sizeof(float), 2 * static_cast<size_t>(N), f);
            fclose(f);
            DebugFirstEstimationWritten = true;
        }
    }
#endif

    double binWidth = SamplingFrequency / N;
    int kMax = static_cast<int>(4.0 * SearchRangeHz * N / SamplingFrequency);
    int kMaxBins = std::min(N / 2, std::max(1, kMax));

    /* Signal is from frequency transposition; no reason to exclude DC. */
    float maxMag2 = 0.f;
    int peakBin = 0;
    for (int k = 0; k <= kMaxBins; k++)
    {
        float re = static_cast<fftwf_complex*>(fftwOut)[k][0];
        float im = static_cast<fftwf_complex*>(fftwOut)[k][1];
        float m2 = re * re + im * im;
        if (m2 > maxMag2) { maxMag2 = m2; peakBin = k; }
    }
    for (int k = std::max(N - kMaxBins, N/2 + 1); k < N; k++)
    {
        float re = static_cast<fftwf_complex*>(fftwOut)[k][0];
        float im = static_cast<fftwf_complex*>(fftwOut)[k][1];
        float m2 = re * re + im * im;
        if (m2 > maxMag2) { maxMag2 = m2; peakBin = k; }
    }

    double freqPeak = (peakBin <= N / 2) ? peakBin * binWidth : (peakBin - N) * binWidth;
    double rawHz = freqPeak / 4.0;

    double maxOffset = SamplingFrequency / 8.0;
    if (rawHz > maxOffset)  rawHz = maxOffset;
    if (rawHz < -maxOffset) rawHz = -maxOffset;

    EstimatedOffsetHz = rawHz;
}
