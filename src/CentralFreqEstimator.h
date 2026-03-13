#pragma once

#define DEBUG_ESTIMATOR  /* when defined: write first FFT output to outfft.bin (I,Q floats); use #ifdef in code */

/** Estimates central frequency (QPSK 4th-power + FFT). No correction applied; returns coherent value in Hz. */
class CentralFreqEstimator
{
public:
    CentralFreqEstimator();
    ~CentralFreqEstimator();

    void SetSamplingFrequency(double fs_Hz) { SamplingFrequency = fs_Hz; }
    void SetInitialOffsetHz(double hz) { EstimatedOffsetHz = hz; }
    void RunEstimation(const float* I, const float* Q, int Length);
    double GetEstimatedOffsetHz() const { return EstimatedOffsetHz; }

private:
    double SamplingFrequency = 21.42e6;
    double EstimatedOffsetHz = 0.0;

    void* fftwPlan = nullptr;
    void* fftwIn   = nullptr;
    void* fftwOut  = nullptr;
    int fftwN = 0;
    bool DebugFirstEstimationWritten = false;
};
