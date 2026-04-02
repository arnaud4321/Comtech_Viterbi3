/**
 * @file AWGNChannel.h
 * @brief Multithreaded channel: AWGN, Es/N0 scaling, optional frequency profile, SCO, IQ to @c BufferShort.
 */
#pragma once
#include "random_generator_new.h"
class Transmitter;
#include <thread>
#include <mutex>
#include <condition_variable>
#include <immintrin.h>
#include "SimpleQueue.h"
#include "BufferShort.h"
#include "FrequencyOffset.h"
#include "ChannelSamplingClockOffset.h"
#include "Params.h"
#include <atomic>
using namespace std;

//#define DEBUG_AWGN

/**
 * @brief Stochastic channel model between @ref Transmitter and @ref Sampler.
 *
 * @details **Noise and scaling:** One thread fills Gaussian buffers; @ref AWGNChannel.cpp (@c GenerateOutput)
 * forms @c mkSig * TX + @c Stdn * noise from Es/N0 in @ref Params, a fixed reference @c Pow0, TX power, and
 * @c Backoff so the scaled signal+noise fits int16 headroom after optional dynamic range and SCO.
 *
 * **Frequency and level profiles:** Four simulation-time segments (stable / ramp / stable / ramp) set both
 * carrier offset (@c InitialFrequencyShift + ramped @c FrequencyShift) and **gain** in dB
 * (@c InitialGainDb + ramped @c DynamicRangeDb), applied as amplitude @f$10^{\mathrm{gainDb}/20}@f$ on floats
 * before the NCO. @ref FrequencyOffset applies rotation per chunk; total SCO ppm is @c ClockMismatchPpm plus a
 * Doppler-like term from offset Hz / (@c CarrierToSymbolRateRatio · @c SymbolRate). Applied values are
 * published via atomics for monitoring.
 *
 * **Sampling-clock offset:** @ref ChannelSamplingClockOffset is **always** started on the output path: each
 * rotated float batch is enqueued there, which resamples according to the current total ppm (channel mismatch
 * + optional Doppler term) then quantizes to shorts. Use @f$\mathrm{ppm}=0@f$ for nominal step @f$1@f$ (no
 * effective rate change, still Lagrange + int16 conversion in that block).
 *
 * **Threads:** noise generation, noisy-sum output, frequency rotation (with an internal queue between
 * stages); exposes thread-safe readout of frequency / ppm / segment for GUIs.
 */
class AWGNChannel
{
private:
    /* data */
    
    #ifdef DEBUG_AWGN
    short *OutAllI = 0;
    unsigned int PtrOutAll = 0;
    #endif
    double EsN0db;
    Params *pParams = 0;
    const double Backoff = 12;
    const int NumBits = 16;
    int NoiseBatchSize;
    random_generator_new oNoiseGen;
    FrequencyOffset objFreqOffset;
    double Stdn; 
    Transmitter *pTx;
    static constexpr int LengthQueue = 4;
    float *Noise[LengthQueue];
    thread NoiseThread, OutputThread, FreqOffsetThread;
    ChannelSamplingClockOffset samplingClockOffset_;
    void GenerateNoise(void);
    void GenerateOutput(void);
    void ApplyFrequencyOffset(void);
    SimpleQueue NoiseQ;
    bool StopAll = false;
    condition_variable CvOutNoise, CvNoiseOut, CvOutUser, CvUserOut;
    // Queue between OutputThread (noise+gain) and FreqOffsetThread (carrier rotation)
    std::mutex mtxFreqQ_;
    std::condition_variable cvFreqQData_;
    std::condition_variable cvFreqQSpace_;
    struct FreqChunk
    {
        std::vector<float> interleaved;
        double freqHz = 0.0;
        double totalPpm = 0.0;
        double gainDb = 0.0;
        int segment = 0; // 0=stable1,1=acc1,2=stable2,3=acc2
        double t_rate = 0.0;
    };
    std::deque<FreqChunk> freqQ_;
    static constexpr int kFreqQDepth = 4;
    __m256i floats_to_shorts_sat_perm_avx2(__m256 x, __m256 y);
    __m256 mkSig;
    double AccelerationPeriod;
    double TotalPeriod, StablePeriod;
    unsigned int TransitionCounter[5];
    double displayPeriodSec_ = 1.0;
    std::mutex mtxOutputBuffer_;
    std::mutex mtxNoiseQ_;

    // Current values actually applied by the channel (updated in ApplyFrequencyOffset thread).
    std::atomic<double> currFreqHz_{0.0};
    std::atomic<double> currTotalPpm_{0.0};
    std::atomic<double> currGainDb_{0.0};
    std::atomic<int> currSegment_{0};
    std::atomic<double> currTRate_{0.0};
public:
    AWGNChannel(unsigned int Seed);
    ~AWGNChannel();
    void StartThreads();
    /** @brief Stop worker threads and join. */
    void StopThreads(void);
    void SetTransmitter(Transmitter *p)
    {
        pTx = p;
    }
    void SetParameters(Params *p)
    {
        pParams = p;
    }
    void GetCvs2User(condition_variable * &CvIn, condition_variable * &CvOut )
    {
        CvOut = &CvOutUser;
        CvIn = &CvUserOut;
    }
    BufferShort OutputBuffer;

    /**
     * @brief Copy samples from @c OutputBuffer (mutex with producer).
     * @param dst Interleaved I,Q shorts.
     * @param nShorts Number of shorts to copy.
     * @param stopAll Stop flag observed while waiting.
     */
    bool CopyOutputSamples(short* dst, int nShorts, bool& stopAll);

    /** @brief Snapshot of channel state for overlays (atomics under internal lock in getter). */
    struct CurrentApplied
    {
        double FreqHz = 0.0;
        double TotalPpm = 0.0;
        double GainDb = 0.0;
        int Segment = 0; // 0=stable1,1=acc1,2=stable2,3=acc2
        double TRate = 0.0;
    };

    CurrentApplied GetCurrentApplied() const
    {
        CurrentApplied s;
        s.FreqHz = currFreqHz_.load(std::memory_order_relaxed);
        s.TotalPpm = currTotalPpm_.load(std::memory_order_relaxed);
        s.GainDb = currGainDb_.load(std::memory_order_relaxed);
        s.Segment = currSegment_.load(std::memory_order_relaxed);
        s.TRate = currTRate_.load(std::memory_order_relaxed);
        return s;
    }
};

