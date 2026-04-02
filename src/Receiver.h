/**
 * @file Receiver.h
 * @brief QPSK receiver orchestration: filtering, resampling, frequency/timing/phase, Viterbi, output.
 */
#pragma once
#include "SimpleQueue.h"
#include "definitions.h"
#include <cstdint>
#include <thread>
#include "Sampler.h"
#include "definitions.h"
#include <cstdio>
#include <immintrin.h>
#include <condition_variable>
#include <mutex>
#include <immintrin.h>
#include "RxFilter.h"
#include "BufferFloat.h"
#include "Viterbi.h"
#include "SelfSyncScrambler_V35.h"
#include "Prbs23.h"
#include "SymbolRateEstimator.h"
#include "ReceiverResampler.h"
#include "ReceiverTimingTracking.h"
#include "ReceiverPhaseTrackingDD.h"
#include "ReceiverFreqCorrector.h"
#include "ConstellationDisplay.h"
#include <memory>
#include <atomic>
using namespace std;
//#define DEBUG1
//#define DEBUG2
#define DEBUG_STATISTICS

// Uncomment to bypass Gardner (TakeEven 2x->1 sps), dump Viterbi input symbols to ../data/viterbi_input_bypass.bin
// (float32 interleaved I,Q per symbol), then compare offline with Gardner on-time dump via compare_viterbi_gardner_dump.m
//#define BYPASS_GARDNER_DUMP_VITERBI_INPUT

class Sampler;

/** @brief Optional debug aggregates for Viterbi metric growth (when @c DEBUG_STATISTICS is on). */
struct DebugStatistics
{
    double FirstDelta;
    double MaxMetricsGrowth;
    double SumMetricsGrowth;
    double MeanMetricsGrowth;
    uint64_t NumBatches;
};


/**
 * @brief Per-decoder-branch metrics from the pre-lock Viterbi alignment search.
 *
 * @details @ref Receiver keeps @c VitSyncResults[3], one entry per @ref Viterbi worker index. Until
 * @c ViterbiSynchronized is true, @c OperateViterbi runs @ref Viterbi::Decode four times on the same symbol
 * block with different IQ permutations (see @c Metrics indices below), resets the trellis between probes, and
 * stores the **final path metric** (or growth proxy) from each decode in @c Metrics. The manager thread waits
 * until all three branches set @c Available, then compares metrics across branches and hypotheses to choose
 * @ref ViterbiParameters and assert lock. After lock, @c Metrics[0] is still updated on the locked decode path
 * for optional statistics.
 *
 * Probe mapping (same order as @c OperateViterbi): @n
 * @c Metrics[0]: @c ExchangeIQ=false,  @c SignI=+1, @c SignQ=+1 @n
 * @c Metrics[1]: @c ExchangeIQ=false,  @c SignI=+1, @c SignQ=-1 @n
 * @c Metrics[2]: @c ExchangeIQ=true,   @c SignI=-1, @c SignQ=+1 @n
 * @c Metrics[3]: @c ExchangeIQ=true,   @c SignI=+1, @c SignQ=+1 @n
 *
 * @see @ref ViterbiParameters, @ref Receiver, @ref Receiver.cpp @c OperateViterbi, @c OperateViterbiManager.
 */
struct ViterbiSyncResults
{
    bool Available = false; ///< Set true after this branch’s four probe decodes complete for the current slot.
    float Metrics[4]{};     ///< One metric per IQ/sign hypothesis; lower is better for best-path selection.
};

/**
 * @brief IQ swap and I/Q sign flips applied by @ref Viterbi::Decode after alignment lock.
 *
 * @details Chosen by @c OperateViterbiManager from @ref ViterbiSyncResults across three decoders. Matches one
 * of the four probe configurations (see @ref ViterbiSyncResults).
 *
 * @see @ref ViterbiSyncResults
 */
struct ViterbiParameters
{
    bool ExchangeIQ; ///< If true, I/Q are exchanged before decoding (spectrum / rotation family).
    float SignI;     ///< Multiply I by +1 or -1 before decode.
    float SignQ;     ///< Multiply Q by +1 or -1 before decode.
};

/**
 * @brief End-to-end receive chain and worker threads (includes **Viterbi alignment**: three decoders, four IQ/sign
 * metric probes until lock, then merged output — see **Viterbi manager thread** below).
 *
 * **Steady-state data path:** matched filter ring → @ref ReceiverResampler → @ref ReceiverFreqCorrector →
 * @ref ReceiverTimingTracking → @ref ReceiverPhaseTrackingDD → soft symbols → @ref Viterbi → descramble/PRBS
 * check → @c OutputQ.
 *
 * @details **Filter thread:** pulls @c short IQ batches from @ref Sampler, converts to float (AVX), runs
 * @ref RxFilter into the filter ring write window (@c SPB complex samples per batch). Until
 * @c SymbolRateDetected, @c AdvancePtrWr is **not** called: batches still drive @ref SymbolRateEstimator
 * during collection windows, but the ring does not advance toward @ref ReceiverResampler (the same slot may be
 * overwritten each iteration). After the first accepted estimate, each batch commits @c SPB samples (plus a
 * one-shot advance if detection occurred without a commit in the same pass). Symbol-rate windows are **not**
 * started while Gardner (@ref ReceiverTimingTracking) is locked; re-estimation is similarly gated.
 *
 * A @ref SymbolRateEstimator periodically estimates symbol rate (@f$|r[n]|@f$ → real FFT, peak search near
 * @f$F_s/2@f$); accepted estimates update @ref ReceiverResampler::UpdateFromSymbolRateHz so downstream
 * stages see ~2 samples per symbol at the tracked rate.
 *
 * **Frequency / timing / phase:** @ref ReceiverFreqCorrector optionally runs @ref CentralFreqEstimatorViva
 * on blocks of 2-sps data, drives an NCO (@ref FrequencyOffset) and AGC before @ref ReceiverTimingTracking
 * (see @ref GardnerTiming). @ref ReceiverPhaseTrackingDD removes residual carrier phase with a
 * decision-directed loop, then symbols are fanned out to three @ref Viterbi decoders.
 *
 * **Viterbi manager thread (alignment / lock):** before lock, each @c OperateViterbi worker computes four metrics
 * (@c ExchangeIQ / sign combinations) into @ref ViterbiSyncResults; the manager picks the lowest-metric hypothesis,
 * sets @ref ViterbiParameters, and sets @c ViterbiSynchronized true. Then: merge the three decoded streams,
 * differential decode, PRBS lock, @ref SelfSyncScrambler_V35 into @c OutputQ.
 * (See @ref Receiver.cpp @c OperateViterbi / @c OperateViterbiManager.)
 */
class Receiver
{
private:
#ifdef DEBUG1
    float *OutAllI = 0;
    float *OutAllQ = 0;
    unsigned int PtrOutAll = 0;
    #define NUM_SAMPLES_R1 19000000
#endif

    Prbs23 oPrbs;
    SelfSyncScrambler_V35 objDescrambler;
    RxFilter objRxFilter;
    SymbolRateEstimator objSymRateEstimator;
    std::atomic<bool> SymbolRateDetected{false};
    std::atomic<double> SymbolRateEstimateHz{0.0};
    int SymRateWindowBatches = 0;
    int SymRateBatchCounter = 0;
    int SymRateAcceptedCount = 0;
    double SymRateEstimatePeriodSec = 1.0;
    double SymRatePeakToMedianThreshold = 8.0;
    double SymRateMaxRelativeJump = 0.02;
    ReceiverResampler resampler_;
    ReceiverFreqCorrector freqCorrector_;
    ReceiverTimingTracking timingTracking_;
    ReceiverPhaseTrackingDD phaseTrackingDD_;
    Sampler *pSampler;
    BufferFloat oBufferFilter;
    BufferFloat oBufferResampled{SPB*256,32*SPB};
    BufferFloat oBufferFreqCorrected{SPB*256,32*SPB};
    Viterbi oViterbi[3] = {Viterbi(0),Viterbi(1),Viterbi(2)};
    bool StopAll;
    float *SplitI[3], *SplitQ[3];
    unsigned char *Outputs[3];
    void OperateFilter(void);
    int IndexViterbi[3] = {0,1,2};
    void OperateViterbiManager(void);
    void OperateViterbi(void *p);
    thread FilterThread, ViterbiManagerThread;
    thread ViterbiThreads[3];
    ViterbiSyncResults VitSyncResults[3];
    float *FilterInI,  *FilterInQ;
    void shorts_to_floats_avx2(__m256i v16, __m256* out0, __m256* out1);
    condition_variable CvFilterUser,CvRx2Out,CvOut2Rx;
    condition_variable CvResampData;
    condition_variable CvFreqCorrData;
    condition_variable CvVitManager2Vit;
    condition_variable CvViterbis2VitManager[3];
    std::mutex mtxFilterRing;
    std::mutex mtxResampRing;
    std::mutex mtxFreqCorrRing;
    // Protects all SimpleQueue (DemodulatorQ, DecodedQ, OutputQ) + associated flags.
    std::mutex mtxQueues;

    void Split3(float *InputI, float *InputQ, int Length, int PtrWr);
    void TakeEvenDebug(float *Input, float *Output, int Length);
    static constexpr int LengthQueue = 32;
    SimpleQueue DecodedQ[3] = {SimpleQueue(LengthQueue),SimpleQueue(LengthQueue),SimpleQueue(LengthQueue)};
    SimpleQueue DemodulatorQ[3] = {SimpleQueue(LengthQueue),SimpleQueue(LengthQueue),SimpleQueue(LengthQueue)};
    SimpleQueue OutputQ;
    // Protected by mtxQueues (all read/write must be done under lock).
    bool ViterbiSynchronized = false, PRBSSynchronized = false;
    ViterbiParameters ViterbiParams;
    unsigned char *ViterbiOutputs[3];
    unsigned char *DiffDec[3];
    unsigned char *OutputAll;
    unsigned char *Merged;
    unsigned char *PrbsOut;
    void DifferentialDecode(unsigned char *In, unsigned char *Out, unsigned int Length, unsigned char &Last);
    bool SyncPRBS(unsigned char *Input, int &PtrStart, unsigned int& Seed, int& NumErrorsAtLock);
    void CountErrors(unsigned char *Input, unsigned char *Template, int Length);
    TxModes RxMode;
    double displayPeriodSec_ = 1.0;
    double rxThroughputMeasurePeriodSec_ = 1.0;
    ConstellationDisplayConfig constellationCfg_{};
    std::unique_ptr<ConstellationDisplay> constellationDisplay_;
    class AWGNChannel* pChannel_ = nullptr;
    std::atomic<double> lastRxFilterThroughputMsps_{0.0};
public:
#ifdef DEBUG_STATISTICS

    DebugStatistics CurrDebugStatistics;
#endif
    Receiver(/* args */);
    ~Receiver();

    /** @brief Attach upstream IQ producer (short samples). */
    void SetSampler(Sampler *p)
    {
        pSampler = p;
    }

    /** @brief Attach channel for constellation overlay / applied-parameter readout. */
    void SetChannel(class AWGNChannel* ch) { pChannel_ = ch; }

    /**
     * @brief Start all RX threads (filter, resampler, NCO, timing, phase, Viterbi workers).
     * @param RollOff RRC rolloff for receive filter.
     * @param RxModeIn @c PRBS_TX or @c FILE_TX.
     * @param symRateCfg Symbol-rate FFT estimator settings.
     * @param displayPeriodSec Console cadence (\>0) or event-only logging (\<=0).
     * @param centralFreqCfg Central-frequency (VIVA + NCO) block; may be default-disabled.
     * @param constellationCfg Optional matplotlib/ascii constellation bridge.
     */
    void StartThreads(double RollOff, TxModes RxModeIn,
                      const SymbolRateEstimatorConfig& symRateCfg = SymbolRateEstimatorConfig{},
                      double displayPeriodSec = 1.0,
                      const ReceiverFreqCorrectorConfig& centralFreqCfg = ReceiverFreqCorrectorConfig{},
                      const ConstellationDisplayConfig& constellationCfg = ConstellationDisplayConfig{});

    /** @brief Signal shutdown and join worker threads. */
    void StopThreads(void);

    /** @brief Non-blocking peek at next descrambled byte triple (under queue mutex). */
    unsigned char *GetOutput()
    {
        unsigned char *RetVal = 0;
        std::lock_guard<std::mutex> lk(mtxQueues);
        if(OutputQ.AvailableRead())
        {
            int PtrRd = static_cast<int>(OutputQ.GetPtrRd());
            RetVal = OutputAll + PtrRd*BatchSize3;
        }
        return RetVal;
    }

    /** @brief Release current output slot after @ref GetOutput consumption. */
    void AdvanceOutput(void)
    {
        {
            std::lock_guard<std::mutex> lk(mtxQueues);
            OutputQ.AdvanceRead();
        }
        CvOut2Rx.notify_one();
    }

    /** @brief Legacy hook: output-side condition variables (FILE_TX path). */
    void Get_Condition_Variables(condition_variable * pcv_rx_out, condition_variable * pcv_out_rx)
    {
        pcv_out_rx = &CvOut2Rx;
        pcv_rx_out = &CvRx2Out;
    }
    uint64_t NumErrorsAll, NumBitsAll;
};
