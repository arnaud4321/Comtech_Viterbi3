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
#include "Transmitter.h"
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
    double EmaMetricsGrowth;
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
 * @brief Public structure regrouping all critical statistics and state flags for the Receiver.
 * 
 * @details This structure is intended to be exported to external users or wrappers (e.g. Web GUI, Python bindings)
 * via @ref Receiver::GetRxStatistics(), providing a clean snapshot of the system's performance at a given time.
 */
struct RxStatistics {
    double t_sim;               ///< Wall-clock simulation time elapsed (seconds).
    double t_rate;              ///< Signal time processed, based on estimated symbol rate (seconds).
    
    double samplerMsps;         ///< Sample rate delivered by the USRP sampler (Msps).
    double rxFilterMsps;        ///< Processing throughput of the matched filter (Msps).
    uint64_t uhdOverflows;      ///< Number of UHD hardware overflows encountered.

    bool isSymRateLocked;       ///< True if the Symbol Rate Estimator has detected a stable rate.
    double symRateMsps;         ///< Current estimated symbol rate (Msps).
    
    bool isCentralFreqReady;    ///< True if VIVA (coarse frequency) is producing valid offsets.
    double centralNcoHz;        ///< Current NCO offset applied (Hz).
    double gainDb;              ///< Current AGC gain applied (dB).

    bool isGardnerLocked;       ///< True if the Gardner timing loop has locked its phase.
    double gardnerPpm;          ///< Fine residual sampling clock offset estimated by Gardner (ppm).

    bool isPhaseLocked;         ///< True if the Decision-Directed Phase PLL is locked.
    double phaseFreqHz;         ///< Fine residual frequency offset estimated by the DD-PLL (Hz).
    
    bool isViterbiLocked;       ///< True if the Viterbi decoder has synchronized to a hypothesis.
    /// Number of times Viterbi went from locked to unlocked (e.g. metric-growth desync).
    uint64_t viterbiUnlockEvents = 0;
    /// Number of times Viterbi went from unlocked to locked (includes the initial acquisition).
    uint64_t viterbiRelockEvents = 0;
    bool isPrbsLocked;          ///< True if the descrambled PRBS sequence is locked.

    double evmRms;              ///< Root Mean Square of the Error Vector Magnitude (EVM).
    double snrFromEvmDb;        ///< Signal-to-Noise Ratio (SNR) in dB estimated directly from the EVM.
    
    uint64_t numBits;           ///< Total number of PRBS bits compared.
    uint64_t numErrors;         ///< Total number of PRBS bit errors.
    double ber;                 ///< Current Bit Error Rate (BER) over the whole history since lock.
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
 * @f$F_s/2@f$); accepted estimates update @ref ReceiverResampler::UpdateFromSymbolRateHz.
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
    /// Symboles IQ après phase DD (entrée Viterbi), pour @ref GetTRateRx — jamais remis à zéro par PRBS/Viterbi.
    std::atomic<uint64_t> rxSymStreamTotal_{0};
    int SymRateWindowBatches = 0;
    int SymRateBatchCounter = 0;
    int SymRateAcceptedCount = 0;
    double SymRateEstimatePeriodSec = 1.0;
    double SymRatePeakToMedianThreshold = 8.0;
    double SymRateMaxRelativeJump = 0.02;
    bool SymRateAllowRefinement = true;
    int prbsHighBerStreak_ = 0;
    ReceiverResampler resampler_;
    ReceiverFreqCorrector freqCorrector_;
    ReceiverTimingTracking timingTracking_;
    ReceiverPhaseTrackingDD phaseTrackingDD_;
    Sampler *pSampler;
    Transmitter* pTx_ = nullptr;
    BufferFloat oBufferFilter;
    BufferFloat oBufferResampled{kRxRingFloatLen, kRxRingFloatExtra};
    BufferFloat oBufferFreqCorrected{kRxRingFloatLen, kRxRingFloatExtra};
    Viterbi oViterbi[3] = {Viterbi(0),Viterbi(1),Viterbi(2)};
    std::atomic<bool> StopAll{false};
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
    std::atomic<uint64_t> viterbiUnlockEvents_{0};
    std::atomic<uint64_t> viterbiRelockEvents_{0};
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
    // Raw (pre-Viterbi) SER/BER from hard decisions vs ideal TX symbols (mapped +/-1).
    uint64_t RawNumSymErrorsAll = 0;
    uint64_t RawNumSymsAll = 0;
    uint64_t RawNumErrorsAll = 0;
    uint64_t RawNumBitsAll = 0;
    // Latest SER sync diagnostic (inter-correlation peak).
    std::atomic<double> RawSyncPeakAbs{0.0};
    std::atomic<double> RawSyncPeakPhaseRad{0.0};
    std::atomic<int> RawSyncLagSym{0};
    std::atomic<double> RawSyncBestAbs{0.0};
    std::atomic<double> RawSyncThrAbs{0.0};
    std::atomic<bool> RawSyncLocked{false};
    // Applied alignment used for SER/rawBER (resolves QPSK quadrant ambiguity).
    std::atomic<double> RawSyncAppliedPhaseRad{0.0};
    std::atomic<int> RawSyncAppliedQuad{0}; // 0..3 for +k*(pi/2) added to raw phase
    Receiver(/* args */);
    ~Receiver();

    /** @brief Attach upstream IQ producer (short samples). */
    void SetSampler(Sampler *p)
    {
        pSampler = p;
    }

    /** @brief Attach channel for constellation overlay / applied-parameter readout. */
    void SetChannel(class AWGNChannel* ch) { pChannel_ = ch; }

    /** @brief Attach transmitter for raw pre-Viterbi comparisons (symbol reference). */
    void SetTransmitter(Transmitter* tx) { pTx_ = tx; }

    /**
     * @brief Start all RX threads (filter, resampler, NCO, timing, phase, Viterbi workers).
     * @param RollOff RRC rolloff for receive filter.
     * @param RxModeIn @c PRBS_TX or @c FILE_TX.
     * @param symRateCfg Symbol-rate FFT estimator settings.
     * @param timingCfg Gardner timing PI (@c NominalOmega, @c Kp, @c Ki, @c UpdatePeriodSymbols).
     * @param displayPeriodSec Console cadence (\>0) or event-only logging (\<=0).
     * @param centralFreqCfg Central-frequency (VIVA + NCO) block; may be default-disabled.
     * @param constellationCfg Optional matplotlib/ascii constellation bridge.
     */
    void StartThreads(double RollOff, TxModes RxModeIn,
                      const SymbolRateEstimatorConfig& symRateCfg = SymbolRateEstimatorConfig{},
                      const TimingTrackingConfig& timingCfg = TimingTrackingConfig{},
                      double displayPeriodSec = 1.0,
                      const ReceiverFreqCorrectorConfig& centralFreqCfg = ReceiverFreqCorrectorConfig{},
                      const PhaseTrackingConfig& phaseCfg = PhaseTrackingConfig{},
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
    double GetPhaseTrackingEvmRms() const { return phaseTrackingDD_.GetLastEvmRms(); }

    /**
     * @brief Collects and returns a snapshot of all system-level performance metrics.
     * @param t_sim_optional Externally measured elapsed wall-clock simulation time to inject (if known).
     * @return Fully populated @ref RxStatistics struct.
     */
    RxStatistics GetRxStatistics(double t_sim_optional = 0.0) const;

    /**
     * @brief Temps « signal » écoulé : symboles reçus (après phase DD) / débit symbole estimé.
     * @details Utilise @ref SymbolRateEstimateHz dès qu’il est \> 1 Hz, sinon @c SymbolRate nominal.
     *          Indépendant des lock/désync Viterbi et PRBS (contrairement à @ref NumBitsAll).
     */
    double GetTRateRx() const
    {
        const uint64_t n = rxSymStreamTotal_.load(std::memory_order_relaxed);
        const double est = SymbolRateEstimateHz.load(std::memory_order_relaxed);
        const double denom = (est > 1.0) ? est : SymbolRate;
        return static_cast<double>(n) / denom;
    }
};
