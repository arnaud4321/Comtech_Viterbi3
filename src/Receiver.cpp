#include <deque>
/**
 * @file Receiver.cpp
 * @brief Receive chain: filter thread, symbol-rate estimate, resampler, NCO/AGC, Gardner, phase DD, three Viterbi, output.
 *
 * @details **Startup (@c StartThreads)** — Configures @ref RxFilter, @ref SymbolRateEstimator, optional
 * @ref ConstellationDisplay; starts @ref ReceiverResampler on @c oBufferFilter → @c oBufferResampled;
 * @ref ReceiverFreqCorrector on resampled → @c oBufferFreqCorrected; @ref ReceiverTimingTracking on that ring;
 * @ref ReceiverPhaseTrackingDD fed by timing; then **FilterThread** (@c OperateFilter), **ViterbiManagerThread**
 * (@c OperateViterbiManager), three **OperateViterbi** workers.
 *
 * **OperateFilter** — @ref Sampler::ReadFilterBatch → AVX shorts-to-float → @ref RxFilter → @c oBufferFilter
 * write window (@c SPB pairs). Until @c SymbolRateDetected, @c AdvancePtrWr is skipped: @ref SymbolRateEstimator
 * still receives filter output during timed collection windows when Gardner is **not** locked; the resampler
 * does not consume. On first acceptance, @c SymbolRateEstimateHz and @c ReceiverResampler::UpdateFromSymbolRateHz;
 * then each batch advances the ring by @c SPB (with catch-up @c AdvancePtrWr if detection flipped mid-pass).
 * Optional smoothed re-estimates when not locked. Notifies @c CvFilterUser / @c CvResampData; throughput and
 * constellation hooks as implemented.
 *
 * **OperateViterbi** (×3) — Waits on @c DemodulatorQ[i]. If not yet synchronized, runs four metric probes
 * (IQ swap/sign combinations) into @c VitSyncResults; else decodes with locked @c ViterbiParams, differential
 * decode into @c DecodedQ[i].
 *
 * **OperateViterbiManager** — @c WaitPopSymbolFrame from phase DD; @c Split3 to three alignments; enqueue
 * @c DemodulatorQ; collects sync metrics, picks best hypothesis, merges bytes, PRBS sync, descramble,
 * @c OutputQ. Optional constellation @c UpdateEx and periodic status.
 *
 * **StopThreads** — Join order: filter → @c resampler_ → @c freqCorrector_ → @c timingTracking_ →
 * @c phaseTrackingDD_ → three Viterbi workers → manager (releases constellation child).
 */

#include "Receiver.h"
#include "AWGNChannel.h"
#include "ConsoleAlert.h"
#include <chrono>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <limits>
#include <sys/stat.h>

extern std::mutex mtxfilethr;

// Forward declaration: used in OperateViterbiManager before the definition further below.
static void InjectPRBSErrorDeterministic(unsigned char* buf, int length);

namespace
{
} // namespace

Receiver::Receiver(/* args */):oBufferFilter(kRxRingFloatLen, kRxRingFloatExtra),OutputQ(LengthQueue)
{
    FilterInI = (float*) _mm_malloc(SPB*2*sizeof(float),32);
    FilterInQ = FilterInI + SPB;
    for(int i = 0; i < 3; i++)
    {
        SplitI[i] = (float*) _mm_malloc(BatchSize1*LengthQueue*sizeof(float),32);
        SplitQ[i] = (float*) _mm_malloc(BatchSize1*LengthQueue*sizeof(float),32);
        ViterbiOutputs[i] = (unsigned char *) _mm_malloc((BatchSize1+32),32);
        DiffDec[i] = (unsigned char *) _mm_malloc((BatchSize1*LengthQueue),32);

    }
    
    
    Outputs[0] = (unsigned char*) _mm_malloc(BatchSize3,32);




    Merged = (unsigned char *) _mm_malloc((BatchSize3),32);
    OutputAll = (unsigned char *) _mm_malloc((BatchSize3)*LengthQueue,32);
    PrbsOut = (unsigned char *) _mm_malloc((BatchSize3),32);

    #ifdef DEBUG1
    OutAllI = (float*) _mm_malloc( (NUM_SAMPLES_R1+100000)*sizeof(float),32);
    OutAllQ = (float*) _mm_malloc( (NUM_SAMPLES_R1+100000)*sizeof(float),32);
#endif



}

Receiver::~Receiver()
{
    _mm_free(FilterInI);
    _mm_free(Merged);
    _mm_free(OutputAll);
    _mm_free(PrbsOut);

    for(int i = 0; i < 3; i++)
    {
        _mm_free(SplitI[i]);
        _mm_free(SplitQ[i]);
        _mm_free(ViterbiOutputs[i]);
        _mm_free(DiffDec[i]);

    }

    #ifdef DEBUG1
        _mm_free(OutAllI);
        _mm_free(OutAllQ);
    #endif

}
void Receiver::StartThreads(double RollOff, TxModes RxModeIn,
                            const SymbolRateEstimatorConfig& symRateCfg,
                            double displayPeriodSec,
                            const ReceiverFreqCorrectorConfig& centralFreqCfg,
                            const ConstellationDisplayConfig& constellationCfg)
{
    StopAll = false;
    RxMode = RxModeIn;
    displayPeriodSec_ = displayPeriodSec;
    constellationCfg_ = constellationCfg;
    // Keep measure cadence at >= 1 s when console is off: short PeriodSec (e.g. 0.5) was tied here and
    // differed from DisplayPeriodSec==1 runs (see sampler tpm in main.cpp too).
    rxThroughputMeasurePeriodSec_ =
        (displayPeriodSec > 0.0) ? displayPeriodSec : 1.0;
    constellationDisplay_.reset();
    if (constellationCfg_.PeriodSec > 0.0 && constellationCfg_.Backend == "matplotlib")
    {
        constellationDisplay_ = std::make_unique<ConstellationDisplay>();
        if (!constellationDisplay_->Start(constellationCfg_))
        {
            std::cout << "[Constellation] \033[31mDisabled\033[0m (failed to start backend)"
                      << " (Backend=" << constellationCfg_.Backend << ")"
                      << std::endl;
            constellationDisplay_.reset();
        }
        else
        {
            std::cout << "[Constellation] \033[32mEnabled\033[0m"
                      << " backend=" << constellationCfg_.Backend
                      << " period=" << constellationCfg_.PeriodSec << " s"
                      << " N=" << constellationCfg_.NumSymbols
                      << " XDisplay=" << (constellationCfg_.XDisplay.empty() ? "(inherit)" : constellationCfg_.XDisplay)
                      << " pid=" << (constellationDisplay_ ? constellationDisplay_->GetChildPid() : -1)
                      << std::endl;
        }
    }
    oBufferFilter.Reset();
    objRxFilter.CreateObjects(RollOff);
    ViterbiSynchronized = false;
    PRBSSynchronized = false;
    RawNumBitsAll = 0;
    RawNumErrorsAll = 0;
    for(int i = 0; i < 3; i++)
    {
        DemodulatorQ[i].Reset();
        DecodedQ[i].Reset();
    }
    SymRatePeakToMedianThreshold = symRateCfg.PeakToMedianThreshold;
    SymRateMaxRelativeJump = symRateCfg.MaxRelativeJump;
    SymRateEstimatePeriodSec = std::max(0.01, symRateCfg.EstimatePeriodSec);
    objSymRateEstimator.Reset(SamplingFrequency,
                              symRateCfg.FftSize,
                              symRateCfg.MaxOffsetHz,
                              symRateCfg.PeakToMedianThreshold);
    SymRateWindowBatches = objSymRateEstimator.GetRequiredBatches(SPB);
    SymRateBatchCounter = 0;
    SymRateAcceptedCount = 0;
    SymbolRateDetected.store(false, std::memory_order_relaxed);
    SymbolRateEstimateHz.store(0.0, std::memory_order_relaxed);
    oBufferResampled.Reset();
    oBufferFreqCorrected.Reset();
    resampler_.Start(&oBufferFilter, &mtxFilterRing, &CvFilterUser,
                     &oBufferResampled, &mtxResampRing, &CvResampData,
                     &CvResampData, &StopAll);
    freqCorrector_.Configure(centralFreqCfg, displayPeriodSec);
    timingTracking_.ResetGardner(2.0, 1.0e-4, 1.0e-6, 64);
    // Central frequency correction (NCO) between resampler (2 sps) and Gardner.
    // Config is provided via Params in main (see StartThreads call).
    freqCorrector_.Start(&oBufferResampled, &mtxResampRing, &CvResampData,
                         &oBufferFreqCorrected, &mtxFreqCorrRing, &CvFreqCorrData,
                         &phaseTrackingDD_, &StopAll);
    timingTracking_.Start(&oBufferFreqCorrected, &mtxFreqCorrRing, &CvFreqCorrData, &StopAll, displayPeriodSec);
    phaseTrackingDD_.Start(&timingTracking_, &StopAll, displayPeriodSec);

    FilterThread = std::thread(&Receiver::OperateFilter, this);
    ViterbiManagerThread = std::thread(&Receiver::OperateViterbiManager, this);

    for(int i = 0; i < 3; i++)
    {
        void *p = (void*) (IndexViterbi+i);
        ViterbiThreads[i] = std::thread(&Receiver::OperateViterbi,this,p );
    }
}
void Receiver::StopThreads(void)
{
    StopAll = true;
    if (pSampler)
        pSampler->NotifyFilterWaiters();
    CvRx2Out.notify_all();
    CvOut2Rx.notify_all();
    CvVitManager2Vit.notify_all();
    CvFilterUser.notify_all();
    CvResampData.notify_all();
    CvFreqCorrData.notify_all();
    if(FilterThread.joinable())
        FilterThread.join();
    resampler_.StopJoin();
    freqCorrector_.StopJoin();
    timingTracking_.StopJoin();
    phaseTrackingDD_.StopJoin();
    if (constellationDisplay_)
        constellationDisplay_->ClosePipeKeepAlive();
    constellationDisplay_.reset();
    for(int i = 0; i < 3;i++)
        CvViterbis2VitManager[i].notify_one();
    for(int i = 0; i < 3;i++)
    {
        if(ViterbiThreads[i].joinable())
            ViterbiThreads[i].join();
    }
    if(ViterbiManagerThread.joinable())
        ViterbiManagerThread.join();
}

/**
 * @brief RX filter thread: IQ shorts → matched filter ring, symbol-rate acquisition, resampler feed.
 */
void Receiver::OperateFilter(void)
{
    auto throughput_window_start = std::chrono::steady_clock::now();
    auto rx_filter_console_window_start = std::chrono::steady_clock::now();
    uint64_t samples_in_window = 0;
    uint64_t samples_total = 0;
    const auto wall_start = std::chrono::steady_clock::now();
    auto symrate_display_start = std::chrono::steady_clock::now();
    auto next_symrate_start = std::chrono::steady_clock::now();
    bool symrate_collecting = false;
    auto lastRxFilterRingSatLog = std::chrono::steady_clock::now();

    #ifdef WRITE_LOG_THR
	mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"Receiver OperateFilter Thread %d\n", gettid());
    fclose(fidthr);
    mtxfilethr.unlock();
    #endif

    alignas(32) short samplerChunk[2 * SPB];
    while (!StopAll)
    {
        if (!pSampler->ReadFilterBatch(samplerChunk, 2 * SPB, StopAll))
            break;
        __m256 mouti;
        __m256 moutq;
        int PtrOut = 0;
        for (int i = 0; i < (2 * SPB); i += 16)
        {
            __m256i mIn = _mm256_loadu_si256((__m256i*)(samplerChunk + i));
           shorts_to_floats_avx2(mIn,&mouti,&moutq);
           _mm256_store_ps(FilterInI+PtrOut,mouti);
           _mm256_store_ps(FilterInQ+PtrOut,moutq);
           PtrOut += 8;
        }
        auto now_sym = std::chrono::steady_clock::now();

        bool filterAdvancedRing = false;
        bool symRateEstimateReady = false;
        {
            std::unique_lock<std::mutex> lk(mtxFilterRing);
            if (SymbolRateDetected.load(std::memory_order_relaxed))
            {
                if (oBufferFilter.AlmostFull())
                {
                    const auto nowSat = std::chrono::steady_clock::now();
                    if (nowSat - lastRxFilterRingSatLog >= std::chrono::seconds(1))
                    {
                        lastRxFilterRingSatLog = nowSat;
                        CONSOLE_ALERT_STMT(std::cout << ConsoleAlert::kRedOpen
                                                     << "[RXFilter] matched-filter float ring almost full; fill="
                                                     << oBufferFilter.GetSizeInBuffer() << "/"
                                                     << oBufferFilter.GetBufferSize() - 1 << ConsoleAlert::kReset
                                                     << std::endl;);
                    }
                }
                CvFilterUser.wait(lk, [&] {
                    return StopAll || !oBufferFilter.AlmostFull();
                });
            }
            if (StopAll)
                break;

            float *FilterOutI, *FilterOutQ;
            oBufferFilter.GetWriteBuffer(FilterOutI, FilterOutQ, SPB);
            objRxFilter.CreateOutputs(FilterInI, FilterInQ, FilterOutI, FilterOutQ, SPB);
            samples_total += static_cast<uint64_t>(SPB);

            if (!symrate_collecting && now_sym >= next_symrate_start && !timingTracking_.IsLocked())
            {
                symrate_collecting = true;
                SymRateBatchCounter = 0;
                next_symrate_start += std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                    std::chrono::duration<double>(SymRateEstimatePeriodSec));
                while (next_symrate_start <= now_sym)
                {
                    next_symrate_start += std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                        std::chrono::duration<double>(SymRateEstimatePeriodSec));
                }
            }

            if (symrate_collecting)
            {
                if (timingTracking_.IsLocked())
                {
                    // Do not build a new symbol-rate estimate while Gardner is locked.
                    symrate_collecting = false;
                    SymRateBatchCounter = 0;
                    symRateEstimateReady = false;
                }
                else
                {
                    symRateEstimateReady = objSymRateEstimator.PushBatch(FilterOutI, FilterOutQ, SPB);
                    SymRateBatchCounter++;
                }
            }

            if (SymbolRateDetected.load(std::memory_order_relaxed))
            {
                assert(oBufferFilter.CanCommitWrite(SPB));
                oBufferFilter.AdvancePtrWr(SPB);
                filterAdvancedRing = true;
            }
        }

        if (symrate_collecting && symRateEstimateReady)
        {
            SymbolRateEstimateResult est = objSymRateEstimator.RunEstimation();
            bool newSymRateEstimate = false;
            if (!SymbolRateDetected.load(std::memory_order_relaxed))
            {
                if (est.Detected)
                {
                    SymbolRateEstimateHz.store(est.SymbolRateHz, std::memory_order_relaxed);
                    SymbolRateDetected.store(true, std::memory_order_relaxed);
                    SymRateAcceptedCount = 1;
                    resampler_.UpdateFromSymbolRateHz(est.SymbolRateHz);
                    newSymRateEstimate = true;
                    const double t_rate = static_cast<double>(samples_total) / SamplingFrequency;
                    const double t_sim =
                        std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
                    std::cout << "[SymbolRateEstimator] \033[32mDetected\033[0m: " << (est.SymbolRateHz / 1e6)
                              << " Msym/s"
                              << " [window " << SymRateBatchCounter << "/" << SymRateWindowBatches << "]"
                              << " t_sim=" << t_sim << " s"
                              << " t_rate=" << t_rate << " s"
                              << " (search kMax=" << est.SearchKMaxBins
                              << ", fMax=" << est.SearchMaxOffsetHz << " Hz)"
                              << " (FFT res " << est.FftResolutionHz
                              << " Hz, peak/median " << est.PeakToMedian
                              << ", threshold " << SymRatePeakToMedianThreshold << ")" << std::endl;
                }
                else
                {
                    const double t_rate = static_cast<double>(samples_total) / SamplingFrequency;
                    const double t_sim =
                        std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
                    std::cout << "[SymbolRateEstimator] \033[31mNot detected\033[0m"
                              << " [window " << SymRateBatchCounter << "/" << SymRateWindowBatches << "]"
                              << " t_sim=" << t_sim << " s"
                              << " t_rate=" << t_rate << " s"
                              << " (reason: " << (est.FailReason ? est.FailReason : "") << ")"
                              << " (search kMax=" << est.SearchKMaxBins
                              << ", fMax=" << est.SearchMaxOffsetHz << " Hz"
                              << ", dF=" << est.FftResolutionHz << " Hz"
                              << ", peak/median=" << est.PeakToMedian
                              << ", threshold=" << SymRatePeakToMedianThreshold << ")"
                              << std::endl;
                }
            }
            else if (est.Detected)
            {
                if (timingTracking_.IsLocked())
                {
                    // Skip re-estimation while Gardner is locked to avoid delocking transients.
                    symrate_collecting = false;
                    SymRateBatchCounter = 0;
                }
                else
                {
                    const double prev = SymbolRateEstimateHz.load(std::memory_order_relaxed);
                    const double rel = std::abs(est.SymbolRateHz - prev) / std::max(1.0, prev);
                    if (est.PeakToMedian >= SymRatePeakToMedianThreshold && rel <= SymRateMaxRelativeJump)
                    {
                        const double alpha = 1.0 / static_cast<double>(std::max(2, SymRateAcceptedCount + 1));
                        const double refined = (1.0 - alpha) * prev + alpha * est.SymbolRateHz;
                        SymbolRateEstimateHz.store(refined, std::memory_order_relaxed);
                        SymRateAcceptedCount++;
                        resampler_.UpdateFromSymbolRateHz(refined);
                        newSymRateEstimate = true;
                        const double t_rate = static_cast<double>(samples_total) / SamplingFrequency;
                        const double t_sim =
                            std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
                        const double delta_ppm = (refined - prev) * 1.0e6 / std::max(1.0, prev);
                        std::cout << "[SymbolRateEstimator] \033[32mRe-estimated\033[0m: " << (refined / 1e6)
                                  << " Msym/s"
                                  << " [window " << SymRateBatchCounter << "/" << SymRateWindowBatches << "]"
                                  << " t_sim=" << t_sim << " s"
                                  << " t_rate=" << t_rate << " s"
                                  << " (delta " << delta_ppm << " ppm, peak/median " << est.PeakToMedian << ")"
                                  << std::endl;
                    }
                }
            }
            else
            {
                const double t_rate = static_cast<double>(samples_total) / SamplingFrequency;
                const double t_sim =
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
                std::cout << "[SymbolRateEstimator] \033[31mRe-estimation not detected\033[0m"
                          << " [window " << SymRateBatchCounter << "/" << SymRateWindowBatches << "]"
                          << " t_sim=" << t_sim << " s"
                          << " t_rate=" << t_rate << " s"
                          << " (reason: " << (est.FailReason ? est.FailReason : "") << ")"
                          << " (search kMax=" << est.SearchKMaxBins
                          << ", fMax=" << est.SearchMaxOffsetHz << " Hz"
                          << ", dF=" << est.FftResolutionHz << " Hz"
                          << ", peak/median=" << est.PeakToMedian << ")"
                          << std::endl;
            }
            symrate_collecting = false;
            SymRateBatchCounter = 0;

            // Only display the "current estimate" line when a new estimate was accepted.
            if (newSymRateEstimate)
            {
                std::cout << "Current symbol rate estimate: "
                          << (SymbolRateEstimateHz.load(std::memory_order_relaxed) / 1e6)
                          << " Msym/s (window " << SymRateWindowBatches << " batches)" << std::endl;
            }
        }

        if (!filterAdvancedRing && SymbolRateDetected.load(std::memory_order_relaxed))
        {
            std::unique_lock<std::mutex> lk(mtxFilterRing);
            assert(oBufferFilter.CanCommitWrite(SPB));
            oBufferFilter.AdvancePtrWr(SPB);
            filterAdvancedRing = true;
        }

        if (std::chrono::duration<double>(now_sym - symrate_display_start).count() >= 1.0)
        {
            if (!SymbolRateDetected.load(std::memory_order_relaxed))
                std::cout << "Waiting for symbol-rate detection..." << std::endl;
            symrate_display_start = now_sym;
        }

        if (!SymbolRateDetected.load(std::memory_order_relaxed))
            continue;

        if (filterAdvancedRing)
        {
            samples_in_window += SPB;

            auto now_tp = std::chrono::steady_clock::now();
            const double dtMeas =
                std::chrono::duration<double>(now_tp - throughput_window_start).count();
            if (dtMeas >= rxThroughputMeasurePeriodSec_)
            {
                const double msps = static_cast<double>(samples_in_window) / dtMeas / 1e6;
                lastRxFilterThroughputMsps_.store(msps, std::memory_order_relaxed);
                throughput_window_start = std::chrono::steady_clock::now();
                samples_in_window = 0;
            }
            const double dtConsole =
                std::chrono::duration<double>(now_tp - rx_filter_console_window_start).count();
            if (displayPeriodSec_ > 0.0 && dtConsole >= displayPeriodSec_)
            {
                std::cout << "[RXFilter] throughput=" << lastRxFilterThroughputMsps_.load(std::memory_order_relaxed)
                          << " Msps" << std::endl;
                rx_filter_console_window_start = std::chrono::steady_clock::now();
            }

            CvFilterUser.notify_one();
        }
    }
}

void Receiver::DifferentialDecode(unsigned char *In, unsigned char *Out, unsigned int Length, unsigned char &Last)
{
    In[0] = Last;
    for(int i = 0; i < Length; i+=32)
    {
        __m256i mIn = _mm256_loadu_si256((__m256i*) (In+i));
        __m256i mIn1 = _mm256_loadu_si256((__m256i*) (In+i+1));
        __m256i mxIn = _mm256_xor_si256(mIn , mIn1);
        _mm256_storeu_si256((__m256i*) (Out+i), mxIn);
    }

    Last = In[Length];
}
/**
 * @brief One Viterbi worker: hypothesis @a p indexes @c DemodulatorQ / @c oViterbi slot.
 */
void Receiver::OperateViterbi(void *p)
{
    int *ip = (int*) p;
    int IndexViterbi = *ip;
    unsigned char Last = 0;


    #ifdef WRITE_LOG_THR
	mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"Receiver OperateViterbi Thread %d %d\n",IndexViterbi, gettid());
    fclose(fidthr);
    mtxfilethr.unlock();
    #endif


    while(!StopAll)
    {
        std::unique_lock<std::mutex> lk(mtxQueues);
        CvVitManager2Vit.wait(lk, [&] {
            return StopAll || DemodulatorQ[IndexViterbi].AvailableRead();
        });
        if(StopAll)
            break;
        int PtrRd = DemodulatorQ[IndexViterbi].GetPtrRd()*BatchSize1;
        if(!ViterbiSynchronized)
        {
            lk.unlock();
            oViterbi[IndexViterbi].reset_decoder();
            float m[4] = {};
            oViterbi[IndexViterbi].Decode(SplitI[IndexViterbi]+PtrRd, SplitQ[IndexViterbi] + PtrRd, BatchSize1, ViterbiOutputs[IndexViterbi]+1, false, 1, 1, m[0]);
            oViterbi[IndexViterbi].reset_decoder();
            oViterbi[IndexViterbi].Decode(SplitI[IndexViterbi]+PtrRd, SplitQ[IndexViterbi] + PtrRd, BatchSize1, ViterbiOutputs[IndexViterbi]+1, false, 1, -1, m[1]);
            oViterbi[IndexViterbi].reset_decoder();
            oViterbi[IndexViterbi].Decode(SplitI[IndexViterbi]+PtrRd, SplitQ[IndexViterbi] + PtrRd, BatchSize1, ViterbiOutputs[IndexViterbi]+1, true, -1, 1, m[2]);
            oViterbi[IndexViterbi].reset_decoder();
            oViterbi[IndexViterbi].Decode(SplitI[IndexViterbi]+PtrRd, SplitQ[IndexViterbi] + PtrRd, BatchSize1, ViterbiOutputs[IndexViterbi]+1, true,   1, 1, m[3]);
            
            lk.lock();
            for (int k = 0; k < 4; ++k)
                VitSyncResults[IndexViterbi].Metrics[k] = m[k];
            VitSyncResults[IndexViterbi].Available = true;
            DemodulatorQ[IndexViterbi].AdvanceRead();
            lk.unlock();

        }
        else
        {
            ViterbiParameters vp = ViterbiParams;
            lk.unlock();
            
            oViterbi[IndexViterbi].Decode(SplitI[IndexViterbi]+PtrRd, SplitQ[IndexViterbi] + PtrRd, BatchSize1, ViterbiOutputs[IndexViterbi]+1, vp.ExchangeIQ, vp.SignI, vp.SignQ, VitSyncResults[IndexViterbi].Metrics[0]);
            
            lk.lock();
            CvViterbis2VitManager[IndexViterbi].wait(lk, [&] {
                return StopAll || DecodedQ[IndexViterbi].AvailableWrite();
            });
            if (StopAll)
                break;
            int PtrWr = BatchSize1 * static_cast<int>(DecodedQ[IndexViterbi].GetPtrWr());
            lk.unlock();
            
            DifferentialDecode(ViterbiOutputs[IndexViterbi],DiffDec[IndexViterbi]+PtrWr, BatchSize1, Last);
            
            lk.lock();
            DemodulatorQ[IndexViterbi].AdvanceRead();
            DecodedQ[IndexViterbi].AdvanceWrite();
            lk.unlock();
        }
        CvViterbis2VitManager[IndexViterbi].notify_one();

    }
}

/**
 * @brief Fan-out symbols to decoders, four-way metric lock, merge, PRBS/descramble, output queue, optional constellation.
 */
void Receiver::OperateViterbiManager(void)
{
    NumBitsAll = 0;
    NumErrorsAll = 0;
    unsigned int PRBSSeed = 0;
    RawNumBitsAll = 0;
    RawNumErrorsAll = 0;
    #ifdef WRITE_LOG_THR
	mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"Receiver OperateViterbiManager Thread %d\n", gettid());
    fclose(fidthr);
    mtxfilethr.unlock();
    #endif
    #ifdef DEBUG_STATISTICS
                    
    CurrDebugStatistics.MaxMetricsGrowth = 0;
    CurrDebugStatistics.SumMetricsGrowth = 0;
    CurrDebugStatistics.NumBatches = 0;
    #endif
    constexpr int kSymFrame = ReceiverInputBatchIQSymbols;
    alignas(32) float pendingI[kSymFrame];
    alignas(32) float pendingQ[kSymFrame];
    auto lastConstellationDisplay = std::chrono::steady_clock::now();
    uint64_t symbols_total = 0;
    const auto wall_start = std::chrono::steady_clock::now();

    while (!StopAll)
    {
        //4 options for Viterbi - taking also care of Spectrum Inversion
        //(I+jQ)
        //I-jQ
        //j(I+jQ) = -Q+jI
        //j(I-jQ)= Q+jI
        //Controls:
        //  false 1 1
        //false 1 -1
        //true -1 1
        //true 1 1

        if (!phaseTrackingDD_.WaitPopSymbolFrame(pendingI, pendingQ, kSymFrame))
            break;
        symbols_total += static_cast<uint64_t>(kSymFrame);

        // RAW (pre-Viterbi) diagnostics (optional).
        // Compile-time options:
        // - ENABLE_RAW_PREVITERBI_METRICS: enables the whole RAW block below.
        // - DUMP_RAW_INTERCORR: additionally writes the exact TX/RX windows used for correlation to ../data/.
#ifdef ENABLE_RAW_PREVITERBI_METRICS
        if (RxMode != FILE_TX && pTx_ != nullptr)
        {
            // Requested behavior:
            // - While PhaseTrackingDD is NOT locked: for each RX frame processed, pop and discard one TX frame.
            // - On the first RX frame after PhaseTrackingDD becomes locked: pop the next TX frame, dump and
            //   compute inter-correlation on those exact RX/TX symbols. If inter-corr fails, exit(-1).
            static bool phaseWasLocked = false;
            static bool rawCorrDone = false;
            static int rawSkipLockedFrames = 0;

            const bool phaseLocked = phaseTrackingDD_.IsLocked();

            alignas(32) float txFrameI[kSymFrame];
            alignas(32) float txFrameQ[kSymFrame];
            static constexpr int kCorrFrames = 4;
            static constexpr int kCorrLen2 = kCorrFrames * kSymFrame;
            // Search lag over +/- N whole frames (pipeline delay can be several frames).
            static constexpr int kMaxLagFrames = 32;
            static constexpr int kMaxLag2 = kMaxLagFrames * kSymFrame;
            static constexpr int kTxLen2 = kCorrLen2 + 2 * kMaxLag2;
            static constexpr int kTxFramesNeeded = (kTxLen2 + kSymFrame - 1) / kSymFrame;

            static std::vector<float> rxIacc;
            static std::vector<float> rxQacc;
            static std::vector<float> txIacc;
            static std::vector<float> txQacc;
            static std::deque<float> txDelayQueueI;
            static std::deque<float> txDelayQueueQ;
            static double savedPhaseAdj = 0.0;
            
            bool haveEnough = false;

            if (!phaseLocked && !rawCorrDone)
            {
                (void)pTx_->WaitPopTxSymbolFrame(txFrameI, txFrameQ, kSymFrame);
                phaseWasLocked = false;
                rawCorrDone = false;
                rawSkipLockedFrames = 0;
                RawSyncLocked.store(false, std::memory_order_relaxed);
                goto skip_raw;
            }

            // Transition: first iteration after lock -> start accumulation window.
            if (!phaseWasLocked)
            {
                phaseWasLocked = true;
                rawCorrDone = false;
                // Requested: drop one extra RX+TX frame after lock before accumulating/dumping.
                rawSkipLockedFrames = 1;
                RawSyncLocked.store(false, std::memory_order_relaxed);
                rxIacc.clear();
                rxQacc.clear();
                txIacc.clear();
                txQacc.clear();
                rxIacc.reserve(static_cast<size_t>(kCorrLen2));
                rxQacc.reserve(static_cast<size_t>(kCorrLen2));
                txIacc.reserve(static_cast<size_t>(kTxLen2));
                txQacc.reserve(static_cast<size_t>(kTxLen2));
                txDelayQueueI.clear();
                txDelayQueueQ.clear();
            }

            if (rawSkipLockedFrames > 0)
            {
                // Discard one TX frame to keep pairing, and also discard this RX frame.
                (void)pTx_->WaitPopTxSymbolFrame(txFrameI, txFrameQ, kSymFrame);
                rawSkipLockedFrames--;
                goto skip_raw;
            }

            // After the one-shot raw correlation/dump, we must keep draining TX frames; otherwise the TX
            // reference ring fills up and backpressures/stalls the transmitter thread.
            // We now maintain a continuous delay line to calculate SER/rawBER over the whole simulation.
            if (rawCorrDone)
            {
                if (pTx_->TryPopTxSymbolFrame(txFrameI, txFrameQ, kSymFrame))
                {
                    txDelayQueueI.insert(txDelayQueueI.end(), txFrameI, txFrameI + kSymFrame);
                    txDelayQueueQ.insert(txDelayQueueQ.end(), txFrameQ, txFrameQ + kSymFrame);

                    // Wait until we have enough TX symbols to cover the current RX frame.
                    if (txDelayQueueI.size() >= static_cast<size_t>(kSymFrame))
                    {
                        const float c = static_cast<float>(std::cos(-savedPhaseAdj));
                        const float s = static_cast<float>(std::sin(-savedPhaseAdj));

                        uint64_t symErrCount = 0;
                        uint64_t bitErrCount = 0;

                        for (int k = 0; k < kSymFrame; ++k)
                        {
                            const float rI0 = pendingI[k];
                            const float rQ0 = pendingQ[k];
                            const float rIr = rI0 * c + rQ0 * s;
                            const float rQr = rQ0 * c - rI0 * s;

                            const float rxIh = (rIr >= 0.0f) ? 1.0f : -1.0f;
                            const float rxQh = (rQr >= 0.0f) ? 1.0f : -1.0f;
                            
                            const float txIh = (txDelayQueueI[k] >= 0.0f) ? 1.0f : -1.0f;
                            const float txQh = (txDelayQueueQ[k] >= 0.0f) ? 1.0f : -1.0f;

                            // Apply the same quadrant ambiguity resolution to TX symbol (actually applied to rx, so rx should match tx)
                            // But wait! `phaseAdj` already includes both the continuous phase offset AND the k*pi/2 quadrant rotation.
                            // So `rIr` and `rQr` are ALREADY rotated to the correct quadrant.
                            // Thus, `rxIh` and `rxQh` should EXACTLY match `txIh` and `txQh`!
                            
                            // Let's check: in the one-shot lock, the phase is found such that rx * exp(-j*phase) matches tx.
                            // phaseAdj = phase + quad * pi/2.
                            // So we just need to do rx_rot = rx * exp(-j*phaseAdj).
                            // Which is exactly what is done above!
                            const bool symErr = (rxIh != txIh) || (rxQh != txQh);
                            symErrCount += static_cast<uint64_t>(symErr);
                            bitErrCount += static_cast<uint64_t>(rxIh != txIh) + static_cast<uint64_t>(rxQh != txQh);
                        }

                        RawNumSymsAll += static_cast<uint64_t>(kSymFrame);
                        RawNumSymErrorsAll += symErrCount;
                        RawNumBitsAll += static_cast<uint64_t>(2 * kSymFrame);
                        RawNumErrorsAll += bitErrCount;

                        txDelayQueueI.erase(txDelayQueueI.begin(), txDelayQueueI.begin() + kSymFrame);
                        txDelayQueueQ.erase(txDelayQueueQ.begin(), txDelayQueueQ.begin() + kSymFrame);
                    }
                }
                goto skip_raw;
            }

            // Accumulate RX symbols (full frames) and matching TX frames.
            // New logic: keep popping both until we have enough for txIacc.
            // Since we need them to be 1:1 aligned, we push to rxIacc until it has the same size as txIacc.
            if (static_cast<int>(txIacc.size()) < kTxLen2)
            {
                for (int k = 0; k < kSymFrame; ++k)
                {
                    rxIacc.push_back(pendingI[k]);
                    rxQacc.push_back(pendingQ[k]);
                }

                if (pTx_->WaitPopTxSymbolFrame(txFrameI, txFrameQ, kSymFrame))
                {
                    for (int k = 0; k < kSymFrame && static_cast<int>(txIacc.size()) < kTxLen2; ++k)
                    {
                        txIacc.push_back(txFrameI[k]);
                        txQacc.push_back(txFrameQ[k]);
                    }
                }
                else
                {
                    break; // Error
                }
            }

            haveEnough = (static_cast<int>(txIacc.size()) >= kTxLen2);
            if (!haveEnough)
            {
                goto skip_raw;
            }

            if (haveEnough)
            {
                // Aliases for clarity below.
                const float* rxI2 = rxIacc.data();
                const float* rxQ2 = rxQacc.data();
                const float* txI2 = txIacc.data();
                const float* txQ2 = txQacc.data();

#ifdef DUMP_RAW_INTERCORR
            mkdir("../data", 0755);
            {
                FILE* ftx = std::fopen("../data/raw_icorr_tx.bin", "wb");
                FILE* frx = std::fopen("../data/raw_icorr_rx.bin", "wb");
                FILE* fmeta = std::fopen("../data/raw_icorr_meta.txt", "wt");
                if (fmeta)
                {
                    std::fprintf(fmeta, "kCorrLen=%d\nkMaxLag=%d\nnTx=%d\nphaseLocked=1\n",
                                 kCorrLen2, kMaxLag2, kTxLen2);
                    std::fclose(fmeta);
                }
                if (ftx)
                {
                    for (int k = 0; k < kTxLen2; ++k)
                    {
                        float iq[2] = {txI2[k], txQ2[k]};
                        std::fwrite(iq, sizeof(float), 2, ftx);
                    }
                    std::fclose(ftx);
                }
                if (frx)
                {
                    for (int k = 0; k < kCorrLen2; ++k)
                    {
                        float iq[2] = {rxI2[k], rxQ2[k]};
                        std::fwrite(iq, sizeof(float), 2, frx);
                    }
                    std::fclose(frx);
                }
            }
#endif

            double es = 0.0;
            for (int k = 0; k < kCorrLen2; ++k)
                es += static_cast<double>(txI2[kMaxLag2 + k] * txI2[kMaxLag2 + k] + txQ2[kMaxLag2 + k] * txQ2[kMaxLag2 + k]);
            es /= static_cast<double>(kCorrLen2);

            double bestMag = -1.0, bestRe = 0.0, bestIm = 0.0;
            int bestLag = 0;
            double sumMag = 0.0;
            int nMag = 0;
            for (int lag = -kMaxLag2; lag <= kMaxLag2; ++lag)
            {
                const int txBase = kMaxLag2 + lag;
                double re = 0.0, im = 0.0;
                for (int k = 0; k < kCorrLen2; ++k)
                {
                    const double ar = rxI2[k];
                    const double ai = rxQ2[k];
                    const double br = txI2[txBase + k];
                    const double bi = txQ2[txBase + k];
                    re += ar * br + ai * bi;
                    im += ai * br - ar * bi;
                }
                const double mag = std::sqrt(re * re + im * im);
                sumMag += mag;
                nMag++;
                if (mag > bestMag)
                {
                    bestMag = mag;
                    bestRe = re;
                    bestIm = im;
                    bestLag = lag;
                }
            }

            const double meanMag = (nMag > 0) ? (sumMag / static_cast<double>(nMag)) : 0.0;
            const double peakToMean = (meanMag > 1e-30) ? (bestMag / meanMag) : 0.0;
            constexpr double kPeakToMeanThr = 50.0;

            // For display: BestAbs=peak/mean, ThrAbs=threshold, PeakAbs=absolute peak.
            RawSyncPeakAbs.store(bestMag, std::memory_order_relaxed);
            RawSyncBestAbs.store(peakToMean, std::memory_order_relaxed);
            RawSyncThrAbs.store(kPeakToMeanThr, std::memory_order_relaxed);
            RawSyncLagSym.store(bestLag, std::memory_order_relaxed);
            RawSyncPeakPhaseRad.store(std::atan2(bestIm, bestRe), std::memory_order_relaxed);

            if (!(peakToMean > kPeakToMeanThr))
            {
                // Reset state and skip so we try again next time PhaseDD is locked
                goto skip_raw;
            }

            RawSyncLocked.store(true, std::memory_order_relaxed);
            rawCorrDone = true;

            std::cout << "[RAW] \033[32mLOCKED\033[0m"
                      << " peakToMean=" << peakToMean
                      << " thr=" << kPeakToMeanThr
                      << " peakAbs=" << bestMag
                      << " phaseRad=" << std::atan2(bestIm, bestRe)
                      << " lagSym=" << bestLag
                      << std::endl;

            // Compute SER + raw BER on the aligned window using hard decisions.
            const double phase = std::atan2(bestIm, bestRe);
            const int txBase = kMaxLag2 + bestLag;

            // corr(lag) = sum rx * conj(tx) gives phase = phi_rx - phi_tx, but QPSK has a quadrant ambiguity.
            // Resolve ambiguity by trying phase + k*pi/2 and picking the smallest SER.
            int bestQuad = 0;
            uint64_t bestSymErr = UINT64_MAX;
            uint64_t bestBitErr = UINT64_MAX;
            for (int quad = 0; quad < 4; ++quad)
            {
                const double phaseAdj = phase + (static_cast<double>(quad) * (M_PI / 2.0));
                const float c = static_cast<float>(std::cos(-phaseAdj));
                const float s = static_cast<float>(std::sin(-phaseAdj));

                uint64_t symErrCount = 0;
                uint64_t bitErrCount = 0;
                for (int k = 0; k < kCorrLen2; ++k)
                {
                    const float rI0 = rxI2[k];
                    const float rQ0 = rxQ2[k];
                    // (rI + j rQ) * exp(-j*phaseAdj)
                    const float rIr = rI0 * c + rQ0 * s;
                    const float rQr = rQ0 * c - rI0 * s;

                    const float rxIh = (rIr >= 0.0f) ? 1.0f : -1.0f;
                    const float rxQh = (rQr >= 0.0f) ? 1.0f : -1.0f;
                    const float txIh = (txI2[txBase + k] >= 0.0f) ? 1.0f : -1.0f;
                    const float txQh = (txQ2[txBase + k] >= 0.0f) ? 1.0f : -1.0f;

                    const bool symErr = (rxIh != txIh) || (rxQh != txQh);
                    symErrCount += static_cast<uint64_t>(symErr);
                    bitErrCount += static_cast<uint64_t>(rxIh != txIh) + static_cast<uint64_t>(rxQh != txQh);
                }

                if (symErrCount < bestSymErr || (symErrCount == bestSymErr && bitErrCount < bestBitErr))
                {
                    bestSymErr = symErrCount;
                    bestBitErr = bitErrCount;
                    bestQuad = quad;
                }
            }

            const double phaseAdj = phase + (static_cast<double>(bestQuad) * (M_PI / 2.0));
            RawSyncAppliedPhaseRad.store(phaseAdj, std::memory_order_relaxed);
            RawSyncAppliedQuad.store(bestQuad, std::memory_order_relaxed);

            RawNumSymsAll += static_cast<uint64_t>(kCorrLen2);
            RawNumSymErrorsAll += bestSymErr;
            RawNumBitsAll += static_cast<uint64_t>(2 * kCorrLen2);
            RawNumErrorsAll += bestBitErr;

            // Prepare continuous delay queue with the unused 'future' TX history 
            savedPhaseAdj = phaseAdj;
            txDelayQueueI.clear();
            txDelayQueueQ.clear();
            
            const int nextRxIdx = kTxLen2;
            const int neededTxIdx = nextRxIdx + txBase;
            
            if (neededTxIdx < kTxLen2)
            {
                // We have the needed TX symbol in txIacc.
                txDelayQueueI.insert(txDelayQueueI.end(), txIacc.begin() + neededTxIdx, txIacc.end());
                txDelayQueueQ.insert(txDelayQueueQ.end(), txQacc.begin() + neededTxIdx, txQacc.end());
            }
            else
            {
                // We need to fetch more TX symbols to catch up.
                int txToFetch = neededTxIdx - kTxLen2;
                while (txToFetch > 0)
                {
                    int fetchNow = std::min(txToFetch, kSymFrame);
                    if (pTx_->WaitPopTxSymbolFrame(txFrameI, txFrameQ, kSymFrame))
                    {
                        if (fetchNow < kSymFrame)
                        {
                            // Keep the remaining part in the delay queue
                            txDelayQueueI.insert(txDelayQueueI.end(), txFrameI + fetchNow, txFrameI + kSymFrame);
                            txDelayQueueQ.insert(txDelayQueueQ.end(), txFrameQ + fetchNow, txFrameQ + kSymFrame);
                        }
                        txToFetch -= fetchNow;
                    }
                    else
                    {
                        break;
                    }
                }
            }

            rawCorrDone = true;
        }
    }
skip_raw:;
#endif // ENABLE_RAW_PREVITERBI_METRICS

        if (constellationDisplay_ && constellationCfg_.PeriodSec > 0.0)
        {
            const auto now = std::chrono::steady_clock::now();
            if (now - lastConstellationDisplay >= std::chrono::duration<double>(constellationCfg_.PeriodSec))
            {
                const double t_rate = static_cast<double>(symbols_total) / SymbolRate;
                const double t_sim = std::chrono::duration<double>(now - wall_start).count();
                // Extra status overlay for the constellation window (space-separated key=value).
                // Keep it short to avoid slowing down the pipe/GUI.
                const double symRateHz = SymbolRateEstimateHz.load(std::memory_order_relaxed);
                const double symRateMsps = symRateHz / 1e6;
                // Symbol-rate estimator error in ppm, referenced to NOMINAL symbol rate.
                const double symRatePpm = (symRateHz > 0.0) ? ((symRateHz - SymbolRate) * 1.0e6 / SymbolRate) : 0.0;
                // Coarse/fixed ppm correction from SymbolRateEstimator (requested).
                // Channel totalPpm affects symbol rate with opposite sign in our estimator convention, so we negate here.
                const double symRateCorrPpm = -symRatePpm;

                // Gardner ppm residual estimate from omega deviation around nominal (fine correction).
                const double omega = timingTracking_.GetGardnerOmega();
                const double omegaNom = timingTracking_.GetGardnerOmegaNom();
                // From omega ≈ Fs/Rs and omegaNom = Fs/Rs_nom => Rs/Rs_nom ≈ omegaNom/omega.
                // totalPpm ≈ (omegaNom/omega - 1) * 1e6
                const double gardnerPpm = (omega > 1e-20 && omegaNom > 1e-20) ? ((omegaNom / omega) - 1.0) * 1.0e6 : 0.0;

                // Requested total corrected ppm: coarse (SymbolRateEstimator) + fine (Gardner).
                const double rxTotalPpmCorr = symRateCorrPpm + gardnerPpm;
                const int gardLocked = timingTracking_.IsLocked() ? 1 : 0;
                const int phaseLocked = phaseTrackingDD_.IsLocked() ? 1 : 0;
                const int vitLocked = ViterbiSynchronized ? 1 : 0;
                const int prbsLocked = PRBSSynchronized ? 1 : 0;

                const int centralReady = freqCorrector_.IsFreqReady() ? 1 : 0;
                const double centralTargetHz = freqCorrector_.GetTargetHz();
                const double centralNcoHz = freqCorrector_.GetNcoHz();
                const double gainDb = freqCorrector_.GetGainDb();
                const double phaseFreqHz = phaseTrackingDD_.GetLastFreqEstHz();

                const unsigned long long numBits = static_cast<unsigned long long>(NumBitsAll);
                const unsigned long long numErrors = static_cast<unsigned long long>(NumErrorsAll);
                const double ber = (numBits > 0ULL) ? (static_cast<double>(numErrors) / static_cast<double>(numBits)) : 0.0;

                const double samplerMsps = pSampler ? pSampler->GetLastThroughputMsps() : 0.0;
                const double rxFilterMsps = lastRxFilterThroughputMsps_.load(std::memory_order_relaxed);

                // Channel applied values (if available).
                double chanFreqHz = 0.0;
                double chanTotalPpm = 0.0;
                double chanGainDb = 0.0;
                int chanSeg = 0;
                double chanTRate = 0.0;
                double chanSnrDb = std::numeric_limits<double>::quiet_NaN();
                if (pChannel_)
                {
                    const auto s = pChannel_->GetCurrentApplied();
                    chanFreqHz = s.FreqHz;
                    chanTotalPpm = s.TotalPpm;
                    chanGainDb = s.GainDb;
                    chanSeg = s.Segment;
                    chanTRate = s.TRate;
                    chanSnrDb = s.SnrDb;
                }

                double evmRms = phaseTrackingDD_.GetLastEvmRms();
                double snrFromEvmDb = (evmRms > 1e-12) ? (-20.0 * std::log10(evmRms)) : std::numeric_limits<double>::quiet_NaN();

                char extra[768];
                std::snprintf(extra, sizeof(extra),
                              "symRateMsps=%.6f symRatePpm=%.3f gardLocked=%d phaseLocked=%d vitLocked=%d prbsLocked=%d "
                              "symRateCorrPpm=%.3f gardnerPpm=%.3f rxTotalPpmCorr=%.3f "
                              "numBits=%llu numErrors=%llu ber=%.6e "
                              "samplerMsps=%.3f rxFilterMsps=%.3f "
                              "centralReady=%d centralTargetHz=%.3f centralNcoHz=%.3f gainDb=%.3f phaseFreqHz=%.3f "
                              "chanFreqHz=%.3f chanTotalPpm=%.6f chanGainDb=%.3f chanSeg=%d chanTRate=%.3f "
                              "chanSnrDb=%.2f evmRms=%.6f snrFromEvmDb=%.2f",
                              symRateMsps, symRatePpm, gardLocked, phaseLocked, vitLocked, prbsLocked,
                              symRateCorrPpm, gardnerPpm, rxTotalPpmCorr,
                              numBits, numErrors, ber,
                              samplerMsps, rxFilterMsps,
                              centralReady, centralTargetHz, centralNcoHz, gainDb, phaseFreqHz,
                              chanFreqHz, chanTotalPpm, chanGainDb, chanSeg, chanTRate, chanSnrDb, evmRms,
                              snrFromEvmDb);

                constellationDisplay_->UpdateEx(pendingI, pendingQ, kSymFrame, t_sim, t_rate, std::string(extra));
                lastConstellationDisplay = now;
            }
        }

        {
                for(int i = 0; i < 3; i++)
                {
                    std::unique_lock<std::mutex> lk(mtxQueues);
                    CvViterbis2VitManager[i].wait(lk, [&] {
                        return StopAll || DemodulatorQ[i].AvailableWrite();
                    });
                    if(StopAll)
                        break;
                }
                if(StopAll)
                    break;

#ifdef DEBUG_GARDNER_OUTPUTS
                {
                    static FILE* gVitFromGardnerDump = nullptr;
                    if (!gVitFromGardnerDump)
                    {
                        mkdir("../data", 0755);
                        gVitFromGardnerDump = std::fopen("../data/viterbi_input_from_gardner.bin", "wb");
                    }
                    if (gVitFromGardnerDump)
                    {
                        for (int i = 0; i < kSymFrame; ++i)
                        {
                            float iq[2] = {pendingI[i], pendingQ[i]};
                            std::fwrite(iq, sizeof(float), 2, gVitFromGardnerDump);
                        }
                        std::fflush(gVitFromGardnerDump);
                    }
                }
#endif

                // Option A: keep the mutex during access to the slot (SplitI/SplitQ) + AdvanceWrite.
                std::unique_lock<std::mutex> lkSlots(mtxQueues);
                int PtrWr = static_cast<int>(DemodulatorQ[0].GetPtrWr());
                Split3(pendingI, pendingQ, kSymFrame, PtrWr*BatchSize1);
                if(StopAll)
                    break;
                if(ViterbiSynchronized == false)
                {
                    for(int i = 0; i < 3; i++)
                    {
                        VitSyncResults[i].Available = false;
                    }
                }
                for(int i = 0; i < 3; i++)
                {
                    DemodulatorQ[i].AdvanceWrite();
                }
                lkSlots.unlock();
                CvVitManager2Vit.notify_all();

                if(StopAll)
                    break;
                if(!ViterbiSynchronized)
                {
                    static bool prevViterbiLocked = false;
                    static auto lastViterbiStatus = std::chrono::steady_clock::now();
                    for(int i = 0; i < 3; i++)
                    {
                        std::unique_lock<std::mutex> lk(mtxQueues);
                        CvViterbis2VitManager[i].wait(lk, [&] {
                            return StopAll || VitSyncResults[i].Available;
                        });
                        if(StopAll)
                            break;
                    }
                    if(StopAll)
                        break;
                    //Test if Synchronization
                    int BestIndex[3];
                    __m128 mThresh = _mm_set1_ps(ViterbiThreshold1);
                    bool Success = true;
                    float minMargin[3] = {0.0f, 0.0f, 0.0f};
                    float bestMetrics[3] = {0.0f, 0.0f, 0.0f};
                    for(int i = 0; i < 3; i++)
                    {
                        BestIndex[i] = 0;
                        float m0 = VitSyncResults[i].Metrics[0];
                        float m1 = VitSyncResults[i].Metrics[1];
                        float m2 = VitSyncResults[i].Metrics[2];
                        float m3 = VitSyncResults[i].Metrics[3];
                        float BestMetrics = m0;
                        for(int j = 1; j < 4; j++)
                        {
                            if(VitSyncResults[i].Metrics[j] < BestMetrics)
                            {
                                BestIndex[i] = j;
                                BestMetrics = VitSyncResults[i].Metrics[j];
                            }
                        }
                        bestMetrics[i] = BestMetrics;
                        // Compute smallest margin to competing hypotheses: min_j!=best (Mj - Mbest).
                        // This is what the threshold check below enforces.
                        const float best = BestMetrics;
                        float margins[4] = {m0 - best, m1 - best, m2 - best, m3 - best};
                        margins[BestIndex[i]] = 1e7f;
                        minMargin[i] = std::min(std::min(margins[0], margins[1]), std::min(margins[2], margins[3]));

                        VitSyncResults[i].Metrics[BestIndex[i]] = 1e7;
                        __m128 mResults = _mm_loadu_ps(VitSyncResults[i].Metrics);
                        __m128 mBest = _mm_set1_ps(BestMetrics);
                        mResults = _mm_sub_ps(mResults,mBest );
                        __m128 mCmp = _mm_cmp_ps(mResults, mThresh, _CMP_GE_OS);
                        int PassThresh =  _mm_movemask_ps(mCmp);
                        if(PassThresh != 15)
                        {
                            Success = false;
                            break;
                        }
                    }

                    if(Success)
                    {
                        if((BestIndex[0]== BestIndex[1]) && (BestIndex[0]== BestIndex[2]))
                        {
                            switch (BestIndex[0]) {
                            case 0:
                                ViterbiParams.ExchangeIQ = false;
                                ViterbiParams.SignI = 1;
                                ViterbiParams.SignQ = 1;
                                break;
                            case 1:
                                ViterbiParams.ExchangeIQ = false;
                                ViterbiParams.SignI = 1;
                                ViterbiParams.SignQ = -1;
                                break;
                            case 2:
                                ViterbiParams.ExchangeIQ = true;
                                ViterbiParams.SignI = -1;
                                ViterbiParams.SignQ = 1;
                                break;
                            case 3:
                                ViterbiParams.ExchangeIQ = true;
                                ViterbiParams.SignI = 1;
                                ViterbiParams.SignQ = 1;
                                break;   
                            }
                            ViterbiSynchronized = true;
                            std::cout << "[ViterbiSync] \033[32mLOCKED\033[0m"
                                      << " thresh=" << ViterbiThreshold1
                                      << " minMargin=(" << minMargin[0] << ", " << minMargin[1] << ", " << minMargin[2] << ")"
                                      << " bestIndex=(" << BestIndex[0] << ", " << BestIndex[1] << ", " << BestIndex[2] << ")"
                                      << std::endl;
                            prevViterbiLocked = true;
                            for(int i = 0; i < 3; i++)
                                oViterbi[i].reset_decoder();
                        }
                        else
                        {
                            // Thresholds passed but hypotheses disagree.
                            prevViterbiLocked = false;
                        }
                    }
                    else
                    {
                        // Failed threshold: at least one block has insufficient margin vs competitors.
                        prevViterbiLocked = false;
                    }

                    // Periodic status (even before first lock) to explain why it's not locked.
                    {
                        const auto now = std::chrono::steady_clock::now();
                        if (displayPeriodSec_ > 0.0 &&
                            now - lastViterbiStatus >= std::chrono::duration<double>(displayPeriodSec_))
                        {
                            const char* reason = "unknown";
                            if (!Success)
                                reason = "margin_below_threshold";
                            else if (!((BestIndex[0] == BestIndex[1]) && (BestIndex[0] == BestIndex[2])))
                                reason = "hypothesis_mismatch";
                            else
                                reason = "locked_pending";

                            std::cout << "[ViterbiSync] status"
                                      << " locked=" << (ViterbiSynchronized ? 1 : 0)
                                      << " reason=" << reason
                                      << " thresh=" << ViterbiThreshold1
                                      << " minMargin=(" << minMargin[0] << ", " << minMargin[1] << ", " << minMargin[2] << ")"
                                      << " bestIndex=(" << BestIndex[0] << ", " << BestIndex[1] << ", " << BestIndex[2] << ")"
                                      << " bestMetrics=(" << bestMetrics[0] << ", " << bestMetrics[1] << ", " << bestMetrics[2] << ")"
                                      << std::endl;
                            lastViterbiStatus = now;
                        }
                    }
                }
                else 
                {
                    for(int i = 0; i < 3; i++)
                    {
                        std::unique_lock<std::mutex> lk(mtxQueues);
                        CvViterbis2VitManager[i].wait(lk, [&] {
                            return StopAll || DecodedQ[i].AvailableRead();
                        });
                        if(StopAll)
                            break;
                    }

                    #ifdef DEBUG_STATISTICS
                    for(int i = 0; i < 3; i++)
                    {
                        if(CurrDebugStatistics.MaxMetricsGrowth < VitSyncResults[i].Metrics[0])
                            CurrDebugStatistics.MaxMetricsGrowth = VitSyncResults[i].Metrics[0];
                        CurrDebugStatistics.SumMetricsGrowth += VitSyncResults[i].Metrics[0];
                    }
                    CurrDebugStatistics.NumBatches += 3;
                    CurrDebugStatistics.MeanMetricsGrowth =  CurrDebugStatistics.SumMetricsGrowth/double(CurrDebugStatistics.NumBatches);
                    
                    double currentGrowth = (VitSyncResults[0].Metrics[0] + VitSyncResults[1].Metrics[0] + VitSyncResults[2].Metrics[0]) / 3.0;
                    if(CurrDebugStatistics.NumBatches <= 3) {
                        CurrDebugStatistics.EmaMetricsGrowth = currentGrowth;
                    } else {
                        CurrDebugStatistics.EmaMetricsGrowth = 0.99 * CurrDebugStatistics.EmaMetricsGrowth + 0.01 * currentGrowth;
                    }
                    
                    // Viterbi Unlock Condition
                    constexpr double kViterbiDesyncGrowthThr = 40.0; // Calibrated for Viterbi cliff (SNR ~1.9dB, BER jumping to >2e-1)
                    if (CurrDebugStatistics.EmaMetricsGrowth > kViterbiDesyncGrowthThr && CurrDebugStatistics.NumBatches > 100)
                    {
                        ViterbiSynchronized = false;
                        PRBSSynchronized = false;
                        NumBitsAll = 0;
                        NumErrorsAll = 0;
                        
                        std::cout << "[ViterbiSync] \033[31mDESYNC\033[0m EmaMetricsGrowth=" << CurrDebugStatistics.EmaMetricsGrowth << " > " << kViterbiDesyncGrowthThr << std::endl;
                        // Reset stats for next sync
                        CurrDebugStatistics.NumBatches = 0;
                        CurrDebugStatistics.SumMetricsGrowth = 0;
                        CurrDebugStatistics.MaxMetricsGrowth = 0;
                    }
                    #endif

                    if(StopAll)
                        break;
                    // Option A: keep the mutex during reading DiffDec[slot] + AdvanceRead.
                    std::unique_lock<std::mutex> lkDec(mtxQueues);
                    int PtrRd = BatchSize1*static_cast<int>(DecodedQ[0].GetPtrRd());
                    int PtrOut = 0;
                    for(int i = 0; i < BatchSize1; i++)
                    {
                        Merged[PtrOut++] = DiffDec[0][PtrRd+i];
                        Merged[PtrOut++] = DiffDec[1][PtrRd+i];
                        Merged[PtrOut++] = DiffDec[2][PtrRd+i];
                    }
                    for(int i = 0; i < 3; i++)
                        DecodedQ[i].AdvanceRead();
                    lkDec.unlock();
                    for(int i = 0; i < 3; i++)
                        CvViterbis2VitManager[i].notify_one();

                    //Descramble
                    std::unique_lock<std::mutex> lkOut(mtxQueues);
                    CvOut2Rx.wait(lkOut, [&] { return StopAll || OutputQ.AvailableWrite(); });
                    if(StopAll)
                        break;
                    const int outSlotByte = BatchSize3 * static_cast<int>(OutputQ.GetPtrWr());
                    objDescrambler.Descramble(Merged, OutputAll + outSlotByte, BatchSize3);
                    OutputQ.AdvanceWrite();
                    lkOut.unlock();

                    if(RxMode == FILE_TX)
                    {
                        CvRx2Out.notify_one();
                    }
                    else {
                        // Option A: keep mtxQueues for all access to the slot OutputAll[outSlotByte..]
                        std::unique_lock<std::mutex> lkOutRead(mtxQueues);
                        int PtrRdOut = outSlotByte;
                        if (!ViterbiSynchronized)
                        {
                            // Viterbi just lost lock on this batch, or we are otherwise not locked.
                            // Don't try to sync PRBS on garbage data.
                        }
                        else if(!PRBSSynchronized)
                        {
                            int PtrStart;
                            int NumErrorsAtLock = 0;
                            // Debug: corrupt the descrambled PRBS stream before lock detection.
                            // When PRBSInjectStride > 0, PRBS sync should fail (stay unlocked).
                            InjectPRBSErrorDeterministic(OutputAll + PtrRdOut, BatchSize3);
                            PRBSSynchronized = SyncPRBS(OutputAll+PtrRdOut, PtrStart, PRBSSeed, NumErrorsAtLock);
                            if(PRBSSynchronized)
                            {
                                std::cout << "[PRBS] \033[32mLOCKED\033[0m"
                                          << " NumErrors=" << NumErrorsAtLock
                                          << " thresh=" << PRBSThreshold
                                          << " PtrStart=" << PtrStart
                                          << " Seed=" << PRBSSeed
                                          << std::endl;
                                oPrbs.CreateOutputs(PRBSSeed, BatchSize3 - PtrStart, PrbsOut);
                                CountErrors(OutputAll+PtrRdOut+PtrStart, PrbsOut, BatchSize3-PtrStart);
                            }
                        }
                        else {
                            oPrbs.CreateOutputs(PRBSSeed, BatchSize3, PrbsOut);
                            uint64_t prevErrors = NumErrorsAll;
                            uint64_t prevBits = NumBitsAll;
                            CountErrors(OutputAll+PtrRdOut, PrbsOut, BatchSize3);
                            uint64_t diffErrors = NumErrorsAll - prevErrors;
                            uint64_t diffBits = NumBitsAll - prevBits;
                            if (diffBits > 0 && (static_cast<double>(diffErrors) / static_cast<double>(diffBits)) > 0.3) {
                                PRBSSynchronized = false;
                                std::cout << "[PRBS] \033[31mDESYNC\033[0m High BER detected (" << diffErrors << "/" << diffBits << ")" << std::endl;
                            }
                        }
                        OutputQ.AdvanceRead();
                        lkOutRead.unlock();
                        CvOut2Rx.notify_one();
                    }
                }
        }
    }

}
void  Receiver:: Split3(float *InputI, float *InputQ, int Length, int PtrWr)
{
    __m256i mIndex[3];
    mIndex[0] = _mm256_setr_epi32(0, 3, 6, 9, 12, 15, 18, 21);
    mIndex[1] = _mm256_setr_epi32(1, 4, 7, 10, 13, 16, 19, 22);
    mIndex[2] = _mm256_setr_epi32(2, 5, 8, 11, 14, 17, 20, 23);
    //mIndex[3] = _mm256_set1_epi32(24);
    int PtrOut = PtrWr;
    for(int i = 0; i<Length; i += 24)
    {
        __m256 I0 = _mm256_i32gather_ps(InputI + i,mIndex[0],4);
        __m256 Q0 = _mm256_i32gather_ps(InputQ + i,mIndex[0],4);
        __m256 I1 = _mm256_i32gather_ps(InputI + i,mIndex[1],4);
        __m256 Q1 = _mm256_i32gather_ps(InputQ + i,mIndex[1],4);
        __m256 I2 = _mm256_i32gather_ps(InputI + i,mIndex[2],4);
        __m256 Q2 = _mm256_i32gather_ps(InputQ + i,mIndex[2],4);
        _mm256_store_ps(SplitI[0]+PtrOut,I0);
        _mm256_store_ps(SplitQ[0]+PtrOut,Q0);
        _mm256_store_ps(SplitI[1]+PtrOut,I1);
        _mm256_store_ps(SplitQ[1]+PtrOut,Q1);
        _mm256_store_ps(SplitI[2]+PtrOut,I2);
        _mm256_store_ps(SplitQ[2]+PtrOut,Q2);
        PtrOut += 8;
    }


}
void Receiver::shorts_to_floats_avx2(__m256i v16, __m256* out0, __m256* out1)
{
    // Extract lower 128 bits (8 shorts)
   __m128i lo = _mm256_castsi256_si128(v16);

    // Extract upper 128 bits (8 shorts)
    __m128i hi = _mm256_extracti128_si256(v16, 1);

    // Sign-extend int16 -> int32
    __m256i lo32 = _mm256_cvtepi16_epi32(lo); //I0Q0I1Q1I2Q2I3Q3
    __m256i hi32 = _mm256_cvtepi16_epi32(hi);//I4Q4I5Q5I6Q6I7Q7
    lo32 = _mm256_shuffle_epi32(lo32,0xD8);//11011000 I0I1Q0Q1I2I3Q2Q3
    hi32 = _mm256_shuffle_epi32(hi32,0xD8);//I45,Q45,I67,Q67
     
    __m256i mi32 = _mm256_unpacklo_epi64 (lo32,hi32);//I01,I45,I23,I67 
    __m256i mq32 = _mm256_unpackhi_epi64 (lo32,hi32);//Q01,Q45,Q23,Q67
    mi32 = _mm256_permute4x64_epi64(mi32,0xD8);//11011000 
    mq32 = _mm256_permute4x64_epi64(mq32,0xD8);//11011000 

    // Convert int32 -> float
    *out0 = _mm256_cvtepi32_ps(mi32);
    *out1 = _mm256_cvtepi32_ps(mq32);
}

void Receiver::TakeEvenDebug(float *x, float *Output, int Length)
{
    //for(int i = 0; i<16;i++)
    //    x[i] = i;
    for (size_t i = 0; i < Length; i += 16) 
    {
    __m256 a0 = _mm256_loadu_ps(x + i);
    __m256 a1 = _mm256_loadu_ps(x + i + 8);

    __m256 e0 = _mm256_shuffle_ps(a0, a1, _MM_SHUFFLE(2,0,2,0));
    __m256d e0d = _mm256_castps_pd(e0);
    e0d = _mm256_permute4x64_pd(e0d, 0xD8);//1000
    e0 = _mm256_castpd_ps(e0d);


    _mm256_storeu_ps(Output + i/2, e0);
    }

}
void Receiver::CountErrors(unsigned char *Input, unsigned char *Template, int Length)
{
    __m256i mOnes = _mm256_set1_epi8(1);
    for(int j = 0; j < Length; j+=32)
    {
        __m256i mA = _mm256_loadu_si256((__m256i*) (Input+j));
        __m256i mB = _mm256_loadu_si256((__m256i*) (Template+j));
        mA = _mm256_xor_si256(mA, mB);
        mA = _mm256_cmpeq_epi8(mA,mOnes);
        int cmp = _mm256_movemask_epi8(mA);
        int NumErrors1 = _mm_popcnt_u32(cmp);
        NumErrorsAll += NumErrors1;
    }
    NumBitsAll += Length;
}


/**
 * @brief Deterministic PRBS error injection used to test PRBS lock behavior.
 *
 * This flips bits directly in the descrambled byte buffer passed to SyncPRBS().
 * Intended for debugging only.
 */
static void InjectPRBSErrorDeterministic(unsigned char* buf, int length)
{
    if (PRBSInjectStride <= 0 || !buf || length <= 0)
        return;

    const int start = std::max(0, PRBSInjectStart);
    for (int i = start; i < length; ++i)
    {
        if (((i - PRBSInjectStart) % PRBSInjectStride) == 0)
            buf[i] ^= 1;
    }
}
bool Receiver::SyncPRBS(unsigned char *In, int &PtrStart,unsigned int &Seed, int& NumErrorsAtLock)
{
    bool Success = false;
    NumErrorsAtLock = 0;

    int Ptr = BatchSize1 - 23;
    while(!Success)
    {
        if((Ptr + 23 + 128) >= BatchSize3)
            break; 
        Seed = 0;
        int i;
        for ( i = 0; i < 23; i++)
        {
            Seed = (Seed << 1) | (unsigned int) In[Ptr+i];
        }
        int Ptr1 = Ptr + i;
        
        oPrbs.CreateOutputs(Seed,128,PrbsOut);

        __m256i mOnes = _mm256_set1_epi8(1);
        int NumErrors = 0;
        for(int j = 0; j < 128; j+=32)
        {
            __m256i mA = _mm256_loadu_si256((__m256i*) (In+Ptr1+j));
            __m256i mB = _mm256_loadu_si256((__m256i*) (PrbsOut+j));
            mA = _mm256_xor_si256(mA, mB);
            mA = _mm256_cmpeq_epi8(mA,mOnes);
            int cmp = _mm256_movemask_epi8(mA);
            int NumErrors1 = _mm_popcnt_u32(cmp);
            NumErrors += NumErrors1;
        }
        if(NumErrors <= PRBSThreshold)
        {
            Success = true;
            NumErrorsAtLock = NumErrors;
            PtrStart = Ptr1 + 128;
        }
        else {
        
            Ptr+= 128;
        }
    }
    return Success;
}
