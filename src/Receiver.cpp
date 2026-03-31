#include "Receiver.h"
#include <chrono>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <sys/stat.h>

extern std::mutex mtxfilethr;

Receiver::Receiver(/* args */):oBufferFilter(SPB*256,32*SPB),OutputQ(LengthQueue)
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
                            const SymbolRateEstimatorConfig& symRateCfg)
{
    StopAll = false;
    RxMode = RxModeIn;
    oBufferFilter.Reset();
    objRxFilter.CreateObjects(RollOff);
    ViterbiSynchronized = false;
    PRBSSynchronized = false;
    for(int i = 0; i < 3; i++)
    {
        DemodulatorQ[i].Reset();
        DecodedQ[i].Reset();
    }
    SymRatePeakToMedianThreshold = symRateCfg.PeakToMedianThreshold;
    SymRateMaxRelativeJump = symRateCfg.MaxRelativeJump;
    SymRateEstimatePeriodSec = std::max(0.01, symRateCfg.EstimatePeriodSec);
    objSymRateEstimator.Reset(SamplingFrequency, symRateCfg.FftSize, symRateCfg.MaxOffsetHz);
    SymRateWindowBatches = objSymRateEstimator.GetRequiredBatches(SPB);
    SymRateBatchCounter = 0;
    SymRateAcceptedCount = 0;
    SymbolRateDetected.store(false, std::memory_order_relaxed);
    SymbolRateEstimateHz.store(0.0, std::memory_order_relaxed);
    oBufferResampled.Reset();
    resampler_.Start(&oBufferFilter, &mtxFilterRing, &CvFilterUser,
                     &oBufferResampled, &mtxResampRing, &CvResampData,
                     &CvResampData, &StopAll);
    timingTracking_.ResetGardner(2.0, 1.0e-4, 1.0e-6, 64);
    timingTracking_.Start(&oBufferResampled, &mtxResampRing, &CvResampData, &StopAll);
    phaseTrackingDD_.Start(&timingTracking_, &StopAll);

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
    if(FilterThread.joinable())
        FilterThread.join();
    resampler_.StopJoin();
    timingTracking_.StopJoin();
    phaseTrackingDD_.StopJoin();
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

void Receiver::OperateFilter(void)
{
    auto throughput_window_start = std::chrono::steady_clock::now();
    uint64_t samples_in_window = 0;
    uint64_t samples_total = 0;
    const auto wall_start = std::chrono::steady_clock::now();
    auto symrate_display_start = std::chrono::steady_clock::now();
    auto next_symrate_start = std::chrono::steady_clock::now();
    bool symrate_collecting = false;
    
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
                    std::cout << "[SymbolRateEstimator] Detected: " << (est.SymbolRateHz / 1e6)
                              << " Msym/s"
                              << " [window " << SymRateBatchCounter << "/" << SymRateWindowBatches << "]"
                              << " t_sim=" << t_sim << " s"
                              << " t_rate=" << t_rate << " s"
                              << " (search kMax=" << est.SearchKMaxBins
                              << ", fMax=" << est.SearchMaxOffsetHz << " Hz)"
                              << " (FFT res " << est.FftResolutionHz
                              << " Hz, peak/median " << est.PeakToMedian << ")" << std::endl;
                }
                else
                {
                    const double t_rate = static_cast<double>(samples_total) / SamplingFrequency;
                    const double t_sim =
                        std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
                    std::cout << "[SymbolRateEstimator] Not detected"
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
            }
            else if (est.Detected)
            {
                if (timingTracking_.IsLocked())
                {
                    // Skip re-estimation while Gardner is locked to avoid delocking transients.
                    symrate_collecting = false;
                    SymRateBatchCounter = 0;
                    continue;
                }
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
                    std::cout << "[SymbolRateEstimator] Re-estimated: " << (refined / 1e6)
                              << " Msym/s"
                              << " [window " << SymRateBatchCounter << "/" << SymRateWindowBatches << "]"
                              << " t_sim=" << t_sim << " s"
                              << " t_rate=" << t_rate << " s"
                              << " (delta " << delta_ppm << " ppm, peak/median " << est.PeakToMedian << ")"
                              << std::endl;
                }
            }
            else
            {
                const double t_rate = static_cast<double>(samples_total) / SamplingFrequency;
                const double t_sim =
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
                std::cout << "[SymbolRateEstimator] Re-estimation not detected"
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
            double dt = std::chrono::duration<double>(now_tp - throughput_window_start).count();
            if (dt >= 1.0)
            {
                double msps = static_cast<double>(samples_in_window) / dt / 1e6;
                std::cout << "RX Filter throughput: " << msps << " Msps" << std::endl;

                throughput_window_start = now_tp;
                samples_in_window = 0;
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
            oViterbi[IndexViterbi].reset_decoder();
            float m[4] = {};
            oViterbi[IndexViterbi].Decode(SplitI[IndexViterbi]+PtrRd, SplitQ[IndexViterbi] + PtrRd, BatchSize1, ViterbiOutputs[IndexViterbi]+1, false, 1, 1, m[0]);
            oViterbi[IndexViterbi].reset_decoder();
            oViterbi[IndexViterbi].Decode(SplitI[IndexViterbi]+PtrRd, SplitQ[IndexViterbi] + PtrRd, BatchSize1, ViterbiOutputs[IndexViterbi]+1, false, 1, -1, m[1]);
            oViterbi[IndexViterbi].reset_decoder();
            oViterbi[IndexViterbi].Decode(SplitI[IndexViterbi]+PtrRd, SplitQ[IndexViterbi] + PtrRd, BatchSize1, ViterbiOutputs[IndexViterbi]+1, true, -1, 1, m[2]);
            oViterbi[IndexViterbi].reset_decoder();
            oViterbi[IndexViterbi].Decode(SplitI[IndexViterbi]+PtrRd, SplitQ[IndexViterbi] + PtrRd, BatchSize1, ViterbiOutputs[IndexViterbi]+1, true,   1, 1, m[3]);
            for (int k = 0; k < 4; ++k)
                VitSyncResults[IndexViterbi].Metrics[k] = m[k];
            VitSyncResults[IndexViterbi].Available = true;
            DemodulatorQ[IndexViterbi].AdvanceRead();
            lk.unlock();

        }
        else
        {
            oViterbi[IndexViterbi].Decode(SplitI[IndexViterbi]+PtrRd, SplitQ[IndexViterbi] + PtrRd, BatchSize1, ViterbiOutputs[IndexViterbi]+1, ViterbiParams.ExchangeIQ, ViterbiParams.SignI, ViterbiParams.SignQ, VitSyncResults[IndexViterbi].Metrics[0]);
            
            CvViterbis2VitManager[IndexViterbi].wait(lk, [&] {
                return StopAll || DecodedQ[IndexViterbi].AvailableWrite();
            });
            if (StopAll)
                break;
            int PtrWr = BatchSize1 * static_cast<int>(DecodedQ[IndexViterbi].GetPtrWr());
            DifferentialDecode(ViterbiOutputs[IndexViterbi],DiffDec[IndexViterbi]+PtrWr, BatchSize1, Last);
            DemodulatorQ[IndexViterbi].AdvanceRead();
            DecodedQ[IndexViterbi].AdvanceWrite();
            lk.unlock();
        }
        CvViterbis2VitManager[IndexViterbi].notify_one();

    }
}

void Receiver::OperateViterbiManager(void)
{
    NumBitsAll = 0;
    NumErrorsAll = 0;
    unsigned int PRBSSeed = 0;
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
                    bool PassThresh[3];
                    __m128 mThresh = _mm_set1_ps(ViterbiThreshold1);
                    bool Success = true;
                    for(int i = 0; i < 3; i++)
                    {
                        BestIndex[i] = 0;
                        float BestMetrics = VitSyncResults[i].Metrics[0];
                        for(int j = 1; j < 4; j++)
                        {
                            if(VitSyncResults[i].Metrics[j] < BestMetrics)
                            {
                                BestIndex[i] = j;
                                BestMetrics = VitSyncResults[i].Metrics[j];
                            }
                        }
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
                            for(int i = 0; i < 3; i++)
                                oViterbi[i].reset_decoder();
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
                        if(!PRBSSynchronized)
                        {
                            int PtrStart;
                            PRBSSynchronized = SyncPRBS(OutputAll+PtrRdOut, PtrStart, PRBSSeed);
                            if(PRBSSynchronized)
                            {
                                oPrbs.CreateOutputs(PRBSSeed, BatchSize3 - PtrStart, PrbsOut);
                                CountErrors(OutputAll+PtrRdOut+PtrStart, PrbsOut, BatchSize3-PtrStart);
                            }
                        }
                        else {
                            oPrbs.CreateOutputs(PRBSSeed, BatchSize3, PrbsOut);
                            CountErrors(OutputAll+PtrRdOut, PrbsOut, BatchSize3);
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
bool Receiver::SyncPRBS(unsigned char *In, int &PtrStart,unsigned int &Seed)
{
    bool Success = false;


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
            PtrStart = Ptr1 + 128;
        }
        else {
        
            Ptr+= 128;
        }
    }
    return Success;
}