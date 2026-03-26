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
#include "GardnerTiming.h"
#include <atomic>
using namespace std;
//#define DEBUG1
//#define DEBUG2
#define DEBUG_STATISTICS

// Uncomment to bypass Gardner (TakeEven 2x->1 sps), dump Viterbi input symbols to ../data/viterbi_input_bypass.bin
// (float32 interleaved I,Q per symbol), then compare offline with Gardner on-time dump via compare_viterbi_gardner_dump.m
//#define BYPASS_GARDNER_DUMP_VITERBI_INPUT

// Insert a fixed sampling-clock offset (SCO) at Gardner input (2 sps stream).
// This simulates a constant rhythm error without using the Resampler module.
#define DEBUG_GARDNER_INPUT_SCO
#define DEBUG_GARDNER_INPUT_SCO_PPM 1000.0

class Sampler;
struct DebugStatistics
{
    double FirstDelta;
    double MaxMetricsGrowth;
    double SumMetricsGrowth;
    double MeanMetricsGrowth;
    uint64_t NumBatches;
};


struct ViterbiSyncResults
{
    bool Available = false;
    float Metrics[4];
};
struct ViterbiParameters
{
    bool ExchangeIQ;
    float SignI, SignQ;
};
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
    GardnerTiming objGardnerTiming;
    Sampler *pSampler;
    BufferFloat oBufferFilter;
    Viterbi oViterbi[3] = {Viterbi(0),Viterbi(1),Viterbi(2)};
    bool StopAll;
    float *SplitI[3], *SplitQ[3];
    unsigned char *Outputs[3];
    float *OneSpsI, *OneSpsQ;
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
    condition_variable CvVitManager2Vit;
    condition_variable CvViterbis2VitManager[3];
    std::mutex mtxFilterRing;
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
    bool SyncPRBS(unsigned char *Input, int &PtrStart, unsigned int& Seed);
    void CountErrors(unsigned char *Input, unsigned char *Template, int Length);
    TxModes RxMode;
public:
#ifdef DEBUG_STATISTICS

    DebugStatistics CurrDebugStatistics;
#endif
    Receiver(/* args */);
    ~Receiver();
    void SetSampler(Sampler *p)
    {
        pSampler = p;
    }
    void StartThreads(double RollOff, TxModes RxModeIn,
                      const SymbolRateEstimatorConfig& symRateCfg = SymbolRateEstimatorConfig{});
    void StopThreads(void);
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
    void AdvanceOutput(void)
    {
        {
            std::lock_guard<std::mutex> lk(mtxQueues);
            OutputQ.AdvanceRead();
        }
        CvOut2Rx.notify_one();
    }
    void Get_Condition_Variables(condition_variable * pcv_rx_out, condition_variable * pcv_out_rx)
    {
        pcv_out_rx = &CvOut2Rx;
        pcv_rx_out = &CvRx2Out;
    }
    uint64_t NumErrorsAll, NumBitsAll;
};
