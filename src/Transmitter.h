#pragma once
#include <immintrin.h>
#include <thread>
#include <condition_variable>
#include <string>
using namespace std;
#include "definitions.h"
#include "Prbs23.h"
#include "ConvEncoder.h"
#include "SelfSyncScrambler_V35.h"
#include "DiffEncode.h"
#include "SimpleQueue.h"
#include "TxFilter.h"

//#define DEBUG_TX


class Transmitter
{
    private:

#ifdef DEBUG_TX
    float *OutAllI = 0;
    unsigned int PtrOutAll = 0;
#endif

    float Map[2] = {-1,1};
    unsigned char *Data = 0;
    unsigned char *ScrambledOut = 0;
    unsigned char *SplitOut[3];
    unsigned char *DiffOut[3];
    static constexpr int LengthQueue = 4;

    unsigned char *ConvOut[LengthQueue*3];
    float *FilterInI = 0, *FilterInQ = 0;
    float *FilterOutI = 0, *FilterOutQ = 0;
    float *TxOut[LengthQueue];
    unsigned int LengthData = 0;

    Prbs23 oPrbs;
    SelfSyncScrambler_V35 oScrambler;
    ConvEncoder oConvEncoder[3];
    DiffEncode oDiffEncode[3];
    TxFilter oTxFilter;
    thread DataGenThread, FilterThread;
    bool StopAll = false;
    SimpleQueue DataQ, FilterQ;
    condition_variable CvFilterData, CvDataFilter, CvOutFilter, CvFilterOut;

    void GenerateData(void);
    void FilterData(void);
    void Split3(void);
    void interleave_iq_f32_to_f32_avx2(const float* xi, const float* xq, float* y, int n) ;
    void interleave_iq_f32_to_i16_avx2(const float* xi, const float* xq, short* y, int n) ;
public:
    Transmitter(/* args */);
    ~Transmitter();
    void StartThreads(TxModes TxMode, float RollOff, string FileName = "");
    void StopThreads(void);
	float *GetOutput(void);
    void AdvanceOut(void)
    {
        FilterQ.AdvanceRead();
    }
    void GetCVOut(condition_variable* &CVIn, condition_variable* &CVOut )
    {
        CVIn = &CvOutFilter;
        CVOut = &CvFilterOut;
    }
    double GetTxPower()
    {
        return oTxFilter.TxPower;
    }
};