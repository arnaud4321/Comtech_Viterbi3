#pragma once
#include "definitions.h"
#include <thread>
#include <condition_variable>
#include <immintrin.h>
#include "RxFilter.h"
#include "BufferFloat.h"
using namespace std;
//#define DEBUG1
class Sampler;
class Receiver
{
private:
#ifdef DEBUG1
    float *OutAllI = 0;
    float *OutAllQ = 0;
    unsigned int PtrOutAll = 0;
    #define NUM_SAMPLES_R1 19000000
#endif
    RxFilter objRxFilter;
    Sampler *pSampler;
    BufferFloat oBufferFilter;
    bool StopAll;
    void OperateFilter(void);
    void OperateViterbi(void);
    thread FilterThread, ViterbiThread;
    float *FilterInI,  *FilterInQ;
    void shorts_to_floats_avx2(__m256i v16, __m256* out0, __m256* out1);
    condition_variable CvFilterUser;
public:
    Receiver(/* args */);
    ~Receiver();
    void SetSampler(Sampler *p)
    {
        pSampler = p;
    }
    void StartThreads(double RollOff);
    void StopThreads(void);


};
