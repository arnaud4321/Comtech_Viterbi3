#pragma once
#include <chrono>
#include <thread>
#include <condition_variable>
#include "definitions.h"
class AWGNChannel;
#include "BufferShort.h"
using namespace std;
//#define DEBUG_SAMPLER

class Sampler
{
private:
    /* data */
    
    #ifdef DEBUG_SAMPLER
    short *OutAllI = 0;
    unsigned int PtrOutAll = 0;
    #define STOP_SAMPLE_S 6000000
    #endif
    
    bool Debug;
    bool StopAll = false;
    double TimeBatch;
    thread SamplingThread;
    AWGNChannel *pChannel;
    BufferShort oBuffer;
    condition_variable  CvSamplerUser;
    void OperateSampler(void);
public:
    Sampler(bool DebugIn);
    ~Sampler();
    void StartThread(void);
    void StopThread(void);
    void SetChannel(AWGNChannel *p)
    {
        pChannel = p;
    }
    short * GetOutput(int Size)
    {
        return oBuffer.GetReadBuffer(Size);
    }

    void AdvanceOut(int Size)
    {
        oBuffer.AdvancePtrRd(Size);
    }
    condition_variable * GetCvOut(void)
    {
        return &CvSamplerUser;
    }
};

