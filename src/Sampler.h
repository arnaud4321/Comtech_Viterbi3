#pragma once
#include <chrono>
#include <thread>
#include <condition_variable>
#include <mutex>
#include "definitions.h"
#include "random_generator_new.h"
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
    std::mutex mtxSamplerBuffer_;
    condition_variable CvSamplerUser;
    void OperateSampler(void);
public:
    uint64_t NumBatches;
    Sampler(bool DebugIn);
    ~Sampler();
    void StartThread(void);
    void StopThread(void);
    void SetChannel(AWGNChannel *p)
    {
        pChannel = p;
    }
    /// Lecture + AdvancePtrRd sous mutex (thread filtre vs thread sampler).
    bool ReadFilterBatch(short* dst, int nShorts, bool& stopAll);
    void NotifyFilterWaiters() { CvSamplerUser.notify_all(); }
    condition_variable* GetCvOut(void)
    {
        return &CvSamplerUser;
    }
};

