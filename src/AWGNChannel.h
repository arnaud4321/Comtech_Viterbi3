#pragma once
#include "random_generator_new.h"
class Transmitter;
#include <thread>
#include <condition_variable>
#include <immintrin.h>
#include "SimpleQueue.h"
#include "BufferShort.h"
#include "definitions.h"
using namespace std;

//#define DEBUG_AWGN

class AWGNChannel
{
private:
    /* data */
    
    #ifdef DEBUG_AWGN
    short *OutAllI = 0;
    unsigned int PtrOutAll = 0;
    #endif
    
    const double Backoff = 12;
    const int NumBits = 16;
    int NoiseBatchSize;
    random_generator_new oNoiseGen;
    double Stdn; 
    Transmitter *pTx;
    static constexpr int LengthQueue = 4;
    float *Noise[LengthQueue];
    thread NoiseThread, OutputThread;
    void GenerateNoise(void);
    void GenerateOutput(void);
    SimpleQueue NoiseQ;
    bool StopAll = false;
    condition_variable CvOutNoise, CvNoiseOut, CvOutUser, CvUserOut;
    __m256i floats_to_shorts_sat_perm_avx2(__m256 x, __m256 y);
    __m256 mkSig;

public:
    AWGNChannel(unsigned int Seed);
    ~AWGNChannel();
    void StartThreads(float EsN0db, double TimeDrift, double FrequencyShift, bool RandomFrequency);
    void StopThreads(void);
    void SetTransmitter(Transmitter *p)
    {
        pTx = p;
    }
    short * GetOutput(int Size);
    void AdvanceOut(int Size);
    void GetCvs2User(condition_variable * &CvIn, condition_variable * &CvOut )
    {
        CvOut = &CvOutUser;
        CvIn = &CvUserOut;
    }
    BufferShort OutputBuffer;

};

