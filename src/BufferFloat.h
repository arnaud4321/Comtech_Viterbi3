#pragma once
#include "cstdlib"
#include <algorithm>    // std::copy, std::max
#include "definitions.h"
#include <iostream>
using namespace std;
// I/Q ring (floats). Multithread contract: all reads/writes of PtrRd/PtrWr must be
// serialized by the caller (e.g. mtxFilterRing). GetWriteBuffer does not check capacity:
// the producer must wait for AlmostFull() == false before CreateOutputs + AdvancePtrWr(advance).
class BufferFloat
{
private:
    /* data */
    
    int BufferSize;
    int GuardSize;
    float *datai = 0, *dataq = 0;
    int PtrWr = 0, PtrRd = 0, Mask;
    bool First = true;
public:
    BufferFloat(int BufferSize, int ExtraBufferSize);
    ~BufferFloat();
    void GetWriteBuffer(float * &OutI, float * &OutQ, int Size);
    void GetReadBuffer(float * &OutI, float * &OutQ, int Size);
    void GetBufferStart(float * &OutI, float * &OutQ)
    {
        OutI = datai;
        OutQ = dataq;
    }
    void GetBufferEnd(float * &OutI, float * &OutQ)
    {
        OutI =  datai+BufferSize;
        OutQ =  dataq+BufferSize;
        
    }
    void AdvancePtrWr(int Advance);
    void AdvancePtrRd(int Advance);
    
    void Reset(void)
    {
        PtrRd = 0;
        PtrWr = 0;
    }
    int GetWrPtr(void)
    {
        return PtrWr;
    }
    int GetRdPtr(void)
    {
        return PtrRd;
    }
    int GetBufferSize(void)
    {
        return BufferSize;
    }
    bool AlmostFull(void);
    int GetSizeInBuffer(void) const;
    // After the physical write [PtrWr, PtrWr+advance), verify that an AdvancePtrWr(advance)
    // would not push the fill level beyond BufferSize-1 (ring without wasted slot).
    bool CanCommitWrite(int advance) const;
};

