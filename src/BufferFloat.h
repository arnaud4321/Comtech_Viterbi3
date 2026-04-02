/**
 * @file BufferFloat.h
 * @brief Lock-free ring buffer contract for interleaved float I/Q between RX pipeline stages.
 */
#pragma once
#include "cstdlib"
#include <algorithm>    // std::copy, std::max
#include <cassert>
#include "definitions.h"
#include <iostream>
using namespace std;

/**
 * @brief Circular I/Q buffer (float). External mutex serializes pointer updates.
 *
 * Producer must wait for @ref AlmostFull() == false before @ref AdvancePtrWr after writing.
 */
class BufferFloat
{
private:
    /* data */
    
    int BufferSize;
    int ExtraBufferSize;
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

