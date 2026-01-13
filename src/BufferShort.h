#pragma once
#include "cstdlib"
#include <algorithm>    // std::copy
#include "definitions.h"
#include <iostream>
using namespace std;
//Everything is managed in shorts
class BufferShort
{
private:
    /* data */
    
    int BufferSize;
    short *data = 0;
    int PtrWr = 0, PtrRd = 0, Mask;
    bool First = true;
public:
    BufferShort(int BufferSize, int ExtraBufferSize);
    ~BufferShort();
    short * GetWriteBuffer(int Size);
    short * GetReadBuffer(int Size);
    short *GetBufferStart(void)
    {
        return data;
    }
    short *GetBufferEnd(void)
    {
        return data+BufferSize;
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
    int GetBufferSize(void)
    {
        return BufferSize;
    }
};

