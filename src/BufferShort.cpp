#include "BufferShort.h"
BufferShort::BufferShort(int BufferSize, int ExtraBufferSize)
{
    data = (short *)std::aligned_alloc(32, (BufferSize + ExtraBufferSize) * sizeof(short));
    Mask = BufferSize - 1;
}

BufferShort::~BufferShort()
{
    std::free(data);
}

short *BufferShort::GetWriteBuffer(int Size) //if at End you need to copy ExtraAtEnd at the Start
{
   
    short *RetVal = data + PtrWr;
    return RetVal;
}
short * BufferShort::GetReadBuffer(int Size)
{
    int ExtraAtEnd;
    short *RetVal = 0;
    int Delta = PtrWr - PtrRd;
    if(Delta >= 0)
    {
        ExtraAtEnd = 0;
    }
    else
    {
        Delta += BufferSize;
        ExtraAtEnd = 1;
    }
    if(Delta >= Size)
    {
        RetVal = data + PtrRd;
        if(ExtraAtEnd > 0)
        {
            ExtraAtEnd = (PtrRd+Size) & Mask;
            std::copy(data,data+ExtraAtEnd,data+BufferSize);
        }
    }

    return RetVal;
}
void BufferShort::AdvancePtrWr(int Advance)
{
    int NewPtrWr = PtrWr + Advance;
    if((NewPtrWr) >= BufferSize)
    {
        std::copy(data + BufferSize,data + NewPtrWr,data);//copy the end to the beggining
        NewPtrWr -= BufferSize;
    }    
    PtrWr = NewPtrWr;
}
    

void BufferShort::AdvancePtrRd(int Advance)
{
    int NewPtrRd = PtrRd + Advance;
    if((NewPtrRd) >= BufferSize)
    {
        NewPtrRd -= BufferSize;
    }
    PtrRd = NewPtrRd;
}

