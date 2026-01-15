#include "BufferFloat.h"
BufferFloat::BufferFloat(int BufferSizeIn, int ExtraBufferSize):BufferSize(BufferSizeIn)
{
    int AllocLen = (BufferSize + ExtraBufferSize);
    datai = (float *)std::aligned_alloc(32,  2*AllocLen* sizeof(datai));
    dataq = datai + AllocLen;
    Mask = BufferSize - 1;    
    GuardSize = BufferSize - (BufferSize>>3);
}

BufferFloat::~BufferFloat()
{
    std::free(datai);
}



void BufferFloat::GetWriteBuffer(float * &OutI, float * &OutQ, int Size) //if at End you need to copy ExtraAtEnd at the Start
{
   
    OutI = datai + PtrWr;
    OutQ = dataq + PtrWr;
}

int BufferFloat::GetSizeInBuffer(void)
{
    int Delta = PtrWr - PtrRd;
    if(Delta < 0)
    {
        Delta += BufferSize;
    }
    return Delta;
}

void BufferFloat::GetReadBuffer(float * &OutI, float * &OutQ, int Size)
{
    int ExtraAtEnd;
    OutI = 0;
    OutQ = 0;
    int Delta = GetSizeInBuffer();

    if(Delta >= Size)
    {
        OutI = datai + PtrRd;
        OutQ = dataq + PtrRd;
        
        {
            int PtrEnd = PtrRd + Size;
            if(PtrEnd > BufferSize)
            {
                ExtraAtEnd = (PtrEnd - BufferSize) ;
                std::copy(datai,datai+ExtraAtEnd,datai+BufferSize);
                std::copy(dataq,dataq+ExtraAtEnd,dataq+BufferSize);
            }
        }
    }
}

void BufferFloat::AdvancePtrWr(int Advance)
{
    int NewPtrWr = PtrWr;
    NewPtrWr += Advance;
    if((NewPtrWr) >= BufferSize)
    {
        if(NewPtrWr > BufferSize)
        {
            std::copy(datai + BufferSize,datai + NewPtrWr,datai);//copy the end to the beggining
            std::copy(dataq + BufferSize,dataq + NewPtrWr,dataq);//copy the end to the beggining
        }
        NewPtrWr -= BufferSize;
    }    
    PtrWr = NewPtrWr;
}
    

void BufferFloat::AdvancePtrRd(int Advance)
{
    int NewPtrRd = PtrRd;
    NewPtrRd += Advance;
    if((NewPtrRd) >= BufferSize)
    {
        NewPtrRd -= BufferSize;
    }
    PtrRd = NewPtrRd;
}

bool BufferFloat::AlmostFull()
{
    if(GetSizeInBuffer() >= GuardSize)
    {
        return true;
    }
    else 
    {
        return false;
    }
}
