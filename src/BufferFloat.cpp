#include "BufferFloat.h"
#include "ConsoleAlert.h"
#include <cstdio>
BufferFloat::BufferFloat(int BufferSizeIn, int ExtraBufferSizeIn)
    : BufferSize(BufferSizeIn), ExtraBufferSize(ExtraBufferSizeIn)
{
    int AllocLen = (BufferSize + ExtraBufferSize);
    datai = (float *)std::aligned_alloc(32,  2*AllocLen* sizeof(datai));
    dataq = datai + AllocLen;
    Mask = BufferSize - 1;
    // "Almost full" threshold: trigger backpressure when ~50% of the ring is occupied
    // (previously: ~87.5%), to absorb a very slow Viterbi (e.g. DEBUG_GARDNER_OUTPUTS).
    GuardSize = std::max(SPB * 16, BufferSize / 2);
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

int BufferFloat::GetSizeInBuffer(void) const
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
                if (ExtraAtEnd > ExtraBufferSize)
                {
                    CONSOLE_ALERT_STMT(std::fprintf(stderr,
                                                    "%sFATAL: BufferFloat wrapped read exceeds ExtraBufferSize "
                                                    "(ExtraAtEnd=%d, ExtraBufferSize=%d, Size=%d, PtrRd=%d, "
                                                    "BufferSize=%d)%s\n",
                                                    ConsoleAlert::kRedOpen, ExtraAtEnd, ExtraBufferSize, Size, PtrRd,
                                                    BufferSize, ConsoleAlert::kReset););
                    std::abort();
                }
                assert(ExtraAtEnd <= ExtraBufferSize && "BufferFloat: ExtraBufferSize too small for wrapped read");
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
            const int overflow = (NewPtrWr - BufferSize);
            if (overflow > ExtraBufferSize)
            {
                CONSOLE_ALERT_STMT(std::fprintf(stderr,
                                                "%sFATAL: BufferFloat wrapped write exceeds ExtraBufferSize "
                                                "(overflow=%d, ExtraBufferSize=%d, Advance=%d, PtrWr=%d, "
                                                "BufferSize=%d)%s\n",
                                                ConsoleAlert::kRedOpen, overflow, ExtraBufferSize, Advance, PtrWr,
                                                BufferSize, ConsoleAlert::kReset););
                std::abort();
            }
            assert(overflow <= ExtraBufferSize && "BufferFloat: ExtraBufferSize too small for wrapped write");
            std::copy(datai + BufferSize,datai + NewPtrWr,datai);//copy the end to the beggining
            std::copy(dataq + BufferSize,dataq + NewPtrWr,dataq);//copy the end to the beggining
        }
        NewPtrWr -= BufferSize;
    }    
    PtrWr = NewPtrWr;
}
    

void BufferFloat::AdvancePtrRd(int Advance)
{
    int NewPtrRd = PtrRd + Advance;
    while (NewPtrRd >= BufferSize)
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

bool BufferFloat::CanCommitWrite(int advance) const
{
    if (advance <= 0)
        return true;
    const int delta = GetSizeInBuffer();
    return delta + advance <= BufferSize - 1;
}
