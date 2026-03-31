#include "BufferShort.h"
BufferShort::BufferShort(int BufferSizeIn, int ExtraBufferSizeIn)
    : BufferSize(BufferSizeIn), ExtraBufferSize(ExtraBufferSizeIn)
{
    data = (short*) _mm_malloc((BufferSize + ExtraBufferSize) * sizeof(short), 32);
    Mask = BufferSize - 1;
    GuardSize = BufferSize - (BufferSize>>3);
}

BufferShort::~BufferShort()
{
   _mm_free(data);
}



short *BufferShort::GetWriteBuffer(int Size) //if at End you need to copy ExtraAtEnd at the Start
{
   
    short *RetVal = data + PtrWr;
    return RetVal;
}

int BufferShort::GetSizeInBuffer(void)
{
    int Delta = PtrWr - PtrRd;
    if(Delta < 0)
    {
        Delta += BufferSize;
    }
    return Delta;
}

short * BufferShort::GetReadBuffer(int Size)
{
    int ExtraAtEnd;
    short *RetVal = 0;
    int Delta = GetSizeInBuffer();
    if(Delta >= Size)
    {
        RetVal = data + PtrRd;
        int PtrEnd = PtrRd +Size;
        if(PtrEnd > BufferSize) 
        {
            ExtraAtEnd = PtrEnd - BufferSize;
            if (ExtraAtEnd > ExtraBufferSize)
            {
                std::fprintf(stderr,
                             "FATAL: BufferShort wrapped read exceeds ExtraBufferSize "
                             "(ExtraAtEnd=%d, ExtraBufferSize=%d, Size=%d, PtrRd=%d, BufferSize=%d)\n",
                             ExtraAtEnd, ExtraBufferSize, Size, PtrRd, BufferSize);
                std::abort();
            }
            assert(ExtraAtEnd <= ExtraBufferSize && "BufferShort: ExtraBufferSize too small for wrapped read");
            std::copy(data,data+ExtraAtEnd,data+BufferSize);
        } 
    }

    return RetVal;
}
void BufferShort::AdvancePtrWr(int Advance)
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
                std::fprintf(stderr,
                             "FATAL: BufferShort wrapped write exceeds ExtraBufferSize "
                             "(overflow=%d, ExtraBufferSize=%d, Advance=%d, PtrWr=%d, BufferSize=%d)\n",
                             overflow, ExtraBufferSize, Advance, PtrWr, BufferSize);
                std::abort();
            }
            assert(overflow <= ExtraBufferSize && "BufferShort: ExtraBufferSize too small for wrapped write");
            std::copy(data + BufferSize, data + NewPtrWr, data);//copy the end to the beggining
        }
        NewPtrWr -= BufferSize;
    }    
    PtrWr = NewPtrWr;
}
    

void BufferShort::AdvancePtrRd(int Advance)
{
    int NewPtrRd = PtrRd + Advance;
    while (NewPtrRd >= BufferSize)
        NewPtrRd -= BufferSize;
    PtrRd = NewPtrRd;
}

bool BufferShort::AlmostFull()
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
