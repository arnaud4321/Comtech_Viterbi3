#include "SimpleQueue.h"
SimpleQueue::SimpleQueue(int LengthIn)
{
    Length = LengthIn;
    Mask = Length - 1;
    Available = new bool [Length];
    for(int i = 0; i < Length; i++)
        Available[i] = false;
}

SimpleQueue::~SimpleQueue(void)
{
    delete [] Available;
}
bool SimpleQueue::AvailableRead(void)
{
    if(Available[PtrRd] == true)
    {
        return true;
    }
    else
        return false;
}
bool SimpleQueue::AvailableWrite(void)
{
    if(Available[PtrWr] == false)
    {
        return true;
    }
    else
        return false;

}
void SimpleQueue::AdvanceWrite(void)
{
    Available[PtrWr] = true;
    int NewPtr = PtrWr;
    NewPtr = (NewPtr +1) & Mask;
    PtrWr = NewPtr;
    
}
void SimpleQueue::AdvanceRead(void)
{
    Available[PtrRd] = false;
    int NewPtr = (PtrRd);
    NewPtr = (NewPtr +1) & Mask;
    PtrRd = NewPtr;

}

void SimpleQueue::Reset(void)
{
    for(int i = 0; i < Length; i++)
        Available[i] = false;
}