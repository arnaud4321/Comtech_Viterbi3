#pragma once

class SimpleQueue
{
    private:

        bool *Available;
        unsigned int  PtrRd = 0;
        unsigned int PtrWr = 0;
        unsigned int Length;
        unsigned int Mask;    
    public:
    SimpleQueue(int LengthIn);
    ~SimpleQueue(void);
    bool AvailableRead(void);
    bool AvailableWrite(void);
    void AdvanceWrite(void);    
    void AdvanceRead(void);
    void Reset(void);
    unsigned int GetPtrRd()
    {
        return PtrRd;
    }
    unsigned int GetPtrWr()
    {
        return PtrWr;
    }

};