#pragma once
#include "ipp.h"

class DiffEncode
{
private:
    /* data */
    Ipp8u Byte2BitTable[256*8+32];
    Ipp8u EncodingTable[512];
    
    void InitBit2Byte(int Mode);
    void InitEncodingTable(void);
    unsigned int LastState = 0;

public:
    void Reset(void);
    DiffEncode(/* args */);
    ~DiffEncode();
    void CreateOutputs(unsigned char *In, unsigned char *Out, unsigned int Length);
};

