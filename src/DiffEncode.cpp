#include "DiffEncode.h"
DiffEncode::DiffEncode(void)
{
    InitBit2Byte(1);
    InitEncodingTable();
}

DiffEncode::~DiffEncode()
{

}

void DiffEncode::Reset()
{
    LastState = 0;
}


void DiffEncode::InitBit2Byte(int Mode)
{
	for (unsigned int jj = 0; jj < 256; jj++)
	{
		unsigned int Tmp  = jj;
		if (Mode == 0)
		{
			for (int kk = 7; kk >= 0; kk--)
			{
				Byte2BitTable[(jj << 3) + kk] = Tmp & 1;
				Tmp = (Tmp >> 1);
			}
		}
		else
		{
			for (int kk = 0; kk < 8; kk++)
			{
				Byte2BitTable[(jj << 3) + kk] = Tmp & 1;
				Tmp = (Tmp >> 1);
			}
		}
	}
}

void DiffEncode::InitEncodingTable(void)
{
    for(unsigned int State = 0; State < 2; State++)
    {
        

        for(int Byte = 0; Byte < 256; Byte++)
        {
            unsigned int StartState = State;
            unsigned int Out = 0;
            for(int NumBit = 0; NumBit < 8; NumBit++)
            {
                unsigned char CurrBit = Byte2BitTable[(Byte<<3)+NumBit];
                unsigned char New = CurrBit ^ (unsigned char ) StartState;
                Out = Out | (New << NumBit);
                StartState = (unsigned int ) New;
            }
            EncodingTable[256*State+Byte] = Out;
        }
    }
}

void DiffEncode::CreateOutputs(unsigned char *In, unsigned char *Out, unsigned int Length)
{
    for(unsigned int i = 0; i < Length; i++)
    {
        Out[i]  =  EncodingTable[256*LastState+In[i]];
        LastState = (Out[i] >> 7) & 1;
    }
}