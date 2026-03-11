#pragma  once
#include "ipp.h"
#include <immintrin.h>
#include <cstring>
#include <cstdint>


//#define DEBUG_CONV_ENC

class ConvEncoder
{
#ifdef DEBUG_THREADING

	unsigned __int64 NumSamplesIn, NumSamplesOut;
	Ipp64u TimeWaiting, TimeAll, TimeIdle;

#endif
#ifdef DEBUG_CONV_ENC

	Ipp8u* OutAll, *OutAll0;
	unsigned int PtrOutAll = 0;
	unsigned int PtrOutAll0 = 0;

#endif
protected:
	unsigned int InitialState = 0; //Start from 0
	unsigned int intCLength; //Constraint Length
	unsigned int intNParam; //the n parameter of 1/n
	unsigned int intNoStates, intMask;
	unsigned int *OutputTable, *NextStates;
	int ivecGenPolys[2] = { 0133,0171 };
	Ipp8u Byte2BitTable[256*8+32];
	//Private Methods
	virtual void calc_encoding_table(void);
	int is_odd_ones(int input);
	virtual void calc_encoding_table_8(void);
	void InitBit2Byte(int Mode);
	Ipp8u *NextStates8; //Byte Input
	Ipp8u *Outputs8;
	
	void Bit2Bytes2R(Ipp8u* In, short *Out);
public:
	ConvEncoder();
	~ConvEncoder(void);
	void Encode(unsigned char *Input, unsigned char *Out, unsigned int InputLength);//Byte input, LSB first
	void Reset(void)
	{
		InitialState = 0;
	#ifdef DEBUG_THREADING
		NumSamplesIn = NumSamplesOut = 0;
		TimeWaiting =  TimeAll = TimeIdle = 0;
#endif


	}
	
//	virtual void CreateOutput();
};

