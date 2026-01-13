#include "ConvEncoder.h"

ConvEncoder::ConvEncoder()
{
	//Derived Constants

	intNParam = 2;
	intCLength = 7;


	intNoStates = 1 << (intCLength -1);
	intMask = intNoStates - 1;

	OutputTable = new unsigned int [intNoStates * intNParam * 2];
	NextStates = new unsigned int [intNoStates << 1]; 
	
	calc_encoding_table();
	InitBit2Byte(1);
	Outputs8 = ippsMalloc_8u(intNoStates * 256 * intNParam * 8);//For handling Byte Inputs. NumStates x 2 x 8 x 256
	NextStates8 = ippsMalloc_8u(256);

	calc_encoding_table_8();
	Reset();
#ifdef DEBUG_CONV_ENC

	OutAll = ippsMalloc_8u(5000000);
	OutAll0 = ippsMalloc_8u(5000000);

	
#endif
}

ConvEncoder::~ConvEncoder(void)
{

	delete [] OutputTable;
	
	ippsFree(NextStates8);
	ippsFree(Outputs8);
	
#ifdef DEBUG_CONV_ENC

	ippsFree(OutAll);

#endif
}

int ConvEncoder::is_odd_ones(int input)
{
	return (_mm_popcnt_u32(input) & 0x1);
	
}

void ConvEncoder::calc_encoding_table()
{
	for(unsigned int ii = 0; ii < intNoStates; ii ++)
	{
		for(unsigned int jj =0; jj < intNParam; jj++)
		{
			//Zero Input
			OutputTable[2*ii*intNParam + jj] = is_odd_ones(ii&ivecGenPolys[jj]);
			//One Input
			OutputTable[2*ii*intNParam + intNParam + jj] = is_odd_ones((intNoStates+ii)&ivecGenPolys[jj]);

		}
		
		NextStates[ii<<1] = ii >> 1;//Next State with Zero Input
		NextStates[(ii<<1)+1] = (ii >> 1) + (intNoStates>>1);//One Input
	}
}


void ConvEncoder::InitBit2Byte(int Mode)
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

void ConvEncoder::calc_encoding_table_8()
{
	alignas(32) Ipp8u Bits[32];

	for(unsigned int ii = 0; ii < intNoStates; ii ++)
	{
		
		
		unsigned int BasePtr0 = (ii << (8+4)); //256 Bytes * 16 Bits
		for(unsigned int Byte = 0; Byte < 256; Byte++)
		{
			//load the Byte

			std::memcpy(Bits, Byte2BitTable + (Byte << 3), 8);

			unsigned int BasePtr = BasePtr0 + (Byte << (4)); 

			unsigned int NextState8 = ii;

			for(unsigned int NumBit = 0; NumBit < 8; NumBit++)
			{
				unsigned char CurrentBit = Bits[NumBit];
				for(unsigned int jj = 0; jj < intNParam; jj++)
				{
					unsigned char CurrentOutputBit = (unsigned char)OutputTable[((NextState8 << 1) + (((unsigned int) CurrentBit))) * intNParam + jj];
					Outputs8[BasePtr+ (NumBit<<1) +jj] = CurrentOutputBit;
				}
				NextState8 = NextStates[(NextState8<<1)+Bits[NumBit]];
			}

			NextStates8[Byte] = NextState8;
		}
	}
}

void ConvEncoder::Encode(unsigned char* Input, unsigned char* Out, unsigned int InputLength)
{
	
	int PtrOut = 0;
	for(unsigned int ii = 0; ii < InputLength; ii ++)
	{	
		//Read the Encoded Values from the Table
		unsigned int Ptr = (InitialState<<12)+(Input[ii]<<4);
		_mm_storeu_si128((__m128i*)(Out + PtrOut), _mm_loadu_si128((__m128i*)(Outputs8 + Ptr)));
		InitialState = NextStates8[Input[ii]];
		PtrOut += 16;
	}
	

}



void ConvEncoder::Bit2Bytes2R(Ipp8u* In, short *Out)
{
	__m128i mIn = _mm_loadu_si128((__m128i*) In);
	mIn = _mm_shuffle_epi8(mIn, _mm_set_epi8(8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7)); //Reverse bits
	mIn = _mm_cmpgt_epi8(mIn, _mm_setzero_si128());
	*Out = (((short)_mm_movemask_epi8(mIn)) );

}