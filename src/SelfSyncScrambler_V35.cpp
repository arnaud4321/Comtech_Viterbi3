#include "SelfSyncScrambler_V35.h"


SelfSyncScrambler_V35::SelfSyncScrambler_V35() 
{
	DiffEnc = false;
	Reset();

	// Init Masks
	MASK = (0x1 << 2) | (0x1 << 19);	// Shift register taps
	MASK2 = 0xFFFFF;					// 2^20-1



}

SelfSyncScrambler_V35::~SelfSyncScrambler_V35()
{
}



void SelfSyncScrambler_V35::Scramble(unsigned char* In, unsigned char* Out, unsigned int Length)
{
	
	
	
	//Ipp16u OutTest2[10000];

	//NumBytesPrev = 20;

	for (unsigned int i = 0; i < Length; i ++)
	{
		
			unsigned int flag = (cnt == 31);


			Out[i] = (((shiftReg>>2) & 1) ^ (shiftReg >> 19) ^ flag ^ In[i])^1;

			// Set cnt
			// cnt = bitand(bitand(cnt + (Delayline(1) == Delayline(9)), 31), 31 * (Delayline(1) == Delayline(9)));
			bool Val19 = ((shiftReg & 0x1) == ((shiftReg >> 8) & 0x1));
			cnt = (++cnt & 0x1F) * Val19;

			// Shift left by one bit
			shiftReg <<= 1;

			// Add new input bit
			shiftReg |= Out[i];

			// Shave off additional bits not found in shift register
			shiftReg &= MASK2;
			
		}

}

void SelfSyncScrambler_V35::Descramble(unsigned char* In, unsigned char* Out, unsigned int Length)
{



	//Ipp16u OutTest2[10000];

	//NumBytesPrev = 20;

	for (unsigned int i = 0; i < Length; i++)
	{

		unsigned int flag = (cnt == 31);


		Out[i] = (((shiftReg >> 2) & 1) ^ (shiftReg >> 19) ^ flag ^ In[i]) ^ 1;

		// Set cnt
		// cnt = bitand(bitand(cnt + (Delayline(1) == Delayline(9)), 31), 31 * (Delayline(1) == Delayline(9)));
		bool Val19 = ((shiftReg & 0x1) == ((shiftReg >> 8) & 0x1));
		cnt = (++cnt & 0x1F) * Val19;

		// Shift left by one bit
		shiftReg <<= 1;

		// Add new input bit
		shiftReg |= In[i];

		// Shave off additional bits not found in shift register
		shiftReg &= MASK2;

	}

}
