#pragma once
#include <ipp.h>
#include <immintrin.h>
using namespace std;

class SelfSyncScrambler_V35 
{
	bool DiffEnc;
	
	Ipp32u MASK, MASK2;
	Ipp32u shiftReg, cnt;


public:
	SelfSyncScrambler_V35();
	~SelfSyncScrambler_V35();
	void Scramble(unsigned char *In, unsigned char *Out, unsigned int Length);
	void Descramble(unsigned char* In, unsigned char* Out, unsigned int Length);

	void Reset()
	{ 
		shiftReg = cnt = 0; 
	};
};

