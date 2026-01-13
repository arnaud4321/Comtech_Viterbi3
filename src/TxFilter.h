#pragma once
#include "ipp.h"
#include <math.h>
//#include "../Common//defs.h"

class TxFilter
{
	
	Ipp32f PI;
	Ipp32f *SqrtrcFilter = 0;
	unsigned int FilterLength;
	Ipp32f RollOff;
	Ipp32u TxNSS;
	IppsFIRSpec_32f *pSpec=0;
	Ipp8u* pBuf=0;
	Ipp32f *pDlySrcI = 0, *pDlySrcQ = 0, *pDlyDstI = 0, *pDlyDstQ = 0;
	unsigned int DataLength;
	
	
	Ipp32u BufferMask;	

	void GenerateSqrtrcFilter(Ipp32f Alpha, unsigned int HalfDelay, Ipp32u SamplingFrequency);
public:
	TxFilter();
	~TxFilter(void);
	void CreateObjects(Ipp32f RollofIn);
	void CreateOutputs(Ipp32f *InputI, Ipp32f *InputQ, Ipp32f *OutI, Ipp32f *OutQ,  unsigned int InputLength);
	void DeleteObjects(void);
	double TxPower;

};

