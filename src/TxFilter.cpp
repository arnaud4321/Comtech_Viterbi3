/**
 * @file TxFilter.cpp
 * @brief Transmit pulse shaping: sqrt raised cosine via IPP FIR, @c TxNSS samples per symbol, configurable rolloff.
 */

#include "TxFilter.h"


TxFilter::TxFilter()
{
	
	PI = IPP_PI;
	TxNSS = 2;
	
}


TxFilter::~TxFilter(void)
{
	
}


void TxFilter::DeleteObjects(void)
{
	ippsFree(SqrtrcFilter);
	ippsFree(pSpec);
	ippsFree(pBuf);
	ippsFree(pDlySrcI);
	ippsFree(pDlySrcQ);
	ippsFree(pDlyDstI);
	ippsFree(pDlyDstQ);
	SqrtrcFilter = 0;
	pSpec = 0;
	pBuf = 0;
	pDlySrcI = 0;
	pDlySrcQ = 0;
	pDlyDstI = 0;
	pDlyDstQ = 0;

}
void TxFilter::CreateObjects(double RollOffIn)
{
	//Calculate Length
	RollOff = RollOffIn;
	DeleteObjects();
	unsigned int Delay;
	Delay = ceil(4.0/RollOff);
	FilterLength = 2*Delay + 1;

	GenerateSqrtrcFilter(RollOff, Delay, TxNSS);


	
	int specSize, bufSize;
	IppStatus status = ippsFIRMRGetSize(FilterLength, TxNSS, 1, ipp32f, &specSize, &bufSize);
	pSpec = (IppsFIRSpec_32f*)ippsMalloc_8u(specSize);
	pBuf = ippsMalloc_8u(bufSize);
	status = ippsFIRMRInit_32f(SqrtrcFilter, FilterLength, TxNSS, 0, 1, 0, pSpec);
	
	pDlySrcI = ippsMalloc_32f(FilterLength);
	pDlySrcQ = ippsMalloc_32f(FilterLength);
	pDlyDstI = ippsMalloc_32f(FilterLength);
	pDlyDstQ = ippsMalloc_32f(FilterLength);
	
	ippsZero_32f(pDlySrcI,FilterLength);
	ippsZero_32f(pDlySrcQ,FilterLength);
	ippsZero_32f(pDlyDstI,FilterLength);
	ippsZero_32f(pDlyDstQ,FilterLength);


}


void TxFilter::CreateOutputs(Ipp32f *InI, Ipp32f *InQ, Ipp32f *OutI, Ipp32f *OutQ, unsigned int LengthIn)
{
	//Allocate Memory if needed
	

	ippsFIRMR_32f(InI, OutI, LengthIn, pSpec, pDlySrcI, pDlyDstI, pBuf);
	ippsFIRMR_32f(InQ, OutQ, LengthIn, pSpec, pDlySrcQ, pDlyDstQ, pBuf);
	//exchange source destination
	float *Tmp = pDlySrcI;
	pDlySrcI = pDlyDstI;
	pDlyDstI = Tmp;
	Tmp = pDlySrcQ;
	pDlySrcQ = pDlyDstQ;
	pDlyDstQ = Tmp;
	
}
	



void TxFilter::GenerateSqrtrcFilter(double Alpha, unsigned int HalfDelay, Ipp32u SamplingFrequency)
{
	unsigned int FilterLength = 2*HalfDelay + 1;

	Ipp64f *NTs = ippsMalloc_64f(FilterLength);	//Sampling Times
	Ipp64f TS = 1.0/Ipp64f(SamplingFrequency);
	SqrtrcFilter = ippsMalloc_32f(FilterLength); //Allocate Memory for the filter

	//Initialize to -HD*Ts:Ts:HD*Ts;
	ippsVectorSlope_64f(NTs, FilterLength, -1.0 * Ipp64f(HalfDelay) * TS, TS);
	//Zero is in the middle
	SqrtrcFilter[HalfDelay] = (PI*(1-Alpha)+4*Alpha)/PI/sqrt(double(SamplingFrequency));

	for(unsigned int ii =0; ii < HalfDelay; ii++)
	{
		Ipp64f PiT = PI * NTs[ii];
		if(abs(NTs[ii]) == 0.25/Alpha)
		{
			//Special treatment
			SqrtrcFilter[ii] = -1/sqrt(double(SamplingFrequency))/2*(cos(0.25/Alpha*PI*(1-Alpha))*PI*(1-Alpha)+4*Alpha*cos(0.25/Alpha*PI*(1+Alpha))-sin(0.25/Alpha*PI*(1+Alpha))*PI*(1+Alpha))/PI;
		}
		else
		{
			Ipp64f h1 = sin(PiT*(1-Alpha)) + 4*Alpha*NTs[ii]*cos(PiT*(1+Alpha));

			Ipp64f h2 =(PiT*(1-(4*Alpha*NTs[ii])*(4*Alpha*NTs[ii])));

			SqrtrcFilter[ii] = (float) (h1/h2/sqrt(float(SamplingFrequency)));
			SqrtrcFilter[FilterLength - 1 - ii] = SqrtrcFilter[ii];

		}
	}
	/*
	float Sums[2] = {0,0};
	for(int i =0;i < FilterLength; i++)
	{
		Sums[i&1] += 	abs(SqrtrcFilter[i]);
	}

	
	float Mx = Sums[0];
	if(Sums[1]>Mx)
		Mx = Sums[1];
	float kx = 0.95*32767/Mx;
	for(int i =0;i < FilterLength; i++)
	{
		SqrtrcFilter[i]*= kx;
	}
		
*/
	TxPower = 1;
	ippsFree(NTs);



}

