#include "RxFilter.h"
RxFilter::RxFilter(void)
{
}

RxFilter::~RxFilter()
{
}

void RxFilter::CreateObjects(Ipp32f RollOffIn)
{
	//Calculate Length
	RollOff = RollOffIn;
	DeleteObjects();
	unsigned int Delay;
	Delay = ceil(4.0/RollOff);
	FilterLength = 2*Delay + 1;

	GenerateSqrtrcFilter(RollOff, Delay, TxNSS);

    for(int i = 0; i < FilterLength; i++)
        SqrtrcFilter[i] *= K;
	
	int specSize, bufSize;
	IppStatus status = ippsFIRSRGetSize(FilterLength, ipp32f, &specSize, &bufSize);
	pSpec = (IppsFIRSpec_32f*)ippsMalloc_8u(specSize);
	pBuf = ippsMalloc_8u(bufSize);
	status = ippsFIRSRInit_32f(SqrtrcFilter, FilterLength, ippAlgAuto,pSpec);
	
	pDlySrcI = ippsMalloc_32f(FilterLength);
	pDlySrcQ = ippsMalloc_32f(FilterLength);
	pDlyDstI = ippsMalloc_32f(FilterLength);
	pDlyDstQ = ippsMalloc_32f(FilterLength);
	
	ippsZero_32f(pDlySrcI,FilterLength);
	ippsZero_32f(pDlySrcQ,FilterLength);
	ippsZero_32f(pDlyDstI,FilterLength);
	ippsZero_32f(pDlyDstQ,FilterLength);


}


void RxFilter::CreateOutputs(Ipp32f *InI, Ipp32f *InQ, Ipp32f *OutI, Ipp32f *OutQ, unsigned int LengthIn)
{
	//Allocate Memory if needed
	

	ippsFIRSR_32f(InI, OutI, LengthIn, pSpec, pDlySrcI, pDlyDstI, pBuf);
	ippsFIRSR_32f(InQ, OutQ, LengthIn, pSpec, pDlySrcQ, pDlyDstQ, pBuf);
	//exchange source destination
	float *Tmp = pDlySrcI;
	pDlySrcI = pDlyDstI;
	pDlyDstI = Tmp;
	Tmp = pDlySrcQ;
	pDlySrcQ = pDlyDstQ;
	pDlyDstQ = Tmp;
	
}
	

