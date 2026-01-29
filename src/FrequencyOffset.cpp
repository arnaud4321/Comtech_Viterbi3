#include "FrequencyOffset.h"
#include <cstdint>



FrequencyOffset::FrequencyOffset()
{
	
	
	OutputLength = 0; 
	BufferLength = 0;
	K48 = 281474976710656.0;
	
	CosTable = ippsMalloc_32f(SCTableLength);
	SinTable = ippsMalloc_32f(SCTableLength);
	PhaseFactor = 1.0 / IPP_2PI*TimeSample;
	CreateObjects();
}


FrequencyOffset::~FrequencyOffset(void)
{
	ippsFree(CosTable);
	ippsFree(SinTable);

}

void FrequencyOffset::CreateObjects()
{
	double Pi = 4*atan(1.0);
	double Ts0 = 1.0 / double(SCTableLength);

	Ipp64f Phase0 = Ts0*2.0*Pi;
	for(unsigned int ii = 0; ii < SCTableLength; ii++)

	{
		double Phase = Phase0 * double(ii);
		CosTable[ii] = cos(Phase);
		SinTable[ii] = sin(Phase);
	}
}



void FrequencyOffset::SetFrequency(Ipp64f FreqOffset)
{
	double DeltaFreq = FreqOffset - LastFreqOffset;
	DeltaFreq = abs(DeltaFreq);
	if (DeltaFreq > 1e-8)
	{
		LastFreqOffset = FreqOffset;
		Ipp64f dblDelta = FreqOffset*TimeSample; //Frequency Offset * Ts  {1/(Bit Rate * Tx NSS )}
		int64_t iDelta = floor(dblDelta + 0.5);
		if (Delta < 0)
		{
			iDelta += 281474976710656;
		}
		Delta = (uint64_t) iDelta;
	}
}
void FrequencyOffset::SetAcceleration(double Acceleration)
{
	DDelta = (int64_t) (TimeSample2*Acceleration+0.5);
}




void FrequencyOffset::AdvanceAccumulator()
{
	Accumulator += Delta;
	Accumulator &= Mask;
	Delta += DDelta;
	Delta &= Mask;
}



void FrequencyOffset::CreateOutputs(Ipp32f *InputI, Ipp32f *InputQ, unsigned int InputLength)
{


	alignas(32) uint32_t Indexes[8];
	for (unsigned int ii = 0; ii < InputLength; ii += 8)
	{
		__m256 xi = _mm256_load_ps(&InputI[ii]);
		__m256 xq = _mm256_load_ps(&InputQ[ii]);

		for(int jj = 0; jj < 8; jj++)
		{
			Indexes[jj] = (uint32_t) (Accumulator>>(48-Log2SCTableLength));
			AdvanceAccumulator();
		}

		MM mInd = _mm256_load_si256((MM*)Indexes);
		__m256 cv = _mm256_i32gather_ps(CosTable, mInd, 4);
		__m256 sv = _mm256_i32gather_ps(SinTable, mInd, 4);



		//(xi+i xq)(cv+i sv) = (xi cv - xq sv) + i(xq cv + xi sv)

		__m256 ic = _mm256_mul_ps(xi, cv);
		__m256 qc = _mm256_mul_ps(xq, cv);


		__m256 outi = _mm256_fnmadd_ps(xq, sv, ic);
		__m256 outq = _mm256_fmadd_ps(xi, sv, qc);



		_mm256_storeu_ps((Ipp32f *)&InputI[ii], outi);
		_mm256_storeu_ps((Ipp32f *)&InputQ[ii], outq);


	}

}
