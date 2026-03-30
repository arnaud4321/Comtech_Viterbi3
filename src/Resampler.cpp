#include "Resampler.h"



Resampler::Resampler()
{
	unsigned int Tmp = (1 << ACCUMULATOR_WIDTH);
	MaskAccumulator2 =  Tmp - 1;
	PowAccumulator = double(Tmp);
	PowAccumulator2 = 1.0 / PowAccumulator;
	CreateObjects();
}


Resampler::~Resampler()
{
}


void Resampler::CreateObjects()
{
	
	double Ts0 = 1.0 / double(LagrangeTableLength);



	double Ts = 0;
	for (int jj = 0; jj < LagrangeTableLength; jj++)
	{
		double Ts2 = Ts*Ts;
		double Ts3 = Ts2*Ts;
		double Temp = Ts3 - Ts;
		double Temp2 = Temp / 6.0;
		LagrangeTable[jj * 4 + 3] = Temp2;
		//-Ts3/2+Ts2/2+Ts;
		Temp2 = (Ts2-Ts3) *0.5 + Ts;
		LagrangeTable[jj * 4 + 2] = Temp2;
		//Ts3/2-Ts2-Ts/2+1;
		Temp2 = (Ts3-Ts)*0.5-Ts2+1.0;
		LagrangeTable[jj * 4 + 1] = Temp2;
		// -Ts3/6+Ts2/2-Ts/3
		Temp2 = Ts2*0.5 - (Ts3 / 6.0 + Ts / 3.0);
		LagrangeTable[jj * 4 ] = Temp2;
		Ts += Ts0;
	}

}




unsigned int Resampler::CreateOutputs(const float *InI, const float *InQ,
                                      float *OutputI, float *OutputQ,
                                      double Advance, unsigned int &Start, double &Frac,
                                      unsigned int Length, unsigned int outMax)
{
	if ((Advance < 0.5))
	{
		return 0;
	}
	Accumulator = floor((Frac * PowAccumulator)+0.5);
	Delta = floor((Advance * PowAccumulator) + 0.5);
	unsigned int PtrOut = 0;

	

	while (Start < (Length-3) && PtrOut < outMax)
	{
		//Prepare Table Index
		unsigned int TableIndex = (Accumulator >> (ACCUMULATOR_WIDTH - Log2LagrangeTableLength)); //There are 4 coefficients per row
		TableIndex = (TableIndex << 2);
		__m128 mInputI = _mm_loadu_ps(&InI[Start]);
		__m128 mInputQ = _mm_loadu_ps(&InQ[Start]);
		__m256 mInputIQ = _mm256_castps128_ps256(mInputI);
		mInputIQ = _mm256_insertf128_ps(mInputIQ, mInputQ, 1);
		//Load Coefficients
		__m128 Coefficients = _mm_loadu_ps(&LagrangeTable[TableIndex]);
		__m256 Coefficients2 = _mm256_castps128_ps256(Coefficients);
		Coefficients2 = _mm256_insertf128_ps(Coefficients2, Coefficients, 1);
		mInputIQ = _mm256_mul_ps(mInputIQ, Coefficients2); //I0,I1,I2,I3,Q0,Q1,Q2,Q3
		__m256 Results = _mm256_hadd_ps(mInputIQ, mInputIQ);//I01,I23,I01,I23, Q01,Q23,Q01,Q23
		__m128 ResultsIa = _mm256_castps256_ps128(Results);//I01,I23,I01,I23
		__m128 ResultsQa = _mm256_extractf128_ps(Results,1);//Q01,Q23,Q01,Q23

		__m128 Results0 = _mm_hadd_ps(ResultsIa,ResultsQa);//I,I,Q,Q
		_mm_store_ss(OutputI+PtrOut,Results0);
		Results0 = _mm_castsi128_ps( _mm_srli_si128(_mm_castps_si128(Results0), 8));
		_mm_store_ss(OutputQ + PtrOut, Results0);

		PtrOut++;
		//Advance Accumulator
		AdvanceAccumulator(Start);


	}

	Frac = double(Accumulator)*PowAccumulator2;
	return(PtrOut);
}


void Resampler::AdvanceAccumulator(unsigned int &Ptr)
{
	
	
	Accumulator += Delta;
	unsigned int One = Accumulator >> ACCUMULATOR_WIDTH;
	Ptr += One;
	Accumulator &= MaskAccumulator2;

}
