#pragma once
#define __INTEL_COMPILER_USE_INTRINSIC_PROTOTYPES 1
#define Log2LagrangeTableLength 12
#define LagrangeTableLength (1<<Log2LagrangeTableLength)
#define ACCUMULATOR_WIDTH 30
#include <immintrin.h>
#include <math.h>
#include <immintrin.h>
#define MM __m128
#define MM2 __m256

class Resampler
{
	unsigned int Accumulator, Delta, MaskAccumulator2;
	double PowAccumulator, PowAccumulator2;
	float LagrangeTable[LagrangeTableLength << 2];
	void CreateObjects(void);
	void AdvanceAccumulator(unsigned int &Ptr);
public:
	Resampler();
	~Resampler();
	unsigned int CreateOutputs(float *InI, float *InQ, float *OutputI, float *OutputQ, double Advance, unsigned int &Start, double Frac, unsigned int Length);

};

