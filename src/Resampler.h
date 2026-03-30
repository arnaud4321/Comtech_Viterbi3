#pragma once
#define __INTEL_COMPILER_USE_INTRINSIC_PROTOTYPES 1
#define Log2LagrangeTableLength 12
#define LagrangeTableLength (1<<Log2LagrangeTableLength)
#define ACCUMULATOR_WIDTH 30
#include <immintrin.h>
#include <math.h>
// Note: avoid short type macros (MM/MM2) to prevent collisions with other headers.

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
	/// Fractional-delay resampling using cubic Lagrange interpolation.
	/// - Advance: input samples consumed per output sample (>= 0.5).
	/// - Start: input index cursor (advanced on return).
	/// - Frac: fractional phase in [0,1) (updated on return).
	/// - Length: available input length (must have Start < Length-3 to interpolate).
	/// - outMax: max output samples to generate (prevents output buffer overflow).
	unsigned int CreateOutputs(const float *InI, const float *InQ,
	                          float *OutputI, float *OutputQ,
	                          double Advance, unsigned int &Start, double &Frac,
	                          unsigned int Length, unsigned int outMax);

};

