/**
 * @file FrequencyOffset.h
 * @brief Digital complex NCO: 48-bit phase accumulator, LUT cos/sin, in-place multiply on I/Q (@c CreateOutputs).
 *
 * @details @c SetFrequency sets per-sample phase increment (from Hz and @c SamplingFrequency). @c SetAcceleration
 * adds a per-sample change to that increment (ramp). @c CreateOutputs applies @f$(I+jQ)(e^{j\phi})@f$ per sample using
 * gathered LUT values; phase advances each sample. Used by @ref AWGNChannel and @ref ReceiverFreqCorrector.
 */
#pragma once
//#define __INTEL_COMPILER_USE_INTRINSIC_PROTOTYPES 1
#include "definitions.h"
#include <cmath>
#include <cstdint>
#include "ipp.h"
#include "random_generator_new.h"
#define  SCTableLength 16384
#define  Log2SCTableLength 14
#include <immintrin.h>
#define MM __m256i

class FrequencyOffset
{
	

	uint64_t Accumulator; //48 bits
	uint64_t Mask = 281474976710655;//48 bits; 
	uint64_t Delta;
	int64_t DDelta;
	double LastFreqOffset = -1e10;
	unsigned int OutputLength, BufferLength;
	double TimeSample, PhaseFactor,TimeSample2;
	Ipp32f *CosTable, *SinTable;
	double K48;
	double SamplingFrequency;
	void CreateObjects(void);
	void AdvanceAccumulator(void);
	

public:
	unsigned int Phs;
	FrequencyOffset();
	~FrequencyOffset(void);
	void CreateOutputs(Ipp32f *InputI, Ipp32f *InputQ, unsigned int InputLength);//in place double Freq correction
	void SetFrequency(Ipp64f FreqOffsetHz);
	void SetAcceleration(double AccelerationHzSec);
	void SetSamplingFrequency(double val)
	{
		SamplingFrequency = val;
		TimeSample = K48/SamplingFrequency;
		TimeSample2 = K48/(SamplingFrequency*SamplingFrequency);
	
	}
	unsigned int GetOutputLength(void)
	{
		return OutputLength;
	}
	
};

