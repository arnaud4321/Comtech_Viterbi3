/**
 * @file Resampler.h
 * @brief Cubic Lagrange fractional interpolator used by @ref ReceiverResampler.
 */
#pragma once
#define __INTEL_COMPILER_USE_INTRINSIC_PROTOTYPES 1
#define Log2LagrangeTableLength 12
#define LagrangeTableLength (1<<Log2LagrangeTableLength)
#define ACCUMULATOR_WIDTH 30
#include <immintrin.h>
#include <math.h>
// Note: avoid short type macros (MM/MM2) to prevent collisions with other headers.

/**
 * @brief Precomputed cubic Lagrange coefficients and fixed-point accumulator for variable-rate resampling.
 *
 * @details For each output sample, four consecutive input samples are weighted from @c LagrangeTable rows
 * indexed by the high bits of a fixed-point @c Accumulator (phase); @c Advance is quantized to @c Delta per
 * @c CreateOutputs call, and @c AdvanceAccumulator adds @c Delta then advances the input index @c Start by the
 * integer overflow of @c Accumulator — average input samples per output ≈ @c Advance. This yields variable
 * rate conversion without an explicit polyphase bank. @c CreateOutputs updates @c Start and @c Frac for
 * continuity across calls (see @ref ReceiverResampler).
 */
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

