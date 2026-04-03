/**
 * @file Lagrange4Simd.h
 * @brief SIMD (SSE/AVX) helpers for cubic Lagrange-4 fractional interpolation.
 */
#pragma once

#include <cmath>
#include <immintrin.h>

namespace Lagrange4Simd
{

/** @brief Cubic Lagrange weights for fractional offset @a mu; matches samples x0..x3 (k-1 .. k+2). */
inline __m128 Coeffs_ps(double mu)
{
    const float c0 = static_cast<float>(-mu * (mu - 1.0) * (mu - 2.0) / 6.0);
    const float c1 = static_cast<float>((mu + 1.0) * (mu - 1.0) * (mu - 2.0) / 2.0);
    const float c2 = static_cast<float>(-(mu + 1.0) * mu * (mu - 2.0) / 2.0);
    const float c3 = static_cast<float>((mu + 1.0) * mu * (mu - 1.0) / 6.0);
    return _mm_set_ps(c3, c2, c1, c0);
}

inline float HSum128_ps(__m128 v)
{
    __m128 shuf = _mm_movehdup_ps(v);
    __m128 sums = _mm_add_ps(v, shuf);
    __m128 shuf2 = _mm_movehl_ps(shuf, sums);
    sums = _mm_add_ss(sums, shuf2);
    return _mm_cvtss_f32(sums);
}

inline float Dot_ps(__m128 samples_x0_to_x3, __m128 coeffs)
{
    return HSum128_ps(_mm_mul_ps(samples_x0_to_x3, coeffs));
}

/** @brief @a x_at_k_minus_1 points to four contiguous samples: x[k-1]..x[k+2]. */
inline float EvalContiguous_ps(const float* x_at_k_minus_1, __m128 coeffs)
{
    return Dot_ps(_mm_loadu_ps(x_at_k_minus_1), coeffs);
}

inline float EvalSet_ps(float x0, float x1, float x2, float x3, __m128 coeffs)
{
    return Dot_ps(_mm_set_ps(x3, x2, x1, x0), coeffs);
}

/**
 * @brief Same @a mu for I and Q: one coefficient build, 256-bit multiply for both dot products (AVX).
 */
inline void EvalIQ_ps(float i0, float i1, float i2, float i3, float q0, float q1, float q2, float q3, double mu,
                      float* outI, float* outQ)
{
    const __m128 c = Coeffs_ps(mu);
#if defined(__AVX__)
    const __m256 cc = _mm256_set_m128(c, c);
    const __m256 v = _mm256_set_ps(q3, q2, q1, q0, i3, i2, i1, i0);
    const __m256 p = _mm256_mul_ps(v, cc);
    *outI = HSum128_ps(_mm256_castps256_ps128(p));
    *outQ = HSum128_ps(_mm256_extractf128_ps(p, 1));
#else
    *outI = EvalSet_ps(i0, i1, i2, i3, c);
    *outQ = EvalSet_ps(q0, q1, q2, q3, c);
#endif
}

/** @brief Contiguous I/Q taps at @a i_k_minus_1 and @a q_k_minus_1 (four floats each). */
inline void EvalIQContiguous_ps(const float* i_k_minus_1, const float* q_k_minus_1, double mu, float* outI,
                                float* outQ)
{
    const __m128 c = Coeffs_ps(mu);
#if defined(__AVX__)
    const __m256 cc = _mm256_set_m128(c, c);
    const __m256 v = _mm256_set_m128(_mm_loadu_ps(q_k_minus_1), _mm_loadu_ps(i_k_minus_1));
    const __m256 p = _mm256_mul_ps(v, cc);
    *outI = HSum128_ps(_mm256_castps256_ps128(p));
    *outQ = HSum128_ps(_mm256_extractf128_ps(p, 1));
#else
    *outI = EvalContiguous_ps(i_k_minus_1, c);
    *outQ = EvalContiguous_ps(q_k_minus_1, c);
#endif
}

} // namespace Lagrange4Simd
