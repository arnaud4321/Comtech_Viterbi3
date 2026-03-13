#include "PhaseTracker.h"
#include <cmath>

namespace {
/* QPSK slice: decision (dI, dQ) = sign(I), sign(Q) normalized to ±1 */
inline void slice_qpsk(float I, float Q, float& dI, float& dQ)
{
    dI = (I >= 0.f) ? 1.f : -1.f;
    dQ = (Q >= 0.f) ? 1.f : -1.f;
}
/* Phase error from z and decision: angle(z) - angle(d) ≈ (I*dQ - Q*dI) / (|z|*|d|) for small error; use (I*dQ - Q*dI) as sin(phase_err) proxy */
inline float phase_error_sample(float I, float Q, float dI, float dQ)
{
    float sin_err = I * dQ - Q * dI;
    float cos_err = I * dI + Q * dQ;
    if (std::abs(cos_err) < 1e-6f)
        return 0.f;
    return std::atan2(sin_err, cos_err);
}
}

PhaseTracker::PhaseTracker() = default;

void PhaseTracker::Update(const float* I, const float* Q, int Length)
{
    if (!I || !Q || Length <= 0)
        return;

    float c = std::cos(phase_est_);
    float s = std::sin(phase_est_);
    float sum_err = 0.f;
    int count = 0;
    for (int i = 0; i < Length; i++)
    {
        float ii = I[i], qq = Q[i];
        float zI = ii * c + qq * s;
        float zQ = -ii * s + qq * c;
        float dI, dQ;
        slice_qpsk(zI, zQ, dI, dQ);
        sum_err += phase_error_sample(zI, zQ, dI, dQ);
        count++;
    }
    if (count <= 0)
        return;

    float avg_err = sum_err / count;
    filtered_error_ = alpha_ * avg_err + (1.f - alpha_) * filtered_error_;
    phase_est_ += gain_ * filtered_error_;

    if (std::abs(filtered_error_) < lock_threshold_rad_)
    {
        lock_count_++;
        if (lock_count_ >= lock_count_required_)
            locked_ = true;
    }
    else
        lock_count_ = 0;
}
