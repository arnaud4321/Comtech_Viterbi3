#pragma once

/** Decision-directed phase tracking loop. Order-1 IIR on phase error; rotation applied before symbol path. */
class PhaseTracker
{
public:
    PhaseTracker();

    /** Process a block of I,Q (freq-corrected); update phase estimate via decision-directed loop + IIR. */
    void Update(const float* I, const float* Q, int Length);
    /** Current phase estimate (rad); apply rotation by -GetPhaseEst() to correct. */
    float GetPhaseEst() const { return phase_est_; }
    /** True when residual phase error is low (loop converged). */
    bool IsLocked() const { return locked_; }
    /** Filtered phase error (rad) for display. */
    float GetFilteredError() const { return filtered_error_; }
    /** Lock threshold (rad) for display. */
    float GetLockThreshold() const { return lock_threshold_rad_; }

    void SetIirAlpha(float alpha) { alpha_ = alpha; }
    void SetLockThreshold(float rad) { lock_threshold_rad_ = rad; }
    void SetLockCount(int n) { lock_count_required_ = n; }

private:
    float phase_est_ = 0.f;
    float filtered_error_ = 0.f;
    float alpha_ = 0.02f;           /* IIR: filtered = alpha*raw + (1-alpha)*filtered */
    float gain_ = 0.5f;             /* phase_est += gain * filtered_error */
    float lock_threshold_rad_ = 0.05f;
    int lock_count_required_ = 10;
    int lock_count_ = 0;
    bool locked_ = false;
};
