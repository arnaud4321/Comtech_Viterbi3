/**
 * @file ReceiverPhaseTrackingDD.cpp
 * @brief QPSK decision-directed PLL: slice → @c dd_phase_error → 2nd-order loop → rotate → frame queue.
 *
 * @details **ThreadMain** — @c WaitPopSymbolFrame from @ref ReceiverTimingTracking; per-sample loop matches
 * @ref ReceiverPhaseTrackingDD (LUT derotation, @c dd_phase_error, @c kp_/ @c ki_ updates, @c outI/outQ = derotated
 * @f$z_r@f$ before next symbol’s update). Queues frames for the Viterbi manager.
 */

#include "ReceiverPhaseTrackingDD.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <sys/stat.h>
extern mutex mtxfilethr;

// Define to use std::atan2 (libm reference, comparisons / profiling) instead of the fast
// polynomial in dd_phase_error. Example: -DPHASE_DD_USE_LIBM_ATAN2 on the compiler command line.
// #define PHASE_DD_USE_LIBM_ATAN2

namespace
{
constexpr float kTwoPi = 6.2831853071795864769f;
constexpr float kInvTwoPi = 1.0f / kTwoPi;

inline float wrap_pm_pi(float x)
{
    // Wrap to (-pi, pi]
    x = std::fmod(x, kTwoPi);
    if (x <= -3.14159265358979323846f)
        x += kTwoPi;
    else if (x > 3.14159265358979323846f)
        x -= kTwoPi;
    return x;
}

inline void qpsk_slicer(float i, float q, float& di, float& dq)
{
    di = (i >= 0.0f) ? 1.0f : -1.0f;
    dq = (q >= 0.0f) ? 1.0f : -1.0f;
}

#ifndef PHASE_DD_USE_LIBM_ATAN2
/**
 * @brief Minimax polynomial atan(x) for x in [0,1] (≈ max error ~1e-4 rad), quadrant fix-up.
 *        Cheaper than libm atan2 on the hot PLL path (called every symbol).
 */
inline float fast_atan2(float y, float x)
{
    if (x == 0.0f && y == 0.0f)
        return 0.0f;
    const float ax = std::fabs(x);
    const float ay = std::fabs(y);
    const float mx = (ax > ay) ? ax : ay;
    const float mn = (ax > ay) ? ay : ax;
    const float a = (mn / (mx + 1e-20f));
    const float a2 = a * a;
    float r = (((-0.0464964749f * a2 + 0.15931422f) * a2 - 0.327622764f) * a2 + 0.99997726f) * a;
    if (ay > ax)
        r = 1.57079637f - r;
    if (x < 0.0f)
        r = 3.14159265f - r;
    if (y < 0.0f)
        r = -r;
    return r;
}
#endif

inline float dd_phase_error(float i, float q, float di, float dq)
{
    // e = angle(z * conj(d)) = atan2(Im, Re), where:
    // z = i + j q, d = di + j dq
    // z * conj(d) = (i*di + q*dq) + j(q*di - i*dq)
    const float re = i * di + q * dq;
    const float im = q * di - i * dq;
#ifdef PHASE_DD_USE_LIBM_ATAN2
    return std::atan2(im, re);
#else
    return fast_atan2(im, re);
#endif
}
} // namespace

ReceiverPhaseTrackingDD::ReceiverPhaseTrackingDD()
{
    for (int k = 0; k < kLutSize; ++k)
    {
        const float a = kTwoPi * static_cast<float>(k) / static_cast<float>(kLutSize);
        cosLut_[k] = std::cos(a);
        sinLut_[k] = std::sin(a);
    }
}

ReceiverPhaseTrackingDD::~ReceiverPhaseTrackingDD()
{
    StopJoin();
}

void ReceiverPhaseTrackingDD::Start(ReceiverTimingTracking* timingTracking, std::atomic<bool>* stopAll,
                                    double displayPeriodSec, std::atomic<double>* estSymbolRateHz,
                                    const PhaseTrackingConfig& config)
{
    timingTracking_ = timingTracking;
    stopAll_ = stopAll;
    estSymbolRateHz_ = estSymbolRateHz;
    displayPeriodSec_ = displayPeriodSec;
    phaseRad_ = 0.0f;
    freqRadPerSym_ = 0.0f;
    errEma_ = 0.0f;
    lockCount_ = 0;
    locked_.store(false, std::memory_order_relaxed);
    kp_ = static_cast<float>(config.Kp);
    ki_ = static_cast<float>(config.Ki);

    threadRunning_ = true;
    thread_ = std::thread(&ReceiverPhaseTrackingDD::ThreadMain, this);
}

void ReceiverPhaseTrackingDD::StopJoin()
{
    if (!threadRunning_)
        return;
    threadRunning_ = false;
    cvFrameReady_.notify_all();
    cvFrameSpace_.notify_all();
    if (thread_.joinable())
        thread_.join();

#ifdef DEBUG_PHASE_DD
    closeDebugFiles();
#endif
}

bool ReceiverPhaseTrackingDD::pushOneFrameToQueue(const float* i, const float* q)
{
    std::unique_lock<std::mutex> lk(mtxFrameQ_);
    cvFrameSpace_.wait(lk, [&] { return !threadRunning_ || frameQCount_ < kFrameQueueDepth; });
    if (!threadRunning_)
        return false;
    FrameSlot& slot = framePool_[frameQTail_];
    std::copy(i, i + kSymFrame, slot.i);
    std::copy(q, q + kSymFrame, slot.q);
    frameQTail_ = (frameQTail_ + 1) % kFrameQueueDepth;
    frameQCount_++;
    lk.unlock();
    cvFrameReady_.notify_one();
    return true;
}

bool ReceiverPhaseTrackingDD::WaitPopSymbolFrame(float* dstI, float* dstQ, int nSym)
{
    if (nSym != kSymFrame)
        return false;
    std::unique_lock<std::mutex> lk(mtxFrameQ_);
    cvFrameReady_.wait(lk, [&] { return !threadRunning_ || frameQCount_ > 0; });
    if (frameQCount_ <= 0)
        return false;
    FrameSlot& slot = framePool_[frameQHead_];
    std::copy(slot.i, slot.i + kSymFrame, dstI);
    std::copy(slot.q, slot.q + kSymFrame, dstQ);
    frameQHead_ = (frameQHead_ + 1) % kFrameQueueDepth;
    frameQCount_--;
    lk.unlock();
    cvFrameSpace_.notify_one();
    return true;
}

#ifdef DEBUG_PHASE_DD
void ReceiverPhaseTrackingDD::openDebugFilesIfNeeded()
{
    if (dbgIn_)
        return;
    mkdir("../data", 0755);
    dbgIn_ = std::fopen("../data/phase_dd_in.bin", "wb");
    dbgOut_ = std::fopen("../data/phase_dd_out.bin", "wb");
    dbgPhase_ = std::fopen("../data/phase_dd_phase.bin", "wb");
    dbgErr_ = std::fopen("../data/phase_dd_err.bin", "wb");
    dbgErrFilt_ = std::fopen("../data/phase_dd_err_filt.bin", "wb");
    dbgThresh_ = std::fopen("../data/phase_dd_thresh.bin", "wb");
}

void ReceiverPhaseTrackingDD::closeDebugFiles()
{
    if (dbgIn_)
    {
        std::fclose(dbgIn_);
        dbgIn_ = nullptr;
    }
    if (dbgOut_)
    {
        std::fclose(dbgOut_);
        dbgOut_ = nullptr;
    }
    if (dbgPhase_)
    {
        std::fclose(dbgPhase_);
        dbgPhase_ = nullptr;
    }
    if (dbgErr_)
    {
        std::fclose(dbgErr_);
        dbgErr_ = nullptr;
    }
    if (dbgErrFilt_)
    {
        std::fclose(dbgErrFilt_);
        dbgErrFilt_ = nullptr;
    }
    if (dbgThresh_)
    {
        std::fclose(dbgThresh_);
        dbgThresh_ = nullptr;
    }
}
#endif

/**
 * @brief Phase tracking loop consuming Gardner frames, producing corrected symbol frames.
 */
void ReceiverPhaseTrackingDD::ThreadMain()
{

    #ifdef WRITE_LOG_THR
	mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"Receiver Phase Tracking  Thread %d\n", gettid());
    fclose(fidthr);
    mtxfilethr.unlock();
    #endif

    alignas(32) float inI[kSymFrame];
    alignas(32) float inQ[kSymFrame];
    alignas(32) float outI[kSymFrame];
    alignas(32) float outQ[kSymFrame];
    uint64_t symbols_total = 0;
    const auto wall_start = std::chrono::steady_clock::now();
    bool prevLocked = false;
    auto lastFreqDisplay = wall_start;

    while (threadRunning_ && (stopAll_ == nullptr || !*stopAll_))
    {
        if (!timingTracking_ || !timingTracking_->WaitPopSymbolFrame(inI, inQ, kSymFrame))
            break;

#ifdef DEBUG_PHASE_DD
        openDebugFilesIfNeeded();
#endif

        // When locked, update the DD-PLL only once every N symbols to reduce CPU.
        // Gains are scaled by N on update symbols to roughly preserve loop bandwidth.
        constexpr int kLockedUpdateDecim = 4;
        bool lockedLocal = locked_.load(std::memory_order_relaxed);

        float sumMagSq = 0.0f;
        float sumMag = 0.0f;

        for (int n = 0; n < kSymFrame; ++n)
        {
            // Apply current NCO rotation: z_rot = z * exp(-j*phase)
            const float phase = phaseRad_;
            int idx = static_cast<int>(std::floor((phase * kInvTwoPi) * static_cast<float>(kLutSize)));
            idx %= kLutSize;
            if (idx < 0)
                idx += kLutSize;
            const float c = cosLut_[idx];
            const float s = sinLut_[idx];

            const float i = inI[n];
            const float q = inQ[n];
            const float zrI = i * c + q * s;
            const float zrQ = -i * s + q * c;

            float magI = std::abs(zrI);
            float magQ = std::abs(zrQ);
            sumMagSq += (magI * magI + magQ * magQ);
            sumMag += (magI + magQ);

            float di, dq;
            qpsk_slicer(zrI, zrQ, di, dq);

            float err = 0.0f;
            const bool doUpdate = (!lockedLocal) || (kLockedUpdateDecim <= 1) || ((n % kLockedUpdateDecim) == 0);
            if (doUpdate)
            {
                err = dd_phase_error(zrI, zrQ, di, dq);

                // 2nd-order PLL update (symbol-rate loop)
                freqRadPerSym_ += ki_ * err;
                phaseRad_ += freqRadPerSym_ + kp_ * err;
                phaseRad_ = wrap_pm_pi(phaseRad_);

                // Lock detection on smoothed |err| (only needed until lock is reached).
                if (!lockedLocal)
                {
                    const float aerr = std::abs(err);
                    errEma_ = errEmaAlpha_ * aerr + (1.0f - errEmaAlpha_) * errEma_;
                    if (errEma_ < lockThresholdRad_)
                        lockCount_++;
                    else
                        lockCount_ = 0;
                    if (lockCount_ >= lockCountRequired_)
                    {
                        locked_.store(true, std::memory_order_relaxed);
                        lockedLocal = true;
                    }
                }
            }
            else
            {
                // When skipping the DD update, still advance the phase using current frequency estimate.
                phaseRad_ += freqRadPerSym_;
                phaseRad_ = wrap_pm_pi(phaseRad_);
            }

            // Output corrected sample (already rotated)
            outI[n] = zrI;
            outQ[n] = zrQ;

#ifdef DEBUG_PHASE_DD
            if (dbgIn_)
            {
                float iqIn[2] = {i, q};
                std::fwrite(iqIn, sizeof(float), 2, dbgIn_);
            }
            if (dbgOut_)
            {
                float iqOut[2] = {zrI, zrQ};
                std::fwrite(iqOut, sizeof(float), 2, dbgOut_);
            }
            if (dbgPhase_)
            {
                const float p = phaseRad_;
                std::fwrite(&p, sizeof(float), 1, dbgPhase_);
            }
            if (dbgErr_)
            {
                std::fwrite(&err, sizeof(float), 1, dbgErr_);
            }
            if (dbgErrFilt_)
            {
                const float e = errEma_;
                std::fwrite(&e, sizeof(float), 1, dbgErrFilt_);
            }
            if (dbgThresh_)
            {
                const float th = lockThresholdRad_;
                std::fwrite(&th, sizeof(float), 1, dbgThresh_);
            }
#endif
        }

        float blockEvmRms = 0.0f;
        if (sumMagSq > 1e-12f && sumMag > 1e-12f)
        {
            float P = sumMagSq / static_cast<float>(kSymFrame);
            float M = sumMag / static_cast<float>(kSymFrame);
            // Classical M2M4-based estimator for QPSK: EVM^2 = (2*P / M^2) - 1
            float evmSq = (2.0f * P) / (M * M) - 1.0f;
            if (evmSq > 0.0f)
                blockEvmRms = std::sqrt(evmSq);
        }
        lastEvmRms_.store(static_cast<double>(blockEvmRms), std::memory_order_relaxed);

        symbols_total += static_cast<uint64_t>(kSymFrame);
        const double rsDenom = (estSymbolRateHz_ != nullptr)
                                   ? ([&]() {
                                           const double e =
                                               estSymbolRateHz_->load(std::memory_order_relaxed);
                                           return (e > 1.0) ? e : SymbolRate;
                                       }())
                                   : SymbolRate;
        const bool nowLocked = locked_.load(std::memory_order_relaxed);
        if (nowLocked != prevLocked)
        {
            const double t_rate = static_cast<double>(symbols_total) / rsDenom;
            const double t_sim =
                std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
            if (displayPeriodSec_ > 0.0)
            {
                std::cout << "[ReceiverPhaseTrackingDD] "
                          << (nowLocked ? "\033[32mLOCKED\033[0m" : "\033[31mUNLOCKED\033[0m")
                          << " t_sim=" << t_sim << " s"
                          << " t_rate=" << t_rate << " s"
                          << " errEma=" << errEma_ << " rad"
                          << " thresh=" << lockThresholdRad_ << " rad"
                          << std::endl;
            }
            prevLocked = nowLocked;
        }

        // Display residual frequency estimate once per second (from phase slope / PLL freq state).
        {
            const auto now = std::chrono::steady_clock::now();
            if (displayPeriodSec_ > 0.0 &&
                now - lastFreqDisplay >= std::chrono::duration<double>(displayPeriodSec_))
            {
                // freqRadPerSym_ is radians per symbol; convert to Hz using symbol rate.
                const double f_hz = static_cast<double>(freqRadPerSym_) * SymbolRate / (2.0 * 3.14159265358979323846);
                lastFreqEstHz_.store(f_hz, std::memory_order_relaxed);
                const double t_rate = static_cast<double>(symbols_total) / rsDenom;
                const double t_sim = std::chrono::duration<double>(now - wall_start).count();
                const bool nowLocked2 = locked_.load(std::memory_order_relaxed);
                double snrEvmDb = (blockEvmRms > 1e-12f) ? (-20.0 * std::log10(blockEvmRms)) : std::numeric_limits<double>::quiet_NaN();
                std::cout << "[ReceiverPhaseTrackingDD] freq_est=" << f_hz << " Hz"
                          << " locked=" << (nowLocked2 ? 1 : 0)
                          << " snrEvmDb=" << snrEvmDb
                          << " errEma=" << errEma_ << " rad"
                          << " thresh=" << lockThresholdRad_ << " rad"
                          << " t_sim=" << t_sim << " s"
                          << " t_rate=" << t_rate << " s"
                          << std::endl;
                lastFreqDisplay = now;
            }
        }

        // Always feed downstream, even when not locked.
        // Consumers can check IsLocked() to interpret status; discarding frames here can stall the pipeline.
        if (!pushOneFrameToQueue(outI, outQ))
            break;
    }
}

