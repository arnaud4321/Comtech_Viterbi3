/**
 * @file GardnerTiming.cpp
 * @brief Oversampled ring, Lagrange-4 interpolation, Gardner TED, PI loop on @f$\omega@f$, lock from @f$\omega@f$ stability.
 *
 * @details **ProcessBlock** — @c PushToRing then loop: one TED term and one on-time output per symbol.
 * @f$\omega@f$ / PI run only every @c updatePeriod_ **output symbols** (~ @f$\mathrm{updatePeriod\_}\,\omega/F_s@f$ s apart).
 * See @ref GardnerTiming.
 */

#include "GardnerTiming.h"
#include "Lagrange4Simd.h"
#include "definitions.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <sys/stat.h>
#include <vector>
#ifdef DEBUG_GARDNER_DEBIT
#include <cstdint>
#include <iostream>
#endif


#ifdef DEBUG_GARDNER_DEBIT
namespace
{
struct GardnerDebitStats
{
    std::uint64_t blocks = 0;
    std::uint64_t sumIn = 0;
    std::uint64_t sumOut = 0;
};
GardnerDebitStats gGardnerDebit;
} // namespace
#endif

GardnerTiming::GardnerTiming()
{
    ringI_.assign(static_cast<size_t>(kRingSize), 0.0f);
    ringQ_.assign(static_cast<size_t>(kRingSize), 0.0f);
    Reset();
}

GardnerTiming::~GardnerTiming()
{
#ifdef DEBUG_GARDNER_OUTPUTS
    CloseDebugFiles();
#endif
}

void GardnerTiming::Reset(double omegaNom, double kp, double ki, int updatePeriod)
{
    omegaNom_ = omegaNom;
    omega_ = omegaNom;
    kp_ = kp;
    ki_ = ki;
    integ_ = 0.0;
    updatePeriod_ = std::max(1, updatePeriod);
    updateCounter_ = 0;
    absWrite_ = 0;
    tAbs_ = 0.0;
    tInit_ = false;

    kOmegaLockSpanThreshold=100.0*ki; // heuristic..

    const double omegaUpdatesPerSec = SymbolRate / static_cast<double>(updatePeriod_);
    const double desiredStride =
        (kOmegaLockWindowSeconds * omegaUpdatesPerSec) / static_cast<double>(kOmegaLockMeasures);
    omegaLockStrideUpdates_ = std::max(1, static_cast<int>(std::llround(desiredStride)));
    omegaLockStrideCounter_ = 0;

    historyOmega_.assign(static_cast<size_t>(kOmegaLockMeasures), omegaNom_);
    historyOmegaWr_ = 0;
    historyOmegaFill_ = 0;
    locked_ = false;
    std::fill(ringI_.begin(), ringI_.end(), 0.0f);
    std::fill(ringQ_.begin(), ringQ_.end(), 0.0f);
#ifdef DEBUG_GARDNER_DEBIT
    gGardnerDebit = GardnerDebitStats{};
#endif
}

/**
 * @brief Ingest @a inLen samples at 2 sps; write up to @a outMax one-sps symbols to @a outI / @a outQ.
 */
int GardnerTiming::ProcessBlock(const float* inI, const float* inQ, int inLen,
                     float* outI, float* outQ, int outMax)
{
    if (!inI || !inQ || !outI || !outQ || inLen < 8 || outMax <= 0)
        return 0;

    PushToRing(inI, inQ, inLen);

    const double oldestSafe = static_cast<double>(OldestAbs() + 2);
    const double newestSafeNow = NewestSafe();
    const double targetLag = static_cast<double>(kRingSize) / 4.0;
    assert(targetLag + 16.0 < static_cast<double>(kRingSize));

    // Strict no-clamp: initialize only when enough history exists.
    if (!tInit_)
    {
        if (absWrite_ > static_cast<long long>(targetLag) + 4)
        {
            tAbs_ = static_cast<double>(absWrite_) - targetLag;
            tInit_ = true;
        }
        else
        {
            return 0;
        }
    }

    // Strict no-clamp: if cursor exits valid interpolation region backwards, re-arm.
    // If it exits forwards (tAbs_ > newestSafeNow), it's normal (consumed all data), just break later.
    if (tAbs_ < oldestSafe)
    {
        tInit_ = false;
        return 0;
    }

    const double tLoopStart = tAbs_;
    int outCount = 0;

#ifdef DEBUG_GARDNER_OUTPUTS
    OpenDebugFilesIfNeeded();
    std::vector<float> dbgEarly;
    std::vector<float> dbgOnTime;
    std::vector<float> dbgLate;
    std::vector<double> dbgErr;
    std::vector<float> dbgOmega;
    std::vector<float> dbgPhi;
    dbgEarly.reserve(static_cast<size_t>(outMax) * 2u);
    dbgOnTime.reserve(static_cast<size_t>(outMax) * 2u);
    dbgLate.reserve(static_cast<size_t>(outMax) * 2u);
    dbgErr.reserve(static_cast<size_t>(outMax));
    dbgOmega.reserve(static_cast<size_t>(outMax));
    dbgPhi.reserve(static_cast<size_t>(outMax));
#endif

    while (outCount < outMax)
    {
        // When locked, compute the Gardner TED / PI update only once every N output symbols to reduce CPU.
        // In between, keep omega_ constant and only interpolate the on-time sample.
        constexpr int kLockedTedDecim = 4;
        const bool lockedLocal = locked_;
        const bool doTed = (!lockedLocal) || (kLockedTedDecim <= 1) || ((outCount % kLockedTedDecim) == 0);

        if (!CanInterp(tAbs_))
            break;

        float yI, yQ;
        InterpLagrange4IQ(tAbs_, &yI, &yQ);

        float eI = 0.0f, eQ = 0.0f, lI = 0.0f, lQ = 0.0f;
        if (doTed)
        {
            // Need availability for t-0.5 and t+0.5 only when the TED is evaluated.
            if (!CanInterp(tAbs_ - 0.5) || !CanInterp(tAbs_ + 0.5))
                break;
            InterpLagrange4IQ(tAbs_ - 0.5, &eI, &eQ);
            InterpLagrange4IQ(tAbs_ + 0.5, &lI, &lQ);
        }

        outI[outCount] = yI;
        outQ[outCount] = yQ;

#ifdef DEBUG_GARDNER_OUTPUTS
        dbgEarly.push_back(eI);
        dbgEarly.push_back(eQ);
        dbgOnTime.push_back(yI);
        dbgOnTime.push_back(yQ); 
        dbgLate.push_back(lI);
        dbgLate.push_back(lQ);
        // 8 samples per symbol at Gardner output (interpolated along one symbol length omega_).
        if (fidOut8Sps_)
        {
            constexpr int kOutSps = 8;
            const double w = omega_;
            double t8[kOutSps];
            bool ok8 = true;
            for (int k = 0; k < kOutSps; ++k)
            {
                t8[k] = tAbs_ - 0.5 * w + (static_cast<double>(k) + 0.5) * (w / static_cast<double>(kOutSps));
                if (!CanInterp(t8[k]))
                    ok8 = false;
            }
            if (ok8)
            {
                for (int k = 0; k < kOutSps; ++k)
                {
                    float ii, qq;
                    InterpLagrange4IQ(t8[k], &ii, &qq);
                    float iq[2] = {ii, qq};
                    std::fwrite(iq, sizeof(float), 2, static_cast<FILE*>(fidOut8Sps_));
                }
            }
        }
#endif

        if (doTed)
        {
            // Gardner TED on complex signal (use current on-time symbol).
            errorAcc_ += (lI - eI) * yI + (lQ - eQ) * yQ;

#ifdef DEBUG_GARDNER_OUTPUTS
            dbgErr.push_back(errorAcc_);
            dbgOmega.push_back(static_cast<float>(omega_));
            // Phase within symbol: phi = frac(t/2) in [0,1).
            const double u = 0.5 * tAbs_;
            const double phi = u - std::floor(u);
            dbgPhi.push_back(static_cast<float>(phi));
#endif

            updateCounter_++;
            if (updateCounter_ >= updatePeriod_)
            {
                updateCounter_ = 0;

                double meanErr = errorAcc_ / updatePeriod_;
                integ_ += ki_ * meanErr;

                omega_ = omegaNom_ + integ_ + kp_ * static_cast<double>(meanErr);
                omega_ = std::max(1.7, std::min(2.3, omega_));
                errorAcc_ = 0.0;
                PushOmegaHistoryOnUpdate();
            }
        }

        tAbs_ += omega_;
        outCount++;
    }

#ifdef DEBUG_GARDNER_OUTPUTS
    if (fidEarly_ && fidOnTime_ && fidLate_ && fidErr_ && fidOmega_ && fidPhi_ && !dbgEarly.empty())
    {
        const size_t nfloat = dbgEarly.size();
        std::fwrite(dbgEarly.data(), sizeof(float), nfloat, static_cast<FILE*>(fidEarly_));
        std::fwrite(dbgOnTime.data(), sizeof(float), nfloat, static_cast<FILE*>(fidOnTime_));
        std::fwrite(dbgLate.data(), sizeof(float), nfloat, static_cast<FILE*>(fidLate_));
        std::fwrite(dbgErr.data(), sizeof(double), static_cast<size_t>(outCount), static_cast<FILE*>(fidErr_));
        std::fwrite(dbgOmega.data(), sizeof(float), static_cast<size_t>(outCount), static_cast<FILE*>(fidOmega_));
        std::fwrite(dbgPhi.data(), sizeof(float), static_cast<size_t>(outCount), static_cast<FILE*>(fidPhi_));
    }
    if (fidOut8Sps_)
        std::fflush(static_cast<FILE*>(fidOut8Sps_));
#endif

    const double tLoopEnd = tAbs_;

#ifdef DEBUG_GARDNER_DEBIT
    const bool cursorClamped = false;
    gGardnerDebit.blocks++;
    gGardnerDebit.sumIn += static_cast<std::uint64_t>(inLen);
    gGardnerDebit.sumOut += static_cast<std::uint64_t>(outCount);

    constexpr std::uint64_t kPrintEvery = 200;
    if (gGardnerDebit.blocks % kPrintEvery == 0)
    {
        const double ratio = (inLen > 0) ? (2.0 * static_cast<double>(outCount) / static_cast<double>(inLen)) : 0.0;
        const double deficitBlk =
            static_cast<double>(outCount) - 0.5 * static_cast<double>(inLen);
        const double globalErr =
            static_cast<double>(gGardnerDebit.sumOut) - 0.5 * static_cast<double>(gGardnerDebit.sumIn);
        std::cout << "[Gardner DEBIT] blk=" << gGardnerDebit.blocks << " inLen=" << inLen << " nSym=" << outCount
                  << " ratio(2*nSym/inLen)=" << ratio << " deficit_blk(nSym-inLen/2)=" << deficitBlk
                  << " tStart=" << tLoopStart << " tEnd=" << tLoopEnd
                  << " dT=" << (tLoopEnd - tLoopStart) << " omega_=" << omega_
                  << " cursorClamped=" << (cursorClamped ? 1 : 0)
                  << " cumErr(sumOut-sumIn/2)=" << globalErr << std::endl;
    }
#endif

    // Keep absWrite_ and the interpolation time base tAbs_ in a moderate range so that
    // (1) tAbs_ keeps enough fractional precision in double for Lagrange mu, and
    // (2) indices passed to GetRingI stay safely within int range (it uses int masking).
    // Subtracting one full ring length from both preserves (absIdx mod kRingSize) and
    // all CanInterp inequalities; use a loop in case a single block pushes inLen large
    // enough to cross more than one wrap boundary (unlikely with current callers).
    while (absWrite_ >= 2LL * kRingSize)
    {
        absWrite_ -= static_cast<long long>(kRingSize);
        tAbs_ -= static_cast<double>(kRingSize);
    }

    return outCount;
}

void GardnerTiming::PushOmegaHistoryOnUpdate()
{
    if (historyOmega_.empty())
        historyOmega_.assign(static_cast<size_t>(kOmegaLockMeasures), omegaNom_);

    omegaLockStrideCounter_++;
    if (omegaLockStrideCounter_ < omegaLockStrideUpdates_)
        return;
    omegaLockStrideCounter_ = 0;

    historyOmega_[static_cast<size_t>(historyOmegaWr_)] = omega_;
    historyOmegaWr_ = (historyOmegaWr_ + 1) % kOmegaLockMeasures;
    historyOmegaFill_ = std::min(historyOmegaFill_ + 1, kOmegaLockMeasures);

    if (historyOmegaFill_ < kOmegaLockMeasures)
    {
        locked_ = false;
        return;
    }

    double vmin = historyOmega_[0];
    double vmax = historyOmega_[0];
    for (int i = 1; i < kOmegaLockMeasures; ++i)
    {
        vmin = std::min(vmin, historyOmega_[static_cast<size_t>(i)]);
        vmax = std::max(vmax, historyOmega_[static_cast<size_t>(i)]);
    }
    const bool prevLocked = locked_;
    lastOmegaLockSpan_ = std::abs(vmax - vmin);
    locked_ = (lastOmegaLockSpan_ <= kOmegaLockSpanThreshold);
    // State-change logging is handled at the receiver level so we can print both:
    // - wall-clock simulation time (t_sim)
    // - rate-based time from the sample counters (t_rate)
}

void GardnerTiming::PushToRing(const float* inI, const float* inQ, int inLen)
{
    for (int i = 0; i < inLen; ++i)
    {
        const long long a = absWrite_ + static_cast<long long>(i);
        ringI_[static_cast<size_t>(static_cast<int>(a) & kRingMask)] = inI[i];
        ringQ_[static_cast<size_t>(static_cast<int>(a) & kRingMask)] = inQ[i];
    }
    absWrite_ += static_cast<long long>(inLen);
}

long long GardnerTiming::OldestAbs() const
{
    return absWrite_ - static_cast<long long>(kRingSize);
}

double GardnerTiming::NewestSafe() const
{
    return static_cast<double>(absWrite_ - 4);
}

bool GardnerTiming::CanInterp(double t) const
{
    const long long k = static_cast<long long>(std::floor(t));
    return (k - 1 >= OldestAbs()) && (k + 2 < absWrite_);
}

float GardnerTiming::GetRingI(long long absIdx) const
{
    return ringI_[static_cast<size_t>(static_cast<int>(absIdx) & kRingMask)];
}

float GardnerTiming::GetRingQ(long long absIdx) const
{
    return ringQ_[static_cast<size_t>(static_cast<int>(absIdx) & kRingMask)];
}

float GardnerTiming::InterpLagrange4I(double t) const
{
    const long long k = static_cast<long long>(std::floor(t));
    const double mu = t - static_cast<double>(k);
    return Lagrange4Simd::EvalSet_ps(GetRingI(k - 1), GetRingI(k + 0), GetRingI(k + 1), GetRingI(k + 2),
                                     Lagrange4Simd::Coeffs_ps(mu));
}

float GardnerTiming::InterpLagrange4Q(double t) const
{
    const long long k = static_cast<long long>(std::floor(t));
    const double mu = t - static_cast<double>(k);
    return Lagrange4Simd::EvalSet_ps(GetRingQ(k - 1), GetRingQ(k + 0), GetRingQ(k + 1), GetRingQ(k + 2),
                                     Lagrange4Simd::Coeffs_ps(mu));
}

void GardnerTiming::InterpLagrange4IQ(double t, float* outI, float* outQ) const
{
    const long long k = static_cast<long long>(std::floor(t));
    const double mu = t - static_cast<double>(k);
    Lagrange4Simd::EvalIQ_ps(GetRingI(k - 1), GetRingI(k + 0), GetRingI(k + 1), GetRingI(k + 2),
                             GetRingQ(k - 1), GetRingQ(k + 0), GetRingQ(k + 1), GetRingQ(k + 2), mu, outI, outQ);
}

#ifdef DEBUG_GARDNER_OUTPUTS
void GardnerTiming::OpenDebugFilesIfNeeded()
{
    if (fidEarly_ && fidOnTime_ && fidLate_ && fidErr_ && fidOmega_ && fidPhi_ && fidOut8Sps_)
        return;
    mkdir("../data", 0755);
    if (!fidEarly_)
        fidEarly_ = std::fopen("../data/gardner_early.bin", "wb");
    if (!fidOnTime_)
        fidOnTime_ = std::fopen("../data/gardner_ontime.bin", "wb");
    if (!fidLate_)
        fidLate_ = std::fopen("../data/gardner_late.bin", "wb");
    if (!fidErr_)
        fidErr_ = std::fopen("../data/gardner_err.bin", "wb");
    if (!fidOmega_)
        fidOmega_ = std::fopen("../data/gardner_omega.bin", "wb");
    if (!fidPhi_)
        fidPhi_ = std::fopen("../data/gardner_phi.bin", "wb");
    if (!fidOut8Sps_)
        fidOut8Sps_ = std::fopen("../data/gardner_output_8sps_iq.bin", "wb");
}

void GardnerTiming::CloseDebugFiles()
{
    if (fidEarly_)
    {
        std::fclose(static_cast<FILE*>(fidEarly_));
        fidEarly_ = nullptr;
    }
    if (fidOnTime_)
    {
        std::fclose(static_cast<FILE*>(fidOnTime_));
        fidOnTime_ = nullptr;
    }
    if (fidLate_)
    {
        std::fclose(static_cast<FILE*>(fidLate_));
        fidLate_ = nullptr;
    }
    if (fidErr_)
    {
        std::fclose(static_cast<FILE*>(fidErr_));
        fidErr_ = nullptr;
    }
    if (fidOmega_)
    {
        std::fclose(static_cast<FILE*>(fidOmega_));
        fidOmega_ = nullptr;
    }
    if (fidPhi_)
    {
        std::fclose(static_cast<FILE*>(fidPhi_));
        fidPhi_ = nullptr;
    }
    if (fidOut8Sps_)
    {
        std::fclose(static_cast<FILE*>(fidOut8Sps_));
        fidOut8Sps_ = nullptr;
    }
}
#endif

