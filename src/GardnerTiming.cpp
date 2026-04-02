/**
 * @file GardnerTiming.cpp
 * @brief Oversampled ring, Lagrange-4 interpolation, Gardner TED, PI loop on @f$\omega@f$, lock from @f$\omega@f$ stability.
 *
 * @details **ProcessBlock** — @c PushToRing then loop: one TED term and one on-time output per symbol.
 * @f$\omega@f$ / PI run only every @c updatePeriod_ **output symbols** (~ @f$\mathrm{updatePeriod\_}\,\omega/F_s@f$ s apart).
 * See @ref GardnerTiming.
 */

#include "GardnerTiming.h"
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

    // Strict no-clamp: if cursor exits valid interpolation region, re-arm and output 0.
    if ((tAbs_ < oldestSafe) || (tAbs_ > newestSafeNow))
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
        // Need availability for t, t-0.5 and t+0.5.
        if (!CanInterp(tAbs_) || !CanInterp(tAbs_ - 0.5) || !CanInterp(tAbs_ + 0.5))
            break;

        const float yI = InterpLagrange4I(tAbs_);
        const float yQ = InterpLagrange4Q(tAbs_);
        // outI[outCount] = yI;
        // outQ[outCount] = yQ;

        const float eI = InterpLagrange4I(tAbs_ - 0.5);
        const float eQ = InterpLagrange4Q(tAbs_ - 0.5);
        const float lI = InterpLagrange4I(tAbs_ + 0.5);
        const float lQ = InterpLagrange4Q(tAbs_ + 0.5);
        
        outI[outCount] = InterpLagrange4I(tAbs_);
        outQ[outCount] = InterpLagrange4Q(tAbs_);

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
                    const float ii = InterpLagrange4I(t8[k]);
                    const float qq = InterpLagrange4Q(t8[k]);
                    float iq[2] = {ii, qq};
                    std::fwrite(iq, sizeof(float), 2, static_cast<FILE*>(fidOut8Sps_));
                }
            }
        }
#endif

        // Gardner TED on complex signal (use current on-time symbol).
        errorAcc_ += (lI - eI) * yI + (lQ - eQ) * yQ;
 //      errorAcc_ += (eI - lI) * yI + (eQ - lQ) * yQ;

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
            errorAcc_=0.0;
            PushOmegaHistoryOnUpdate();
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
    const float x0 = GetRingI(k - 1);
    const float x1 = GetRingI(k + 0);
    const float x2 = GetRingI(k + 1);
    const float x3 = GetRingI(k + 2);
    const double c0 = -mu * (mu - 1.0) * (mu - 2.0) / 6.0;
    const double c1 = (mu + 1.0) * (mu - 1.0) * (mu - 2.0) / 2.0;
    const double c2 = -(mu + 1.0) * mu * (mu - 2.0) / 2.0;
    const double c3 = (mu + 1.0) * mu * (mu - 1.0) / 6.0;
    return static_cast<float>(c0 * x0 + c1 * x1 + c2 * x2 + c3 * x3);
}

float GardnerTiming::InterpLagrange4Q(double t) const
{
    const long long k = static_cast<long long>(std::floor(t));
    const double mu = t - static_cast<double>(k);
    const float x0 = GetRingQ(k - 1);
    const float x1 = GetRingQ(k + 0);
    const float x2 = GetRingQ(k + 1);
    const float x3 = GetRingQ(k + 2);
    const double c0 = -mu * (mu - 1.0) * (mu - 2.0) / 6.0;
    const double c1 = (mu + 1.0) * (mu - 1.0) * (mu - 2.0) / 2.0;
    const double c2 = -(mu + 1.0) * mu * (mu - 2.0) / 2.0;
    const double c3 = (mu + 1.0) * mu * (mu - 1.0) / 6.0;
    return static_cast<float>(c0 * x0 + c1 * x1 + c2 * x2 + c3 * x3); 
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

