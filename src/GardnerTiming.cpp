#include "GardnerTiming.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
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
    tCursor_ = 0.0; // grid 0,2,4,... aligned with TakeEven (2x -> 1x)
    hasOverlap_ = false;
    for (int i = 0; i < kOverlap; ++i)
    {
        overlapI_[i] = 0.0f;
        overlapQ_[i] = 0.0f;
    }
#ifdef DEBUG_GARDNER_DEBIT
    gGardnerDebit = GardnerDebitStats{};
#endif
}

int GardnerTiming::ProcessBlock(const float* inI, const float* inQ, int inLen,
                                float* outI, float* outQ, int outMax)
{
    if (!inI || !inQ || !outI || !outQ || inLen < 8 || outMax <= 0)
        return 0;

    // Build a continuous input by prepending the overlap from previous block.
    const int prefix = hasOverlap_ ? kOverlap : 0;
    std::vector<float> workI(static_cast<size_t>(prefix + inLen));
    std::vector<float> workQ(static_cast<size_t>(prefix + inLen));
    if (prefix > 0)
    {
        for (int i = 0; i < kOverlap; ++i)
        {
            workI[i] = overlapI_[i];
            workQ[i] = overlapQ_[i];
        }
    }
    for (int i = 0; i < inLen; ++i)
    {
        workI[prefix + i] = inI[i];
        workQ[prefix + i] = inQ[i];
    }

    const int workLen = prefix + inLen;
    double t = tCursor_ + static_cast<double>(prefix);
    if (t < 0.0)
        t = 0.0;

    const double tLoopStart = t;
    int outCount = 0;

#ifdef DEBUG_GARDNER_OUTPUTS
    OpenDebugFilesIfNeeded();
    std::vector<float> dbgEarly;
    std::vector<float> dbgOnTime;
    std::vector<float> dbgLate;
    std::vector<float> dbgErr;
    std::vector<float> dbgOmega;
    dbgEarly.reserve(static_cast<size_t>(outMax) * 2u);
    dbgOnTime.reserve(static_cast<size_t>(outMax) * 2u);
    dbgLate.reserve(static_cast<size_t>(outMax) * 2u);
    dbgErr.reserve(static_cast<size_t>(outMax));
    dbgOmega.reserve(static_cast<size_t>(outMax));
#endif

    while (outCount < outMax)
    {
        // Last usable even index in work: workLen-2 (inLen/2 symbols per block).
        // With mu recomputed after clamp in InterpLagrange4, t±0.5 stays inside the stencil.
        if (t > static_cast<double>(workLen - 2))
            break;

        const float yI = InterpLagrange4(workI.data(), workLen, t);
        const float yQ = InterpLagrange4(workQ.data(), workLen, t);
        outI[outCount] = yI;
        outQ[outCount] = yQ;

        const float eI = InterpLagrange4(workI.data(), workLen, t - 0.5);
        const float eQ = InterpLagrange4(workQ.data(), workLen, t - 0.5);
        const float lI = InterpLagrange4(workI.data(), workLen, t + 0.5);
        const float lQ = InterpLagrange4(workQ.data(), workLen, t + 0.5);

#ifdef DEBUG_GARDNER_OUTPUTS
        dbgEarly.push_back(eI);
        dbgEarly.push_back(eQ);
        dbgOnTime.push_back(yI);
        dbgOnTime.push_back(yQ);
        dbgLate.push_back(lI);
        dbgLate.push_back(lQ);
#endif

        // Gardner TED on complex signal (use current on-time symbol).
        const float err = (lI - eI) * yI + (lQ - eQ) * yQ;
//        const float err = (eI - lI) * yI + (eQ - lQ) * yQ;

#ifdef DEBUG_GARDNER_OUTPUTS
        dbgErr.push_back(err);
        dbgOmega.push_back(static_cast<float>(omega_));
#endif

        updateCounter_++;
        if (updateCounter_ >= updatePeriod_)
        {
            updateCounter_ = 0;
            integ_ += ki_ * static_cast<double>(err);
            omega_ = omegaNom_ + integ_ + kp_ * static_cast<double>(err);
            omega_ = std::max(1.7, std::min(2.3, omega_));
        }

        t += omega_;
        outCount++;
    }

#ifdef DEBUG_GARDNER_OUTPUTS
    if (fidEarly_ && fidOnTime_ && fidLate_ && fidErr_ && fidOmega_ && !dbgEarly.empty())
    {
        const size_t nfloat = dbgEarly.size();
        std::fwrite(dbgEarly.data(), sizeof(float), nfloat, static_cast<FILE*>(fidEarly_));
        std::fwrite(dbgOnTime.data(), sizeof(float), nfloat, static_cast<FILE*>(fidOnTime_));
        std::fwrite(dbgLate.data(), sizeof(float), nfloat, static_cast<FILE*>(fidLate_));
        std::fwrite(dbgErr.data(), sizeof(float), static_cast<size_t>(outCount), static_cast<FILE*>(fidErr_));
        std::fwrite(dbgOmega.data(), sizeof(float), static_cast<size_t>(outCount), static_cast<FILE*>(fidOmega_));
    }
#endif

    // Carry continuous timing cursor to next block (in "new block" coordinates).
    const double tLoopEnd = t;
    tCursor_ = tLoopEnd - static_cast<double>(workLen);
    if (tCursor_ < -2.0)
        tCursor_ = -2.0;
    if (tCursor_ > static_cast<double>(kOverlap))
        tCursor_ = static_cast<double>(kOverlap);

    // Update overlap with the last kOverlap input samples for next call.
    for (int i = 0; i < kOverlap; ++i)
    {
        overlapI_[i] = inI[inLen - kOverlap + i];
        overlapQ_[i] = inQ[inLen - kOverlap + i];
    }
    hasOverlap_ = true;

#ifdef DEBUG_GARDNER_DEBIT
    const bool cursorClamped = (tCursor_ <= -2.0) || (tCursor_ >= static_cast<double>(kOverlap));
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
        std::cerr << "[Gardner DEBIT] blk=" << gGardnerDebit.blocks << " prefix=" << prefix
                  << " workLen=" << workLen << " inLen=" << inLen << " nSym=" << outCount
                  << " ratio(2*nSym/inLen)=" << ratio << " deficit_blk(nSym-inLen/2)=" << deficitBlk
                  << " tStart=" << tLoopStart << " tEnd=" << tLoopEnd
                  << " dT=" << (tLoopEnd - tLoopStart) << " omega_=" << omega_
                  << " cursorClamped=" << (cursorClamped ? 1 : 0)
                  << " cumErr(sumOut-sumIn/2)=" << globalErr << std::endl;
    }
#endif
    return outCount;
}

#ifdef DEBUG_GARDNER_OUTPUTS
void GardnerTiming::OpenDebugFilesIfNeeded()
{
    if (fidEarly_ && fidOnTime_ && fidLate_ && fidErr_ && fidOmega_)
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
}
#endif

inline int GardnerTiming::ClampIndex(int idx, int lo, int hi)
{
    return std::max(lo, std::min(hi, idx));
}

inline float GardnerTiming::InterpLagrange4(const float* x, int len, double t)
{
    int k = static_cast<int>(std::floor(t));
    k = ClampIndex(k, 1, len - 3);
    const double mu = t - static_cast<double>(k);
    const float x0 = x[k - 1];
    const float x1 = x[k + 0];
    const float x2 = x[k + 1];
    const float x3 = x[k + 2];

    const double c0 = -mu * (mu - 1.0) * (mu - 2.0) / 6.0;
    const double c1 = (mu + 1.0) * (mu - 1.0) * (mu - 2.0) / 2.0;
    const double c2 = -(mu + 1.0) * mu * (mu - 2.0) / 2.0;
    const double c3 = (mu + 1.0) * mu * (mu - 1.0) / 6.0;

    return static_cast<float>(c0 * x0 + c1 * x1 + c2 * x2 + c3 * x3);
}

