/**
 * @file Sampler.h
 * @brief Real-time-paced thread: pulls IQ shorts from @ref AWGNChannel into a ring for @ref Receiver.
 */
#pragma once
#include <thread>
#include <condition_variable>
#include <mutex>
class AWGNChannel;
#include "BufferShort.h"
#include <atomic>
using namespace std;
//#define DEBUG_SAMPLER

/**
 * @brief Producer thread converting channel output to a @c BufferShort consumed by the RX filter.
 *
 * @details Pulls interleaved IQ shorts from @ref AWGNChannel::CopyOutputSamples, copies into @c oBuffer under
 * @c mtxSamplerBuffer_ (waits if @c AlmostFull), and notifies @ref ReadFilterBatch consumers.
 *
 * **Pacing (always):** after each committed batch, if wall time is **ahead** of
 * @f$\texttt{NumBatches}\cdot\texttt{TimeBatch}@f$, the thread sleeps for the remainder.
 * @c TimeBatch = @c SPB / @c SamplingFrequency at construction, multiplied by **10** when @c Debug is true
 * (constructor flag), so debug runs slower in wall time with the same sample schedule.
 *
 * **Telemetry:** @c lastThroughputMsps_ is averaged over @c throughputMeasurePeriodSec_; optional @c [Sampler]
 * console lines use @c displayPeriodSec_.
 */
class Sampler
{
private:
    /* data */
    
    #ifdef DEBUG_SAMPLER
    short *OutAllI = 0;
    unsigned int PtrOutAll = 0;
    #define STOP_SAMPLE_S 6000000
    #endif
    
    bool Debug;
    std::atomic<bool> StopAll{false};
    double TimeBatch;
    double displayPeriodSec_ = 1.0;
    /// Wall-time window for lastThroughputMsps_ (constellation overlay); console line uses displayPeriodSec_ when >0.
    double throughputMeasurePeriodSec_ = 1.0;
    thread SamplingThread;
    AWGNChannel *pChannel;
    BufferShort oBuffer;
    std::mutex mtxSamplerBuffer_;
    condition_variable CvSamplerUser;
    void OperateSampler(void);
    std::atomic<double> lastThroughputMsps_{0.0};
public:
    uint64_t NumBatches;
    Sampler(bool DebugIn);
    ~Sampler();

    /** @brief Console status period; \<=0 disables periodic @c [Sampler] throughput lines. */
    void SetDisplayPeriodSec(double sec) { displayPeriodSec_ = sec; }

    /** @brief Wall-time window for @ref GetLastThroughputMsps (independent of console when decoupled). */
    void SetThroughputMeasurePeriodSec(double sec) { throughputMeasurePeriodSec_ = (sec > 1e-9) ? sec : 1.0; }

    void StartThread(void);
    void StopThread(void);
    void SetChannel(AWGNChannel *p)
    {
        pChannel = p;
    }
    /**
     * @brief Copy one batch from the ring and advance read pointer (mutex vs producer).
     * @param dst Destination interleaved I,Q shorts.
     * @param nShorts Number of shorts to read (typically @c 2*SPB).
     * @param stopAll Global stop flag from caller context.
     */
    bool ReadFilterBatch(short* dst, int nShorts, std::atomic<bool>& stopAll);

    /** @brief Wake threads waiting on ring space / data. */
    void NotifyFilterWaiters() { CvSamplerUser.notify_all(); }

    condition_variable* GetCvOut(void)
    {
        return &CvSamplerUser;
    }

    /** @brief Last measured complex-Msps over the window set by @ref SetThroughputMeasurePeriodSec. */
    double GetLastThroughputMsps() const { return lastThroughputMsps_.load(std::memory_order_relaxed); }
};

