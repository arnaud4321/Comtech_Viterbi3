/**
 * @file Sampler.cpp
 * @brief Producer thread: channel shorts → @c BufferShort ring for the receiver matched-filter thread.
 *
 * @details **OperateSampler** loop: @c CopyOutputSamples from @ref AWGNChannel into a local batch, copy into
 * @c oBuffer under @c mtxSamplerBuffer_ (waits if @c AlmostFull — backpressure). Updates complex-Msps over
 * @c throughputMeasurePeriodSec_ and optional @c [Sampler] console line every @c displayPeriodSec_.
 * After each batch, if wall time lags nominal @c TimeBatch × @c NumBatches, sleeps to approximate real-time
 * pacing (@c Debug multiplies batch duration). Notifies @c CvSamplerUser so @c ReadFilterBatch can proceed.
 */

#include "Sampler.h"
#include "AWGNChannel.h"
#include "ConsoleAlert.h"
#include "definitions.h"
#include <chrono>
#include <iostream>
#include <cstring>
#include <thread>
#include <atomic>



extern std::atomic<bool> Finish;
Sampler::Sampler(bool DebugIn):Debug(DebugIn),oBuffer(kSamplerShortRingSize, kSamplerShortRingExtra)
{ 
#ifdef DEBUG_SAMPLER
    OutAllI = (short *) _mm_malloc( 20000000*sizeof(short),32);
#endif

    TimeBatch = double(SPB)/SamplingFrequency;
    if(Debug)
        TimeBatch *= 10;//debug mode, work slower
}

Sampler::~Sampler()
{
    #ifdef DEBUG_SAMPLER
    _mm_free(OutAllI);
    #endif
}

void Sampler::StartThread(void)
{
    oBuffer.Reset();
    SamplingThread = std::thread(&Sampler::OperateSampler, this);

}
void Sampler::StopThread(void)
{
    CvSamplerUser.notify_all();
    StopAll = true;
    if(SamplingThread.joinable())
        SamplingThread.join();
}
bool Sampler::ReadFilterBatch(short* dst, int nShorts, bool& stopAll)
{
    std::unique_lock<std::mutex> lk(mtxSamplerBuffer_);
    short* p = oBuffer.GetReadBuffer(nShorts);
    while (p == nullptr && !stopAll)
    {
        CvSamplerUser.wait_for(lk, std::chrono::duration<double>(1e-3));
        p = oBuffer.GetReadBuffer(nShorts);
    }
    if (stopAll)
        return false;
    if (!p)
        return false;
    std::memcpy(dst, p, sizeof(short) * static_cast<size_t>(nShorts));
    oBuffer.AdvancePtrRd(nShorts);
    // Consumer released space: wake the producer
    // (producer waits when AlmostFull()==true).
    lk.unlock();
    CvSamplerUser.notify_one();
    return true;
}

/**
 * @brief Thread body: pull IQ shorts from channel, push ring, throughput + optional pacing.
 */
void Sampler::OperateSampler(void)
{
    NumBatches = 0;
    auto Start = std::chrono::high_resolution_clock::now();

    auto ThroughputWindowStart = std::chrono::steady_clock::now();
    auto SamplerConsoleWindowStart = std::chrono::steady_clock::now();
    uint64_t SamplesInWindow = 0;
    auto lastSamplerRingSatLog = std::chrono::steady_clock::now();

    while(!StopAll)
    {
        alignas(32) short chTmp[SPB * 2];
        if (!pChannel->CopyOutputSamples(chTmp, SPB * 2, StopAll))
            break;
        {
            std::unique_lock<std::mutex> lk(mtxSamplerBuffer_);
            // Backpressure: never let the producer advance when the ring is almost full.
            while (!StopAll && oBuffer.AlmostFull())
            {
                const auto nowSat = std::chrono::steady_clock::now();
                if (nowSat - lastSamplerRingSatLog >= std::chrono::seconds(1))
                {
                    lastSamplerRingSatLog = nowSat;
                    CONSOLE_ALERT_STMT(std::cout << ConsoleAlert::kRedOpen
                                                 << "[Sampler] channel→RX short ring almost full; fill="
                                                 << oBuffer.GetSizeInBuffer() << "/" << oBuffer.GetBufferSize() - 1
                                                 << ConsoleAlert::kReset << std::endl;);
                }
                CvSamplerUser.wait(lk);
            }

            short* BufOut = oBuffer.GetWriteBuffer(SPB * 2);
            std::copy(chTmp, chTmp + 2 * SPB, BufOut);
            oBuffer.AdvancePtrWr(SPB * 2);
        }
        SamplesInWindow += SPB; // complex samples per batch

        auto nowTp = std::chrono::steady_clock::now();
        const double dtMeas =
            std::chrono::duration<double>(nowTp - ThroughputWindowStart).count();
        if (dtMeas >= throughputMeasurePeriodSec_)
        {
            const double msps = static_cast<double>(SamplesInWindow) / dtMeas / 1e6;
            lastThroughputMsps_.store(msps, std::memory_order_relaxed);
            ThroughputWindowStart = std::chrono::steady_clock::now();
            SamplesInWindow = 0;
        }
        const double dtConsole =
            std::chrono::duration<double>(nowTp - SamplerConsoleWindowStart).count();
        if (displayPeriodSec_ > 0.0 && dtConsole >= displayPeriodSec_)
        {
            std::cout << "[Sampler] throughput=" << lastThroughputMsps_.load(std::memory_order_relaxed)
                      << " Msps" << std::endl;
            SamplerConsoleWindowStart = std::chrono::steady_clock::now();
        }

        # ifdef DEBUG_SAMPLER
        std::copy(chTmp, chTmp + SPB * 2, OutAllI + PtrOutAll);
        PtrOutAll += SPB*2;
        if(PtrOutAll >= STOP_SAMPLE_S)
        {    
            FILE *fid = fopen("SamplerOut.bin","wb");
            fwrite(OutAllI,sizeof(short),PtrOutAll,fid);
            fclose(fid);
            Finish.store(true, std::memory_order_relaxed);
            cout<<"Exiting"<<endl;
            exit(-1);
        }
        #endif
        
        
        auto Now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = Now - Start;


        double ExpectedTime = NumBatches * TimeBatch;
        if (elapsed.count() < ExpectedTime) 
        {
            std::chrono::duration<double> remaining(ExpectedTime - elapsed.count());
            std::this_thread::sleep_for(remaining);
        }
        CvSamplerUser.notify_one();
        NumBatches++;
        //cout<<"Sampler Batch "<<NumBatches<<endl;
    }
}

