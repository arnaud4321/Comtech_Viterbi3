#include "Sampler.h"
#include "AWGNChannel.h"
#include <chrono>
#include <cstring>
#include <thread>



extern bool Finish;
Sampler::Sampler(bool DebugIn):Debug(DebugIn),oBuffer(128*SPB,2*SPB)
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

void Sampler::OperateSampler(void)
{
    NumBatches = 0;
    auto Start = std::chrono::high_resolution_clock::now();

    auto ThroughputWindowStart = std::chrono::steady_clock::now();
    uint64_t SamplesInWindow = 0;

    while(!StopAll)
    {
        alignas(32) short chTmp[SPB * 2];
        if (!pChannel->CopyOutputSamples(chTmp, SPB * 2, StopAll))
            break;
        {
            std::unique_lock<std::mutex> lk(mtxSamplerBuffer_);
            // Backpressure: never let the producer advance
            // when the ring is almost full.
            while (!StopAll && oBuffer.AlmostFull())
                CvSamplerUser.wait(lk);

            short* BufOut = oBuffer.GetWriteBuffer(SPB * 2);
            std::copy(chTmp, chTmp + 2 * SPB, BufOut);
            oBuffer.AdvancePtrWr(SPB * 2);
        }
        SamplesInWindow += SPB; // complex samples per batch

        auto nowTp = std::chrono::steady_clock::now();
        double dt = std::chrono::duration<double>(nowTp - ThroughputWindowStart).count();
        if (dt >= 1.0)
        {
            double msps = static_cast<double>(SamplesInWindow) / dt / 1e6;
            std::cout << "Sampler throughput: " << msps << " Msps" << std::endl;

            ThroughputWindowStart = nowTp;
            SamplesInWindow = 0;
        }

        # ifdef DEBUG_SAMPLER
        std::copy(chTmp, chTmp + SPB * 2, OutAllI + PtrOutAll);
        PtrOutAll += SPB*2;
        if(PtrOutAll >= STOP_SAMPLE_S)
        {    
            FILE *fid = fopen("SamplerOut.bin","wb");
            fwrite(OutAllI,sizeof(short),PtrOutAll,fid);
            fclose(fid);
            Finish = true;
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

