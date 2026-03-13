#include "Sampler.h"
#include "AWGNChannel.h"
#include <thread>
#include <iostream>
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
    CvSamplerUser.notify_one();
    StopAll = true;
    if(SamplingThread.joinable())
        SamplingThread.join();
}
void Sampler::OperateSampler(void)
{
    NumBatches = 0;
    condition_variable *pCvFromCh, *pCvToCh;
    pChannel->GetCvs2User(pCvToCh, pCvFromCh);//the function is defined from the point of view of Channel
    auto Start = std::chrono::high_resolution_clock::now();
    auto LastDisplayTime = Start;
    while(!StopAll)
    {
        short *ChOut = pChannel->GetOutput(SPB*2);
        while((ChOut == 0) && (StopAll == false))
        {
            std::mutex mtx;
            std::unique_lock<std::mutex> lck(mtx);
            pCvFromCh->wait(lck);
            ChOut = pChannel->GetOutput(SPB*2);
        }
        //cout<<"PtrRd "<<pChannel->OutputBuffer.GetRdPtr()<<endl;
        //cout<<ChOut[0]<<" "<<ChOut[1]<<endl;
        if(StopAll)
            break;
        short *BufOut = oBuffer.GetWriteBuffer(SPB*2);
        std::copy(ChOut,ChOut+2*SPB,BufOut);
        pChannel->AdvanceOut(SPB*2);
        pCvToCh->notify_one();//Signal the channel to try and advance

        oBuffer.AdvancePtrWr(SPB*2);
        
        # ifdef DEBUG_SAMPLER
        std::copy(BufOut,BufOut+SPB*2,OutAllI+PtrOutAll);
        PtrOutAll += SPB*2;
        if(PtrOutAll >= STOP_SAMPLE_S)
        {    
            FILE *fid = fopen("SamplerOut.bin","wb");
            fwrite(OutAllI,sizeof(short),PtrOutAll,fid);
            fclose(fid);
            Finish = true;
            while(!Finish)
                std::this_thread::sleep_for(1ms);
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
        auto TimeEnd = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> sinceLastDisplay = TimeEnd - LastDisplayTime;
        if (sinceLastDisplay.count() >= 1.0)
        {
            std::chrono::duration<double> elapsed = TimeEnd - Start;
            double totalSamples = NumBatches * static_cast<double>(SPB);
            double avgMsps = (elapsed.count() > 0) ? (totalSamples / elapsed.count() / 1e6) : 0;
            std::cout << "Sampler: average " << avgMsps << " Msps" << std::endl;
            LastDisplayTime = TimeEnd;
        }
    }
}

