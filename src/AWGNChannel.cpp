#include "AWGNChannel.h"
#include "Transmitter.h"
#include <cstring>
extern bool Finish;
extern mutex mtxfilethr;
AWGNChannel::AWGNChannel(unsigned int Seed):NoiseQ(LengthQueue),OutputBuffer(128*TxOutputBatchSize,2*TxOutputBatchSize)
{
    oNoiseGen.set_seed(Seed);
    NoiseBatchSize = (int) ((double) TxOutputBatchSize * 1.001); //make sure there will be enough noise even if time is slowed down
    int Tmp = (NoiseBatchSize >> 4)<<4; //make a length of 16
    if(NoiseBatchSize > Tmp)
        Tmp += 16; 
    NoiseBatchSize = Tmp;
    Noise[0] = (float *) _mm_malloc(sizeof(float) * NoiseBatchSize * LengthQueue,32);
    for(int i = 1; i < LengthQueue; i++)
    {
        Noise[i] = Noise[i-1] + TxOutputBatchSize;
    }
    #ifdef DEBUG_AWGN
        OutAllI = (short *) _mm_malloc( 20000000*sizeof(short),32);
    #endif
}

AWGNChannel::~AWGNChannel()
{
    _mm_free(Noise[0]);
    #ifdef DEBUG_AWGN
        _mm_free(OutAllI);
    #endif
}
void AWGNChannel::StartThreads(void)
{
    StopAll = false;
    double TxPower = pTx->GetTxPower();
    double kTx = sqrt(TxPower);
    double Pow0 = 2.0;//constant
    EsN0db = pParams->EsN0;
    double N0 = Pow0*pow(10.0,-0.1*EsN0db);
    double TargetPower = pow(2.0,2*NumBits-1)*pow(10.0,-0.1*Backoff);
    double TotalPower = TxPower + N0;
    double kAll = sqrt(TargetPower/TotalPower);
    double kSig = kAll/kTx;
    Stdn = sqrt(N0*0.5)*kAll;
    mkSig = _mm256_set1_ps(kSig);
    double n2 = round(pParams->StablePeriod/TxOutputBatchDuration);
    double n1 = round(pParams->AccelerationPeriod/TxOutputBatchDuration);

    AccelerationPeriod = n1 *  TxOutputBatchDuration;
    StablePeriod = n2 * TxOutputBatchDuration;
    TotalPeriod = 2*(AccelerationPeriod+StablePeriod);
    TransitionCounter[0] = 0; 
    TransitionCounter[1] = n2;//End of Stable 1
    TransitionCounter[2] = TransitionCounter[1] + n1;//End of Acc 1
    TransitionCounter[3] = TransitionCounter[2] + n2;//End of Stable 2
    TransitionCounter[4] = TransitionCounter[3] + n1;//End of Acc 2
    
    
    
    NoiseQ.Reset();
    NoiseThread = std::thread(&AWGNChannel::GenerateNoise, this);
    OutputThread = std::thread(&AWGNChannel::GenerateOutput, this);

}



void AWGNChannel::StopThreads(void)
{
    StopAll = true;
    CvOutNoise.notify_one();
    CvOutUser.notify_one();
    CvNoiseOut.notify_one();
    CvUserOut.notify_one();
    if(NoiseThread.joinable())
        NoiseThread.join();
    if(OutputThread.joinable())
        OutputThread.join();
}


void AWGNChannel::GenerateNoise(void)
{


    #ifdef WRITE_LOG_THR
	mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"Generate Noise  Thread %d\n", gettid());
    fclose(fidthr);
    mtxfilethr.unlock();
    #endif
    while(!StopAll)
    {
        {
            std::unique_lock<std::mutex> lk(mtxNoiseQ_);
            CvOutNoise.wait(lk, [&] { return StopAll || NoiseQ.AvailableWrite(); });
        }

        
        if(StopAll)
            break;
        std::unique_lock<std::mutex> lk(mtxNoiseQ_);
        unsigned int PtrWr = NoiseQ.GetPtrWr();
        oNoiseGen.randn(Noise[PtrWr],NoiseBatchSize);
        NoiseQ.AdvanceWrite();
        lk.unlock();
        CvNoiseOut.notify_one();
    }
}


void AWGNChannel::GenerateOutput(void)
{
    #ifdef WRITE_LOG_THR
	mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"Channel Output Thread %d\n", gettid());
    fclose(fidthr);
    mtxfilethr.unlock();
    #endif
    int NumSamples = 0;
   
    while(!StopAll)
    {
        {
            std::unique_lock<std::mutex> lk(mtxNoiseQ_);
            CvNoiseOut.wait(lk, [&] { return StopAll || NoiseQ.AvailableRead(); });
        }

        if(StopAll)
            break;
        std::unique_lock<std::mutex> lkNoise(mtxNoiseQ_);
        unsigned int PtrRdNoise = NoiseQ.GetPtrRd();

        alignas(32) float txChunk[TxOutputBatchSize];
        if (!pTx->CopyOutputSamples(txChunk, TxOutputBatchSize, StopAll))
            break;

        {
            std::unique_lock<std::mutex> lk(mtxOutputBuffer_);
            while (OutputBuffer.AlmostFull() && !StopAll)
                CvUserOut.wait(lk); 
            if (StopAll)
                break;
            short* Output = OutputBuffer.GetWriteBuffer(TxOutputBatchSize);

            unsigned int PtrOut = 0;
            __m256 mStdn = _mm256_set1_ps(Stdn);
            for (int i = 0; i < TxOutputBatchSize;)
            {
                __m256 mInl = _mm256_loadu_ps(txChunk + i);
                __m256 mOutl = _mm256_loadu_ps(Noise[PtrRdNoise] + i);
                i += 8;
                __m256 mInh = _mm256_loadu_ps((txChunk + i));
                __m256 mOuth = _mm256_loadu_ps(Noise[PtrRdNoise] + i);
                i += 8;
                mOutl = _mm256_mul_ps(mStdn, mOutl);
                mOuth = _mm256_mul_ps(mStdn, mOuth);

                mOutl = _mm256_fmadd_ps(mInl, mkSig, mOutl);
                mOuth = _mm256_fmadd_ps(mInh, mkSig, mOuth);

                __m256i mOuti = floats_to_shorts_sat_perm_avx2(mOutl, mOuth);
                _mm256_storeu_si256((__m256i*)(Output + PtrOut), mOuti);
                PtrOut += 16;
            }
            NumSamples += PtrOut;

#ifdef DEBUG_AWGN
            std::copy(Output, Output + TxOutputBatchSize, OutAllI + PtrOutAll);
            PtrOutAll += TxOutputBatchSize;
            cout << "Collected " << PtrOutAll << endl;
            if (PtrOutAll >= 1100000)
            {
                FILE* fid = fopen("AWGNOut.bin", "wb");
                fwrite(OutAllI, sizeof(short), PtrOutAll, fid);
                fclose(fid);
                cout << "Saved AWGN " << PtrOutAll << endl;
                std::this_thread::sleep_for(10ms);

                Finish = true;
            }
#endif
            OutputBuffer.AdvancePtrWr(TxOutputBatchSize);
        }
        NoiseQ.AdvanceRead();
        lkNoise.unlock();
        CvOutNoise.notify_one();
        CvOutUser.notify_one();
    }
}

bool AWGNChannel::CopyOutputSamples(short* dst, int nShorts, bool& stopAll)
{
    std::unique_lock<std::mutex> lk(mtxOutputBuffer_);
    short* p = OutputBuffer.GetReadBuffer(nShorts);
    while (p == nullptr && !stopAll)
    {
        CvOutUser.wait_for(lk, std::chrono::milliseconds(1));
        p = OutputBuffer.GetReadBuffer(nShorts);
    }
    if (stopAll)
        return false;
    if (!p)
        return false;
    std::memcpy(dst, p, sizeof(short) * static_cast<size_t>(nShorts));
    OutputBuffer.AdvancePtrRd(nShorts);
    lk.unlock();
    CvUserOut.notify_one();
    return true;
}




// Convert two __m256 floats into 16 int16 with saturation
// Output layout: [x0..x7 y0..y7]
__m256i AWGNChannel::floats_to_shorts_sat_perm_avx2(__m256 x, __m256 y)
{
    // 1) float -> int32 (round to nearest)
    __m256i xi = _mm256_cvtps_epi32(x);
    __m256i yi = _mm256_cvtps_epi32(y);

    // 2) pack int32 -> int16 (lane-wise!)
    __m256i v = _mm256_packs_epi32(xi, yi);

    // 3) fix lane ordering
    // immediate = 0b11011000 = 0xD8
    v = _mm256_permute4x64_epi64(v, 0xD8);

    return v;
}
