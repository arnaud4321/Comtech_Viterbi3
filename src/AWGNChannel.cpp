#include "AWGNChannel.h"
#include "Transmitter.h"
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
void AWGNChannel::StartThreads(float EsN0db, double TimeDrift, double FrequencyShift, bool RandomFrequency)
{
    StopAll = false;
    double TxPower = pTx->GetTxPower();
    double kTx = sqrt(TxPower);
    double Pow0 = 2.0;//constant
    double N0 = Pow0*pow(10.0,-0.1*EsN0db);
    double TargetPower = pow(2.0,2*NumBits-1)*pow(10.0,-0.1*Backoff);
    double TotalPower = TxPower + N0;
    double kAll = sqrt(TargetPower/TotalPower);
    double kSig = kAll/kTx;
    Stdn = sqrt(N0*0.5)*kAll;
    mkSig = _mm256_set1_ps(kSig);
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
        while(!NoiseQ.AvailableWrite())
        {
            std::mutex mtx;
            std::unique_lock<std::mutex> lck(mtx);
            CvOutNoise.wait(lck);
        }

        
        if(StopAll)
            break;
        unsigned int PtrWr = NoiseQ.GetPtrWr();
        
        oNoiseGen.randn(Noise[PtrWr],NoiseBatchSize);

        NoiseQ.AdvanceWrite();
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
    condition_variable *pCvTxChn, *pCvChnTx;

    pTx->GetCVOut(pCvChnTx, pCvTxChn );

    int NumSamples = 0;

    while(!StopAll)
    {
        while((!NoiseQ.AvailableRead())&& (StopAll==false))
        {
            std::mutex mtx;
            std::unique_lock<std::mutex> lck(mtx);
            CvNoiseOut.wait(lck);
        }

        if(StopAll)
            break;
        unsigned int PtrRdNoise = NoiseQ.GetPtrRd();

        float *TxOut = pTx->GetOutput();

        while((TxOut == 0) && (StopAll == false))
        {
            std::mutex mtx;
            std::unique_lock<std::mutex> lck(mtx);
            pCvTxChn->wait(lck);
            TxOut = pTx->GetOutput();
        }
        if(StopAll)
            break;

        while(OutputBuffer.AlmostFull() && (StopAll == false))//If buffer is almost full - wait
        {
            std::mutex mtx;
            std::unique_lock<std::mutex> lck(mtx);
            CvUserOut.wait(lck);
        }
        if(StopAll)
            break;
        short *Output = OutputBuffer.GetWriteBuffer(TxOutputBatchSize);
        
        unsigned int PtrOut = 0;
        __m256 mStdn = _mm256_set1_ps(Stdn);
        for(int i = 0; i < TxOutputBatchSize; )
        {
            __m256i mInl = _mm256_loadu_si256((__m256i*)(TxOut+i));
            __m256 mOutl = _mm256_loadu_ps(Noise[PtrRdNoise]+i);
            i+=8;
            __m256i mInh = _mm256_loadu_si256((__m256i*)(TxOut+i));
            __m256 mOuth = _mm256_loadu_ps(Noise[PtrRdNoise]+i);
            i+= 8;
            mOutl = _mm256_mul_ps(mStdn,mOutl);
            mOuth = _mm256_mul_ps(mStdn,mOuth);
            
            mOutl = _mm256_fmadd_ps(mInl, mkSig , mOutl);
            mOuth = _mm256_fmadd_ps(mInh, mkSig , mOuth);
            
            __m256i mOuti =  floats_to_shorts_sat_perm_avx2(mOutl,mOuth);
            _mm256_storeu_si256((__m256i*)(Output+PtrOut),mOuti);
            PtrOut += 16;
        }
        NumSamples += PtrOut;
        //        cout<<"AWGN "<<NumSamples<<" "<<Output[0]<<" "<<Output[1]<<endl;

 # ifdef DEBUG_AWGN
        std::copy(Output,Output+TxOutputBatchSize,OutAllI+PtrOutAll);
        PtrOutAll += TxOutputBatchSize;
        cout<<"Collected "<<PtrOutAll<<endl;
        if(PtrOutAll >= 1100000)
        {    
            FILE *fid = fopen("AWGNOut.bin","wb");
            fwrite(OutAllI,sizeof(short),PtrOutAll,fid);
            fclose(fid);
            cout<<"Saved AWGN "<<PtrOutAll<<endl;
            std::this_thread::sleep_for(10ms);

            Finish = true;
        }
        #endif
        pTx->AdvanceOut();
        pCvChnTx->notify_one();
        NoiseQ.AdvanceRead();
        CvOutNoise.notify_one();
        OutputBuffer.AdvancePtrWr(TxOutputBatchSize);
        CvOutUser.notify_one();
    }
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

 short * AWGNChannel::GetOutput(int Size)
 {
    return OutputBuffer.GetReadBuffer(Size);
 }
void AWGNChannel::AdvanceOut(int Size)
{
    OutputBuffer.AdvancePtrRd(Size);
}