#include "Receiver.h"
#include "Sampler.h"
#include "definitions.h"
#include <cstdio>
#include <immintrin.h>
#include <thread>
Receiver::Receiver(/* args */):oBufferFilter(SPB*256,32*SPB)
{
    FilterInI = (float*) _mm_malloc(SPB*2*sizeof(float),32);
    FilterInQ = FilterInI + SPB;
    OneSpsI = (float*) _mm_malloc(2*ReceiverInputBatchIQSymbols*sizeof(float),32);
    OneSpsQ = OneSpsI + ReceiverInputBatchIQSymbols;
    SplitI[0] = (float*) _mm_malloc(BatchSize3*sizeof(float),32);
    SplitQ[0] = (float*) _mm_malloc(BatchSize3*sizeof(float),32);
    Outputs[0] = (unsigned char*) _mm_malloc(BatchSize3,32);

    for(int i = 1; i < 3; i++)
    {
        SplitI[i] = SplitI[i-1] + BatchSize1;
        SplitQ[i] = SplitQ[i-1] + BatchSize1;
        Outputs[i] = Outputs[i-1] + BatchSize1;
    }


#ifdef DEBUG1
    OutAllI = (float*) _mm_malloc( (NUM_SAMPLES_R1+100000)*sizeof(float),32);
    OutAllQ = (float*) _mm_malloc( (NUM_SAMPLES_R1+100000)*sizeof(float),32);
#endif



}

Receiver::~Receiver()
{
    _mm_free(FilterInI);
    _mm_free(OneSpsI);
    _mm_free(SplitI[0]);
    _mm_free(SplitQ[0]);
    

    #ifdef DEBUG1
        _mm_free(OutAllI);
        _mm_free(OutAllQ);
    #endif

}
void Receiver::StartThreads(double RollOff)
{
    oBufferFilter.Reset();
    objRxFilter.CreateObjects(RollOff);
    FilterThread = std::thread(&Receiver::OperateFilter, this);
    ViterbiThread = std::thread(&Receiver::OperateViterbi, this);

}
void Receiver::StopThreads(void)
{
    CvFilterUser.notify_all();
    StopAll = true;
    if(FilterThread.joinable())
        FilterThread.join();
}

void Receiver::OperateFilter(void)
{

    condition_variable *pCvSamplerIn = pSampler->GetCvOut();
    while(!StopAll)
    {
        short *SamplerOut = pSampler->GetOutput(SPB*2);
        while((SamplerOut == 0) && (StopAll == false))
        {
            std::mutex mtx;
            std::unique_lock<std::mutex> lck(mtx);
            pCvSamplerIn->wait(lck);
            SamplerOut = pSampler->GetOutput(SPB*2);
        }
        if(StopAll)
            break;
        __m256 mouti;
        __m256 moutq;
        int PtrOut = 0;
        for(int i = 0; i < (2*SPB); i+= 16)
        {
            __m256i mIn = _mm256_loadu_si256((__m256i*) (SamplerOut+i));
           shorts_to_floats_avx2(mIn,&mouti,&moutq);
           _mm256_store_ps(FilterInI+PtrOut,mouti);
           _mm256_store_ps(FilterInQ+PtrOut,moutq);
           PtrOut += 8;
        }
        pSampler->AdvanceOut(2*SPB);
        float *FilterOutI, *FilterOutQ;
        oBufferFilter.GetWriteBuffer(FilterOutI,FilterOutQ,SPB);
        objRxFilter.CreateOutputs(FilterInI,FilterInQ,FilterOutI,FilterOutQ,SPB);
        oBufferFilter.AdvancePtrWr(SPB);
        CvFilterUser.notify_one();
    }
}


void Receiver::OperateViterbi(void)
{
    bool First = true;
    bool Calibrated = false;
    while(!StopAll)
    {
        //4 options for Viterbi - taking also care of Spectrum Inversion
        //(I+jQ)
        //I-jQ
        //j(I+jQ) = -Q+jI
        //j(I-jQ)= Q+jI
        //Controls:
        //  false 1 1
        //false 1 -1
        //true -1 1
        //true 1 1

        float *FilterOutI, *FilterOutQ;
        oBufferFilter.GetReadBuffer(FilterOutI,FilterOutQ,ReceiverInputBatchIQSamples);
        while((FilterOutI == 0) && (StopAll == false))
        {
            std::mutex mtx;
            std::unique_lock<std::mutex> lck(mtx);
            CvFilterUser.wait(lck);
            oBufferFilter.GetReadBuffer(FilterOutI,FilterOutQ,BatchSize3*2);
        }

        if(StopAll)
            break;
        #ifdef DEBUG1
        std::copy(FilterOutI,FilterOutI+BatchSize3*2,OutAllI+PtrOutAll);
        std::copy(FilterOutQ,FilterOutQ+BatchSize3*2,OutAllQ+PtrOutAll);
        PtrOutAll += BatchSize3*2;

        if(PtrOutAll >= NUM_SAMPLES_R1)
        {    
            cout<<"Collected "<<PtrOutAll<<endl;
            FILE *fid = fopen("FilterOut.bin","wb");
            fwrite(OutAllI,sizeof(float),PtrOutAll,fid);
            fwrite(OutAllQ,sizeof(float),PtrOutAll,fid);
            fclose(fid);
            exit(-1);
        }
        #endif


        if(First)
            First = false;
        else
        {
            TakeEvenDebug(FilterOutI, OneSpsI, ReceiverInputBatchIQSamples);
            TakeEvenDebug(FilterOutQ, OneSpsQ, ReceiverInputBatchIQSamples);
            Split3(OneSpsI, OneSpsQ, ReceiverInputBatchIQSymbols);
            float MetricsGrowth;
            oViterbi.Decode(SplitI[0], SplitQ[0], ReceiverInputBatchIQSymbols, Outputs[0], false, 1, 1,MetricsGrowth );

        }
        oBufferFilter.AdvancePtrRd(BatchSize3*2);

    }

}
void  Receiver::Split3(float *InputI, float *InputQ, int Length)
{
    __m256i mIndex[3];
    mIndex[0] = _mm256_setr_epi32(0, 3, 6, 9, 12, 15, 18, 21);
    mIndex[1] = _mm256_setr_epi32(1, 4, 7, 10, 13, 16, 19, 22);
    mIndex[2] = _mm256_setr_epi32(2, 5, 8, 11, 14, 17, 20, 23);
    //mIndex[3] = _mm256_set1_epi32(24);
    int PtrOut = 0;
    for(int i = 0; i<Length; i += 24)
    {
        __m256 I0 = _mm256_i32gather_ps(InputI + i,mIndex[0],4);
        __m256 Q0 = _mm256_i32gather_ps(InputQ + i,mIndex[0],4);
        __m256 I1 = _mm256_i32gather_ps(InputI + i,mIndex[1],4);
        __m256 Q1 = _mm256_i32gather_ps(InputQ + i,mIndex[1],4);
        __m256 I2 = _mm256_i32gather_ps(InputI + i,mIndex[2],4);
        __m256 Q2 = _mm256_i32gather_ps(InputQ + i,mIndex[2],4);
        _mm256_store_ps(SplitI[0]+PtrOut,I0);
        _mm256_store_ps(SplitQ[0]+PtrOut,Q0);
        _mm256_store_ps(SplitI[1]+PtrOut,I1);
        _mm256_store_ps(SplitQ[1]+PtrOut,Q1);
        _mm256_store_ps(SplitI[2]+PtrOut,I2);
        _mm256_store_ps(SplitQ[2]+PtrOut,Q2);
        PtrOut += 8;
    }


}
void Receiver::shorts_to_floats_avx2(__m256i v16, __m256* out0, __m256* out1)
{
    // Extract lower 128 bits (8 shorts)
   __m128i lo = _mm256_castsi256_si128(v16);

    // Extract upper 128 bits (8 shorts)
    __m128i hi = _mm256_extracti128_si256(v16, 1);

    // Sign-extend int16 -> int32
    __m256i lo32 = _mm256_cvtepi16_epi32(lo); //I0Q0I1Q1I2Q2I3Q3
    __m256i hi32 = _mm256_cvtepi16_epi32(hi);//I4Q4I5Q5I6Q6I7Q7
    lo32 = _mm256_shuffle_epi32(lo32,0xD8);//11011000 I0I1Q0Q1I2I3Q2Q3
    hi32 = _mm256_shuffle_epi32(hi32,0xD8);//I45,Q45,I67,Q67
     
    __m256i mi32 = _mm256_unpacklo_epi64 (lo32,hi32);//I01,I45,I23,I67 
    __m256i mq32 = _mm256_unpackhi_epi64 (lo32,hi32);//Q01,Q45,Q23,Q67
    mi32 = _mm256_permute4x64_epi64(mi32,0xD8);//11011000 
    mq32 = _mm256_permute4x64_epi64(mq32,0xD8);//11011000 

    // Convert int32 -> float
    *out0 = _mm256_cvtepi32_ps(mi32);
    *out1 = _mm256_cvtepi32_ps(mq32);
}

void Receiver::TakeEvenDebug(float *x, float *Output, int Length)
{
    //for(int i = 0; i<16;i++)
    //    x[i] = i;
    for (size_t i = 0; i < Length; i += 16) 
    {
    __m256 a0 = _mm256_loadu_ps(x + i);
    __m256 a1 = _mm256_loadu_ps(x + i + 8);

    __m256 e0 = _mm256_shuffle_ps(a0, a1, _MM_SHUFFLE(2,0,2,0));
    __m256d e0d = _mm256_castps_pd(e0);
    e0d = _mm256_permute4x64_pd(e0d, 0xD8);//1000
    e0 = _mm256_castpd_ps(e0d);


    _mm256_storeu_ps(Output + i/2, e0);
    }

}
