#include "Receiver.h"
extern std::mutex mtxfilethr;
Receiver::Receiver(/* args */):oBufferFilter(SPB*256,32*SPB),OutputQ(LengthQueue)
{
    FilterInI = (float*) _mm_malloc(SPB*2*sizeof(float),32);
    FilterInQ = FilterInI + SPB;
    OneSpsI = (float*) _mm_malloc(2*ReceiverInputBatchIQSymbols*sizeof(float),32);
    OneSpsQ = OneSpsI + ReceiverInputBatchIQSymbols;
    for(int i = 0; i < 3; i++)
    {
        SplitI[i] = (float*) _mm_malloc(BatchSize1*LengthQueue*sizeof(float),32);
        SplitQ[i] = (float*) _mm_malloc(BatchSize1*LengthQueue*sizeof(float),32);
        ViterbiOutputs[i] = (unsigned char *) _mm_malloc((BatchSize1+32),32);
        DiffDec[i] = (unsigned char *) _mm_malloc((BatchSize1*LengthQueue),32);

    }
    
    
    Outputs[0] = (unsigned char*) _mm_malloc(BatchSize3,32);




    Merged = (unsigned char *) _mm_malloc((BatchSize3),32);
    OutputAll = (unsigned char *) _mm_malloc((BatchSize3)*LengthQueue,32);
    PrbsOut = (unsigned char *) _mm_malloc((BatchSize3),32);

    #ifdef DEBUG1
    OutAllI = (float*) _mm_malloc( (NUM_SAMPLES_R1+100000)*sizeof(float),32);
    OutAllQ = (float*) _mm_malloc( (NUM_SAMPLES_R1+100000)*sizeof(float),32);
#endif



}

Receiver::~Receiver()
{
    _mm_free(FilterInI);
    _mm_free(OneSpsI);
    _mm_free(Merged);
    _mm_free(OutputAll);
    _mm_free(PrbsOut);

    for(int i = 0; i < 3; i++)
    {
        _mm_free(SplitI[i]);
        _mm_free(SplitQ[i]);
        _mm_free(ViterbiOutputs[i]);
        _mm_free(DiffDec[i]);

    }

    #ifdef DEBUG1
        _mm_free(OutAllI);
        _mm_free(OutAllQ);
    #endif

}
void Receiver::StartThreads(double RollOff, TxModes RxModeIn)
{
    RxMode = RxModeIn;
    oBufferFilter.Reset();
    objRxFilter.CreateObjects(RollOff);
    ViterbiSynchronized = false;
    PRBSSynchronized = false;
    for(int i = 0; i < 3; i++)
    {
        DemodulatorQ[i].Reset();
        DecodedQ[i].Reset();
    }

    FilterThread = std::thread(&Receiver::OperateFilter, this);
    ViterbiManagerThread = std::thread(&Receiver::OperateViterbiManager, this);

    for(int i = 0; i < 3; i++)
    {
        void *p = (void*) (IndexViterbi+i);
        ViterbiThreads[i] = std::thread(&Receiver::OperateViterbi,this,p );
    }
}
void Receiver::StopThreads(void)
{
    StopAll = true;
    CvRx2Out.notify_all();
    CvOut2Rx.notify_all();
    CvVitManager2Vit.notify_all();
    CvFilterUser.notify_all();
    if(FilterThread.joinable())
    FilterThread.join();
    for(int i = 0; i < 3;i++)
        CvViterbis2VitManager[i].notify_one();
    for(int i = 0; i < 3;i++)
    {
        if(ViterbiThreads[i].joinable())
            ViterbiThreads[i].join();
    }
    if(ViterbiManagerThread.joinable())
        ViterbiManagerThread.join();
}

void Receiver::OperateFilter(void)
{

    #ifdef WRITE_LOG_THR
	mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"Receiver OperateFilter Thread %d\n", gettid());
    fclose(fidthr);
    mtxfilethr.unlock();
    #endif

    condition_variable *pCvSamplerIn = pSampler->GetCvOut();
    while(!StopAll)
    {
        short *SamplerOut = pSampler->GetOutput(SPB*2);
        while((SamplerOut == 0) && (StopAll == false))
        {
            std::mutex mtx;
            std::unique_lock<std::mutex> lck(mtx);
            pCvSamplerIn->wait_for(lck,chrono::duration<double>(1e-3));
            #ifdef DEBUG2
                cout<<"WP A"<<endl;
            #endif
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

void Receiver::DifferentialDecode(unsigned char *In, unsigned char *Out, unsigned int Length, unsigned char &Last)
{
    In[0] = Last;
    for(int i = 0; i < Length; i+=32)
    {
        __m256i mIn = _mm256_loadu_si256((__m256i*) (In+i));
        __m256i mIn1 = _mm256_loadu_si256((__m256i*) (In+i+1));
        __m256i mxIn = _mm256_xor_si256(mIn , mIn1);
        _mm256_storeu_si256((__m256i*) (Out+i), mxIn);
    }

    Last = In[Length];
}
void Receiver::OperateViterbi(void *p)
{
    int *ip = (int*) p;
    int IndexViterbi = *ip;
    unsigned char Last = 0;


    #ifdef WRITE_LOG_THR
	mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"Receiver OperateViterbi Thread %d %d\n",IndexViterbi, gettid());
    fclose(fidthr);
    mtxfilethr.unlock();
    #endif


    while(!StopAll)
    {
        int CtrB = 0;
        while((!DemodulatorQ[IndexViterbi].AvailableRead()) && (StopAll == false))
        {
            #ifdef DEBUG2

            if(CtrB > 0)
                cout<<"WP B"<<endl;
            #endif
            std::mutex mtx;
            std::unique_lock<std::mutex> lck(mtx);
            CvVitManager2Vit.wait_for(lck,chrono::duration<double>(1e-3));
            CtrB++;
           
        }
        if(StopAll)
            break;
        int PtrRd = DemodulatorQ[IndexViterbi].GetPtrRd()*BatchSize1;
        if(!ViterbiSynchronized)
        {
            oViterbi[IndexViterbi].reset_decoder();
            oViterbi[IndexViterbi].Decode(SplitI[IndexViterbi]+PtrRd, SplitQ[IndexViterbi] + PtrRd, BatchSize1, ViterbiOutputs[IndexViterbi]+1, false, 1, 1, VitSyncResults[IndexViterbi].Metrics[0]);
            oViterbi[IndexViterbi].reset_decoder();
            oViterbi[IndexViterbi].Decode(SplitI[IndexViterbi]+PtrRd, SplitQ[IndexViterbi] + PtrRd, BatchSize1, ViterbiOutputs[IndexViterbi]+1, false, 1, -1, VitSyncResults[IndexViterbi].Metrics[1]);
            oViterbi[IndexViterbi].reset_decoder();
            oViterbi[IndexViterbi].Decode(SplitI[IndexViterbi]+PtrRd, SplitQ[IndexViterbi] + PtrRd, BatchSize1, ViterbiOutputs[IndexViterbi]+1, true, -1, 1, VitSyncResults[IndexViterbi].Metrics[2]);
            oViterbi[IndexViterbi].reset_decoder();
            oViterbi[IndexViterbi].Decode(SplitI[IndexViterbi]+PtrRd, SplitQ[IndexViterbi] + PtrRd, BatchSize1, ViterbiOutputs[IndexViterbi]+1, true,   1, 1, VitSyncResults[IndexViterbi].Metrics[3]);
            VitSyncResults[IndexViterbi].Available = true;
            DemodulatorQ[IndexViterbi].AdvanceRead();

        }
        else
        {
            oViterbi[IndexViterbi].Decode(SplitI[IndexViterbi]+PtrRd, SplitQ[IndexViterbi] + PtrRd, BatchSize1, ViterbiOutputs[IndexViterbi]+1, ViterbiParams.ExchangeIQ, ViterbiParams.SignI, ViterbiParams.SignQ, VitSyncResults[IndexViterbi].Metrics[0]);
            
            while((!DecodedQ[IndexViterbi].AvailableWrite()) && (StopAll == false))
            {
                std::mutex mtx;
                std::unique_lock<std::mutex> lck(mtx);
                std::this_thread::sleep_for((std::chrono::duration<double>(1e-4)));
            }

            int PtrWr = BatchSize1*DecodedQ[IndexViterbi].GetPtrRd();

            DifferentialDecode(ViterbiOutputs[IndexViterbi],DiffDec[IndexViterbi]+PtrWr, BatchSize1, Last);
            DemodulatorQ[IndexViterbi].AdvanceRead();

            DecodedQ[IndexViterbi].AdvanceWrite();
        }
        CvViterbis2VitManager[IndexViterbi].notify_one();

    }
}

void Receiver::OperateViterbiManager(void)
{
    NumBitsAll = 0;
    NumErrorsAll = 0;
    bool First = true;
    unsigned int PRBSSeed = 0;
    #ifdef WRITE_LOG_THR
	mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"Receiver OperateViterbiManager Thread %d\n", gettid());
    fclose(fidthr);
    mtxfilethr.unlock();
    #endif
    #ifdef DEBUG_STATISTICS
                    
    CurrDebugStatistics.MaxMetricsGrowth = 0;
    CurrDebugStatistics.SumMetricsGrowth = 0;
    CurrDebugStatistics.NumBatches = 0;
    #endif

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
            CvFilterUser.wait_for(lck,chrono::duration<double>(1e-3));
            oBufferFilter.GetReadBuffer(FilterOutI,FilterOutQ,BatchSize3*2);
        #ifdef DEBUG2
            if(FilterOutI == 0)
                cout<<"WP C"<<endl;
            #endif
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

            for(int i = 0; i < 3; i++)
            {
                int CtrD = 0;
                while((!DemodulatorQ[i].AvailableWrite()) && (StopAll==false))
                {
                
                #ifdef DEBUG2

                    CtrD++;
                    if(CtrD  >0)
                        cout<<" WP D "<<i<<endl;
                #endif
                    std::mutex mtx;
                    std::unique_lock<std::mutex> lck(mtx);
                    CvViterbis2VitManager[i].wait_for(lck,chrono::duration<double>(1e-3));
                }
                if(StopAll)
                break;
            }
            if(StopAll)
                break;
            int PtrWr = DemodulatorQ[0].GetPtrWr();
            Split3(OneSpsI, OneSpsQ, ReceiverInputBatchIQSymbols,PtrWr*BatchSize1);
            if(StopAll)
                break;
            if(ViterbiSynchronized == false)
            {
                for(int i = 0; i < 3; i++)
                {
                    VitSyncResults[i].Available = false;
                }
            }
            for(int i = 0; i < 3; i++)
            {
                DemodulatorQ[i].AdvanceWrite();
            }
            CvVitManager2Vit.notify_all();
            
            if(StopAll)
                break;
            if(!ViterbiSynchronized)
            {
                for(int i = 0; i < 3; i++)
                {
                    int CtrE = 0;
                    while((!VitSyncResults[i].Available) && (StopAll == false))
                    {
                        #ifdef DEBUG2
                        CtrE++;
                        if(CtrE > 0)
                            cout<<" WP E "<<i<<endl;
                        #endif
                        std::mutex mtx;
                        std::unique_lock<std::mutex> lck(mtx);
                        CvViterbis2VitManager[i].wait_for(lck,chrono::duration<double>(1e-3));
                    }
                    if(StopAll)
                        break;
                }
                if(StopAll)
                    break;
                //Test if Synchronization
                int BestIndex[3];
                bool PassThresh[3];
                __m128 mThresh = _mm_set1_ps(ViterbiThreshold1);
                bool Success = true;
                for(int i = 0; i < 3; i++)
                {
                    BestIndex[i] = 0;
                    float BestMetrics = VitSyncResults[i].Metrics[0];
                    for(int j = 1; j < 4; j++)
                    {
                        if(VitSyncResults[i].Metrics[j] < BestMetrics)
                        {
                            BestIndex[i] = j;
                            BestMetrics = VitSyncResults[i].Metrics[j];
                        }
                    }
                    VitSyncResults[i].Metrics[BestIndex[i]] = 1e7;
                    __m128 mResults = _mm_loadu_ps(VitSyncResults[i].Metrics);
                    __m128 mBest = _mm_set1_ps(BestMetrics);
                    mResults = _mm_sub_ps(mResults,mBest );
                    __m128 mCmp = _mm_cmp_ps(mResults, mThresh, _CMP_GE_OS);
                    int PassThresh =  _mm_movemask_ps(mCmp);
                    if(PassThresh != 15)
                    {
                        Success = false;
                        break;
                    }

                    
                }

                if(Success)
                {
                    if((BestIndex[0]== BestIndex[1]) && (BestIndex[0]== BestIndex[2]))
                    {
                        switch (BestIndex[0]) {
                        case 0:
                            ViterbiParams.ExchangeIQ = false;
                            ViterbiParams.SignI = 1;
                            ViterbiParams.SignQ = 1;
                            break;
                        case 1:
                            ViterbiParams.ExchangeIQ = false;
                            ViterbiParams.SignI = 1;
                            ViterbiParams.SignQ = -1;
                            break;
                        case 2:
                            ViterbiParams.ExchangeIQ = true;
                            ViterbiParams.SignI = -1;
                            ViterbiParams.SignQ = 1;
                            break;
                        case 3:
                            ViterbiParams.ExchangeIQ = true;
                            ViterbiParams.SignI = 1;
                            ViterbiParams.SignQ = 1;
                            break;   
                        }
                        ViterbiSynchronized = true;
                        for(int i = 0; i < 3; i++)
                            oViterbi[i].reset_decoder();
                    }
                }
            }
            else 
            {
                for(int i = 0; i < 3; i++)
                {
                    int CtrE = 0;
                    while((!DecodedQ[i].AvailableRead()) && (StopAll==false))
                    {
                    #ifdef DEBUG2

                        CtrE++;
                        if(CtrE > 0)
                            cout<<" WP F "<<i<<endl;
                    #endif
                        std::mutex mtx;
                        std::unique_lock<std::mutex> lck(mtx);
                        CvViterbis2VitManager[i].wait_for(lck,chrono::duration<double>(1e-3));
                    }
                    if(StopAll)
                    break;
                }

                #ifdef DEBUG_STATISTICS
                for(int i = 0; i < 3; i++)
                {
                    if(CurrDebugStatistics.MaxMetricsGrowth < VitSyncResults[i].Metrics[0])
                        CurrDebugStatistics.MaxMetricsGrowth = VitSyncResults[i].Metrics[0];
                    CurrDebugStatistics.SumMetricsGrowth += VitSyncResults[i].Metrics[0];
                }
                CurrDebugStatistics.NumBatches += 3;
                CurrDebugStatistics.MeanMetricsGrowth =  CurrDebugStatistics.SumMetricsGrowth/double(CurrDebugStatistics.NumBatches);
                #endif

                if(StopAll)
                    break;
                int PtrRd = BatchSize1*DecodedQ[0].GetPtrRd();
                int PtrOut = 0;
                for(int i = 0; i < BatchSize1; i++)
                {
                    Merged[PtrOut++] = DiffDec[0][PtrRd+i];
                    Merged[PtrOut++] = DiffDec[1][PtrRd+i];
                    Merged[PtrOut++] = DiffDec[2][PtrRd+i];
                }

                for(int i = 0; i < 3; i++)
                {
                    DecodedQ[i].AdvanceRead();
                }

                 //Descramble
                 int CtrE = 0;
                while((!OutputQ.AvailableWrite()) && (StopAll==false))
                {
                #ifdef DEBUG2

                    CtrE++;
                    if(CtrE > 0)
                        cout<<" WP G "<<endl;
                #endif
                    std::mutex mtx;
                    std::unique_lock<std::mutex> lck(mtx);
                    CvOut2Rx.wait_for(lck,chrono::duration<double>(1e-3));
                }
                if(StopAll)
                    break;
                PtrWr = BatchSize3*OutputQ.GetPtrWr();
                objDescrambler.Descramble(Merged, OutputAll+PtrWr, BatchSize3);
                OutputQ.AdvanceWrite();

                if(RxMode == FILE_TX)
                {
                    CvRx2Out.notify_one();
                }
                else {
                
                    int PtrRd = PtrWr;
                    if(!PRBSSynchronized)
                    {
                        int PtrStart;
                        PRBSSynchronized = SyncPRBS(OutputAll+PtrRd, PtrStart, PRBSSeed);
                        if(PRBSSynchronized)
                        {
                            oPrbs.CreateOutputs(PRBSSeed, BatchSize3 - PtrStart, PrbsOut);
                            CountErrors(OutputAll+PtrRd+PtrStart, PrbsOut, BatchSize3-PtrStart);
                        }
                    }
                    else {
                        
                            oPrbs.CreateOutputs(PRBSSeed, BatchSize3, PrbsOut);
                            CountErrors(OutputAll+PtrRd, PrbsOut, BatchSize3);

                    }
                    OutputQ.AdvanceRead();
                }

            }

                 

      
        }
        oBufferFilter.AdvancePtrRd(BatchSize3*2);

    }

}
void  Receiver:: Split3(float *InputI, float *InputQ, int Length, int PtrWr)
{
    __m256i mIndex[3];
    mIndex[0] = _mm256_setr_epi32(0, 3, 6, 9, 12, 15, 18, 21);
    mIndex[1] = _mm256_setr_epi32(1, 4, 7, 10, 13, 16, 19, 22);
    mIndex[2] = _mm256_setr_epi32(2, 5, 8, 11, 14, 17, 20, 23);
    //mIndex[3] = _mm256_set1_epi32(24);
    int PtrOut = PtrWr;
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
void Receiver::CountErrors(unsigned char *Input, unsigned char *Template, int Length)
{
    __m256i mOnes = _mm256_set1_epi8(1);
    for(int j = 0; j < Length; j+=32)
    {
        __m256i mA = _mm256_loadu_si256((__m256i*) (Input+j));
        __m256i mB = _mm256_loadu_si256((__m256i*) (Template+j));
        mA = _mm256_xor_si256(mA, mB);
        mA = _mm256_cmpeq_epi8(mA,mOnes);
        int cmp = _mm256_movemask_epi8(mA);
        int NumErrors1 = _mm_popcnt_u32(cmp);
        NumErrorsAll += NumErrors1;
    }
    NumBitsAll += Length;
}
bool Receiver::SyncPRBS(unsigned char *In, int &PtrStart,unsigned int &Seed)
{
    bool Success = false;


    int Ptr = BatchSize1 - 23;
    while(!Success)
    {
        if((Ptr + 23 + 128) >= BatchSize3)
            break; 
        Seed = 0;
        int i;
        for ( i = 0; i < 23; i++)
        {
            Seed = (Seed << 1) | (unsigned int) In[Ptr+i];
        }
        int Ptr1 = Ptr + i;
        
        oPrbs.CreateOutputs(Seed,128,PrbsOut);

        __m256i mOnes = _mm256_set1_epi8(1);
        int NumErrors = 0;
        for(int j = 0; j < 128; j+=32)
        {
            __m256i mA = _mm256_loadu_si256((__m256i*) (In+Ptr1+j));
            __m256i mB = _mm256_loadu_si256((__m256i*) (PrbsOut+j));
            mA = _mm256_xor_si256(mA, mB);
            mA = _mm256_cmpeq_epi8(mA,mOnes);
            int cmp = _mm256_movemask_epi8(mA);
            int NumErrors1 = _mm_popcnt_u32(cmp);
            NumErrors += NumErrors1;
        }
        if(NumErrors <= PRBSThreshold)
        {
            Success = true;
            PtrStart = Ptr1 + 128;
        }
        else {
        
            Ptr+= 128;
        }
    }
    return Success;
}