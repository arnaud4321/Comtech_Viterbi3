#include "Transmitter.h"
extern mutex mtxfilethr;
Transmitter::Transmitter(/* args */):DataQ(LengthQueue),FilterQ(LengthQueue)
{
    ScrambledOut = (unsigned char*) _mm_malloc(BatchSize3+32,32);
    SplitOut[0] = (unsigned char*) _mm_malloc(BatchSize3>>3,32);
    SplitOut[1] = SplitOut[0] + (BatchSize1>>3);
    SplitOut[2] = SplitOut[1] + (BatchSize1>>3);
    DiffOut[0] = (unsigned char*) _mm_malloc(BatchSize3>>3,32);
    DiffOut[1] = DiffOut[0] + (BatchSize1>>3);
    DiffOut[2] = DiffOut[1] + (BatchSize1>>3);
    ConvOut[0] = (unsigned char*) _mm_malloc(BatchConvSize1*3*LengthQueue,32);
    for(int i = 1; i < (3*LengthQueue); i++)
        ConvOut[i] = ConvOut[i-1] + BatchConvSize1;
    FilterInI = (float*) _mm_malloc(BatchSize3*sizeof(float),32);
    FilterInQ = (float*) _mm_malloc(BatchSize3*sizeof(float),32);
    FilterOutI = (float*) _mm_malloc(BatchSize3*2*sizeof(float),32);
    FilterOutQ = (float*) _mm_malloc(BatchSize3*2*sizeof(float),32);
    TxOut[0] = (float*) _mm_malloc(LengthQueue*BatchSize3*4*sizeof(float),32);
    for(int i = 1; i < LengthQueue; i++)
    {
        TxOut[i] = TxOut[i-1] + BatchSize3*4;        
    }  
    #ifdef DEBUG_TX
        OutAllI = (float *) _mm_malloc( 20000000*sizeof(float),32);
    #endif

}

Transmitter::~Transmitter()
{
    _mm_free(ScrambledOut);
    _mm_free(SplitOut[0]);
    _mm_free(DiffOut[0]);
    _mm_free(ConvOut[0]);
    _mm_free(FilterInI);
    _mm_free(FilterInQ);
    _mm_free(FilterOutI);
    _mm_free(FilterOutQ);
    _mm_free(TxOut[0]);
    #ifdef DEBUG_TX
    _mm_free(OutAllI);

    #endif
}


void Transmitter::StartThreads(TxModes TxMode, float RollOff, string FileName )
{
    StopAll = false;
    if(TxMode == PRBS_TX)
    {
        Data = (unsigned char*) _mm_malloc(8388608,32);
        unsigned int state = 1;
        oPrbs.CreateOutputs(state,8388607,Data);
        LengthData = 8388607;
    }
    else
    {
        FILE *fid = fopen(FileName.c_str(),"rb");
        fseek(fid,0,SEEK_END);
        LengthData = ftell(fid);
        fseek(fid,0,SEEK_SET);
        fread(Data,1,LengthData,fid);
    }
    DataQ.Reset();
    DataGenThread = std::thread(&Transmitter::GenerateData, this);
    oTxFilter.CreateObjects(RollOff);
    FilterQ.Reset();
    FilterThread = std::thread(&Transmitter::FilterData, this);
}

void Transmitter::StopThreads(void)
{
    StopAll = true;
    CvFilterData.notify_one();
    CvDataFilter.notify_one();
    CvOutFilter.notify_one();
    CvFilterOut.notify_one();
    
    if(DataGenThread.joinable())
        DataGenThread.join();
    if(FilterThread.joinable())
        FilterThread.join();
    _mm_free(Data);
}

void Transmitter::GenerateData(void)
{

     #ifdef WRITE_LOG_THR
	mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"Transmitter GenerateData Thread %d\n", gettid());
    fclose(fidthr);
    mtxfilethr.unlock();
    #endif

    unsigned int PtrData = 0;
    while(!StopAll)
    {

        //Scramble Data

        unsigned int PtrEnd = PtrData + BatchSize3;
        if(PtrEnd < LengthData)
        {
            oScrambler.Scramble(Data+PtrData,ScrambledOut,BatchSize3);
            PtrData += BatchSize3;
        }
        else
        {
            int First = LengthData - PtrData;
            int Second = BatchSize3 - First;
            oScrambler.Scramble(Data+PtrData,ScrambledOut,First);
            oScrambler.Scramble(Data,ScrambledOut+First,Second);
            PtrData = Second;
        }
        Split3();//Split the scrambled data into 3
        for(int i = 0; i < 3; i++)
        {
            oDiffEncode[i].CreateOutputs(SplitOut[i],DiffOut[i],BatchSize1Bytes);
        }
        if(StopAll)
            break;

        while((!DataQ.AvailableWrite()) && (StopAll==false))
        {
            std::mutex mtx;
            std::unique_lock<std::mutex> lck(mtx);
            CvFilterData.wait(lck);
        }
        if(StopAll)
            break;
        unsigned int PtrWr = DataQ.GetPtrWr();
        for(int i = 0; i < 3; i++)
        {
            oConvEncoder[i].Encode(DiffOut[i],ConvOut[3*PtrWr+i],BatchSize1Bytes);
        }
        DataQ.AdvanceWrite();
        CvDataFilter.notify_one();

    }
        

}

void Transmitter::FilterData(void)
{
 #ifdef WRITE_LOG_THR
	mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"Transmitter FilterData Thread %d\n", gettid());
    fclose(fidthr);
    mtxfilethr.unlock();
    #endif

    while(!StopAll)
    {
        while((!DataQ.AvailableRead()) && (StopAll==false))
        {
            std::mutex mtx;
            std::unique_lock<std::mutex> lck(mtx);
            CvDataFilter.wait(lck);
        }
        if(StopAll)
            break;
        
       

        //Interleave the data and map to +/- 1    
        unsigned int PtrRd = DataQ.GetPtrRd();
        unsigned int PtrOut = 0;
        for(int i = 0; i < (2*BatchSize1); i+=2 )
        {
            for(int j = 0; j < 3; j++)
            {
                FilterInI[PtrOut] = Map[ConvOut[3*PtrRd+j][i]];
                FilterInQ[PtrOut++] = Map[ConvOut[3*PtrRd+j][i+1]];
            }
        }
        
        DataQ.AdvanceRead();
        CvFilterData.notify_one();

        oTxFilter.CreateOutputs(FilterInI,FilterInQ,FilterOutI,FilterOutQ,BatchSize3);
        
        while((!FilterQ.AvailableWrite()) && (StopAll==false))
        {
            std::mutex mtx;
            std::unique_lock<std::mutex> lck(mtx);
            CvOutFilter.wait(lck);
        }
        if(StopAll)
            break;

        unsigned int PtrWr = FilterQ.GetPtrWr();
        interleave_iq_f32_to_f32_avx2(FilterOutI,FilterOutQ,TxOut[PtrWr],BatchSize3*2);
        #ifdef DEBUG_TX
        std::copy(TxOut[PtrWr],TxOut[PtrWr]+BatchSize3*4,OutAllI+PtrOutAll);
        PtrOutAll += BatchSize3*4;
        if(PtrOutAll >= 5000000)
        {    
            FILE *fid = fopen("TxOut.bin","wb");
            int size1 = sizeof(float);
            fwrite(OutAllI,sizeof(float),PtrOutAll,fid);
            fclose(fid);
            exit(-1);
        }
        #endif
        FilterQ.AdvanceWrite();
        CvFilterOut.notify_one();
    
    }

    oTxFilter.DeleteObjects();
}



void Transmitter::Split3(void)
{
    alignas (32) unsigned char Temp[32];
    int PtrOut = 0;
    unsigned int k = 0;
    for(unsigned int i = 0; i < BatchSize3; i+= 24)
    {
        
        for(int j = 0; j < 8 ; j++)
        {
            Temp[j ] = ScrambledOut[k++];
            Temp[8 + j] = ScrambledOut[k++];
            Temp[16 +j] = ScrambledOut[k++];            
        }

        __m256i mIn = _mm256_load_si256(( __m256i*) Temp);
        mIn = _mm256_cmpgt_epi8(mIn, _mm256_setzero_si256());
        int x = _mm256_movemask_epi8(mIn);
        SplitOut[0][PtrOut] = x & 0xff;
        SplitOut[1][PtrOut] = (x>>8) & 0xff;
        SplitOut[2][PtrOut] = (x>>16) & 0xff;
        PtrOut ++;
    }
}

float* Transmitter::GetOutput(void)
{
    float *Out = 0;

	if(FilterQ.AvailableRead())
    {
        unsigned int PtrRd = FilterQ.GetPtrRd();
        Out =TxOut[PtrRd];
        
    }
   
    return Out;
}

void Transmitter::interleave_iq_f32_to_f32_avx2(const float* xi, const float* xq, float* y,int n)   // n = number of I/Q pairs, multiple of 8
{
    int PtrOut = 0;
    for (int i = 0; i < n; i += 8)
    {
        // 1) Load
        __m256 vi = _mm256_loadu_ps(xi + i);
        __m256 vq = _mm256_loadu_ps(xq + i);
        //interleave
        __m256 lo = _mm256_unpacklo_ps(vi, vq); // I0 Q0 I1 Q1 I4 Q4 I5 Q5
        __m256 hi = _mm256_unpackhi_ps(vi, vq); // I2 Q42 I3 Q3 I6 Q6 I7 Q7

        __m256 lo1 = _mm256_permute2f128_ps(lo,hi,0x20);
        __m256 hi1 = _mm256_permute2f128_ps(lo,hi,0x31);
        
        // 5) Store
        _mm256_storeu_ps(y + PtrOut, lo1);
        PtrOut += 8;
        _mm256_storeu_ps(y + PtrOut, hi1);
        PtrOut += 8;
        
    }
}
void Transmitter::interleave_iq_f32_to_i16_avx2(const float* xi, const float* xq, short* y,int n)   // n = number of I/Q pairs, multiple of 8
{
    for (int i = 0; i < n; i += 8) {
        // 1) Load
        __m256 vi = _mm256_loadu_ps(xi + i);
        __m256 vq = _mm256_loadu_ps(xq + i);

        // 2) Convert float -> int32 (round-to-nearest)
        __m256i ii = _mm256_cvtps_epi32(vi);
        __m256i iq = _mm256_cvtps_epi32(vq);

        // 3) Interleave int32
        __m256i lo = _mm256_unpacklo_epi32(ii, iq); // I0 Q0 I1 Q1 I2 Q2 I3 Q3
        __m256i hi = _mm256_unpackhi_epi32(ii, iq); // I4 Q4 I5 Q5 I6 Q6 I7 Q7

        // 4) Pack int32 -> int16
        __m256i packed = _mm256_packs_epi32(lo, hi);

        // 5) Store
        _mm256_storeu_si256((__m256i*)(y + 2*i), packed);
    }
}



