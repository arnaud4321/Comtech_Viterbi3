/**
 * @file Transmitter.cpp
 * @brief Multithreaded TX path: scramble → three parallel branches (split/diff/conv) → RRC pulse shaping.
 *
 * @details **Threads**
 * - **GenerateData** — Ring @c DataQ: V35 scrambler on source bytes, @c Split3 into three packed streams,
 *   @c DiffEncode and @c ConvEncode per branch, enqueue one coded slot per batch for the filter thread.
 * - **FilterData** — Dequeue coded slot, map bits to @f$\pm1@f$, build interleaved I/Q for @ref TxFilter
 *   (root raised cosine, @c TxNSS samples per symbol), enqueue float batches in @c FilterQ consumed via
 *   @c CopyOutputSamples by @ref AWGNChannel.
 *
 * **StartThreads** builds PRBS or loads @c FileName into @c Data, calls @c oTxFilter.CreateObjects(RollOff),
 * resets queues, starts both workers. **StopThreads** signals @c StopAll, joins threads, frees @c Data.
 */

#include "Transmitter.h"
#include <iostream>
#include <cstring>
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

#ifdef ENABLE_RAW_PREVITERBI_METRICS
    txSymPool_.resize(static_cast<size_t>(kTxSymQueueDepth));
#endif

}

#ifdef ENABLE_RAW_PREVITERBI_METRICS
bool Transmitter::pushTxSymFrameBlocking_(const float* i, const float* q)
{
    std::unique_lock<std::mutex> lk(mtxTxSymQ_);
    while (!StopAll && txSymCount_ >= kTxSymQueueDepth)
    {
        cvTxSymSpace_.wait_for(lk, std::chrono::milliseconds(1));
    }
    if (StopAll)
        return false;
    TxSymSlot& slot = txSymPool_[static_cast<size_t>(txSymTail_)];
    std::memcpy(slot.i, i, sizeof(float) * static_cast<size_t>(kTxSymFrame));
    std::memcpy(slot.q, q, sizeof(float) * static_cast<size_t>(kTxSymFrame));
    txSymTail_ = (txSymTail_ + 1) % kTxSymQueueDepth;
    txSymCount_++;
    lk.unlock();
    cvTxSymReady_.notify_one();
    return true;
}

bool Transmitter::TryPopTxSymbolFrame(float* dstI, float* dstQ, int nSym)
{
    if (nSym != kTxSymFrame)
        return false;
    std::lock_guard<std::mutex> lk(mtxTxSymQ_);
    if (StopAll || txSymCount_ <= 0)
        return false;
    const TxSymSlot& slot = txSymPool_[static_cast<size_t>(txSymHead_)];
    std::memcpy(dstI, slot.i, sizeof(float) * static_cast<size_t>(kTxSymFrame));
    std::memcpy(dstQ, slot.q, sizeof(float) * static_cast<size_t>(kTxSymFrame));
    txSymHead_ = (txSymHead_ + 1) % kTxSymQueueDepth;
    txSymCount_--;
    cvTxSymSpace_.notify_one();
    return true;
}

bool Transmitter::WaitPopTxSymbolFrame(float* dstI, float* dstQ, int nSym)
{
    if (nSym != kTxSymFrame)
        return false;
    std::unique_lock<std::mutex> lk(mtxTxSymQ_);
    while (!StopAll && txSymCount_ <= 0)
    {
        cvTxSymReady_.wait_for(lk, std::chrono::milliseconds(1));
    }
    if (StopAll)
        return false;
    const TxSymSlot& slot = txSymPool_[static_cast<size_t>(txSymHead_)];
    std::memcpy(dstI, slot.i, sizeof(float) * static_cast<size_t>(kTxSymFrame));
    std::memcpy(dstQ, slot.q, sizeof(float) * static_cast<size_t>(kTxSymFrame));
    txSymHead_ = (txSymHead_ + 1) % kTxSymQueueDepth;
    txSymCount_--;
    lk.unlock();
    cvTxSymSpace_.notify_one();
    return true;
}
#endif // ENABLE_RAW_PREVITERBI_METRICS

#ifndef ENABLE_RAW_PREVITERBI_METRICS
bool Transmitter::TryPopTxSymbolFrame(float* dstI, float* dstQ, int nSym)
{
    (void)dstI;
    (void)dstQ;
    (void)nSym;
    return false;
}

bool Transmitter::WaitPopTxSymbolFrame(float* dstI, float* dstQ, int nSym)
{
    (void)dstI;
    (void)dstQ;
    (void)nSym;
    return false;
}
#endif // !ENABLE_RAW_PREVITERBI_METRICS



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
#ifdef ENABLE_RAW_PREVITERBI_METRICS
    {
        std::lock_guard<std::mutex> lk(mtxTxSymQ_);
        txSymHead_ = 0;
        txSymTail_ = 0;
        txSymCount_ = 0;
    }
#endif
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
    std::cout << "[Transmitter] StopThreads called." << std::endl;
    StopAll = true;
    std::cout << "[Transmitter] Notifying CVs..." << std::endl;
    CvFilterData.notify_all();
    CvDataFilter.notify_all();
    CvOutFilter.notify_all();
    CvFilterOut.notify_all();
#ifdef ENABLE_RAW_PREVITERBI_METRICS
    cvTxSymSpace_.notify_all();
    cvTxSymReady_.notify_all();
#endif
    
    std::cout << "[Transmitter] Joining DataGenThread..." << std::endl;
    if(DataGenThread.joinable())
        DataGenThread.join();
    std::cout << "[Transmitter] Joining FilterThread..." << std::endl;
    if(FilterThread.joinable())
        FilterThread.join();
    std::cout << "[Transmitter] Threads joined, freeing Data..." << std::endl;
    _mm_free(Data);
    std::cout << "[Transmitter] StopThreads finished." << std::endl;
}

/**
 * @brief Scramble/split/diff/conv pipeline → @c DataQ for @c FilterData.
 */
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

        std::unique_lock<std::mutex> lk(mtxDataQ_);
        while (!StopAll && !DataQ.AvailableWrite())
        {
            CvFilterData.wait_for(lk, std::chrono::milliseconds(1));
        }
        if(StopAll)
            break;
        unsigned int PtrWr = DataQ.GetPtrWr();
        for(int i = 0; i < 3; i++)
        {
            oConvEncoder[i].Encode(DiffOut[i],ConvOut[3*PtrWr+i],BatchSize1Bytes);
        }
        DataQ.AdvanceWrite();
        lk.unlock();
        CvDataFilter.notify_all();
    }
        

}

/**
 * @brief Map coded bits to symbols, @ref TxFilter, write @c FilterQ / @c TxOut FIFO.
 */
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
        std::unique_lock<std::mutex> lkData(mtxDataQ_);
        while (!StopAll && !DataQ.AvailableRead())
        {
            CvDataFilter.wait_for(lkData, std::chrono::milliseconds(1));
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

        // Store the ideal (noiseless) mapped symbols into the TX ring (pre-Viterbi SER / raw BER).
        // This is before pulse shaping and channel impairments.
#ifdef ENABLE_RAW_PREVITERBI_METRICS
        (void)pushTxSymFrameBlocking_(FilterInI, FilterInQ);
#endif
        
        DataQ.AdvanceRead();
        lkData.unlock();
        CvFilterData.notify_all();

        oTxFilter.CreateOutputs(FilterInI,FilterInQ,FilterOutI,FilterOutQ,BatchSize3);
        
        std::unique_lock<std::mutex> lkFilt(mtxFilterQ_);
        while (!StopAll && !FilterQ.AvailableWrite())
        {
            CvOutFilter.wait_for(lkFilt, std::chrono::milliseconds(1));
        }
        if(StopAll)
            break;

        unsigned int PtrWr = FilterQ.GetPtrWr();
        interleave_iq_f32_to_f32_avx2(FilterOutI,FilterOutQ,TxOut[PtrWr],BatchSize3*2);
        FilterQ.AdvanceWrite();
        lkFilt.unlock();
        CvFilterOut.notify_all();
    }

    oTxFilter.DeleteObjects();
}

bool Transmitter::CopyOutputSamples(float* dst, int nFloats, std::atomic<bool>& stopAll)
{
    std::unique_lock<std::mutex> lk(mtxFilterQ_);
    while (!stopAll && !StopAll && !FilterQ.AvailableRead())
    {
        CvFilterOut.wait_for(lk, std::chrono::milliseconds(1));
    }
    if (stopAll || StopAll)
        return false;
    if (!FilterQ.AvailableRead())
        return false;
    unsigned int PtrRd = FilterQ.GetPtrRd();
    float* p = TxOut[PtrRd];
    std::memcpy(dst, p, sizeof(float) * static_cast<size_t>(nFloats));
    FilterQ.AdvanceRead();
    lk.unlock();
    CvOutFilter.notify_one();
    return true;
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
    std::lock_guard<std::mutex> lk(mtxFilterQ_);
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



