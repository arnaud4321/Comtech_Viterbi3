/**
 * @file Transmitter.h
 * @brief Multithreaded QPSK-oriented transmitter: source data, coding, scrambling, pulse shaping.
 */
#pragma once
#include <immintrin.h>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <string>
#include <vector>
using namespace std;
#include "definitions.h"
#include "Prbs23.h"
#include "ConvEncoder.h"
#include "SelfSyncScrambler_V35.h"
#include "DiffEncode.h"
#include "SimpleQueue.h"
#include "TxFilter.h"

//#define DEBUG_TX

/**
 * @brief Transmit chain from bits to band-limited IQ floats consumed by @ref AWGNChannel.
 *
 * @details **Data path (conceptual):**
 * - **Source:** PRBS-23 or file bytes; mapped to uncoded bits for the parallel trellis branches.
 * - **Scrambling:** ITU-T V.35 self-synchronous scrambler (@ref SelfSyncScrambler_V35) on the data stream.
 * - **Split:** @c Split3() reorganizes scrambled bytes into three parallel bit-packed streams so each
 *   @ref ConvEncoder / @ref Viterbi branch corresponds to a different alignment hypothesis at the receiver.
 * - **Differential encoding:** @ref DiffEncode on each branch before convolution.
 * - **Convolutional coding:** Three @ref ConvEncoder instances (rate 1/2, K=7-style polynomials); outputs
 *   are multiplexed to the pulse shaper as interleaved I/Q symbol triples.
 * - **Pulse shaping:** @ref TxFilter (root raised cosine, configurable roll-off); output is complex
 *   baseband at the simulation sample rate, queued for the channel.
 *
 * **Threads (see @ref Transmitter.cpp):** @c GenerateData — scrambler, @c Split3, @c DiffEncode per branch,
 * @c ConvEncoder, hand-off via @c DataQ. @c FilterData — maps coded bits to @f$\pm1@f$, @ref TxFilter, interleaved
 * float batches into @c FilterQ for @c CopyOutputSamples / @ref AWGNChannel.
 */
class Transmitter
{
    private:

#ifdef DEBUG_TX
    float *OutAllI = 0;
    unsigned int PtrOutAll = 0;
#endif

    float Map[2] = {-1,1};
    unsigned char *Data = 0;
    unsigned char *ScrambledOut = 0;
    unsigned char *SplitOut[3];
    unsigned char *DiffOut[3];
    static constexpr int LengthQueue = 4;

    unsigned char *ConvOut[LengthQueue*3];
    float *FilterInI = 0, *FilterInQ = 0;
    float *FilterOutI = 0, *FilterOutQ = 0;
    float *TxOut[LengthQueue];
    unsigned int LengthData = 0;

    Prbs23 oPrbs;
    SelfSyncScrambler_V35 oScrambler;
    ConvEncoder oConvEncoder[3];
    DiffEncode oDiffEncode[3];
    TxFilter oTxFilter;
    thread DataGenThread, FilterThread;
    std::atomic<bool> StopAll{false};
    SimpleQueue DataQ, FilterQ;
    condition_variable CvFilterData, CvDataFilter, CvOutFilter, CvFilterOut;
    std::mutex mtxDataQ_;
    std::mutex mtxFilterQ_;

    // TX symbol frame queue (pre-TxFilter, mapped to +/-1) with backpressure.
    // Only built when ENABLE_RAW_PREVITERBI_METRICS is enabled.
#ifdef ENABLE_RAW_PREVITERBI_METRICS
    static constexpr int kTxSymFrame = ReceiverInputBatchIQSymbols; // == BatchSize3
    static constexpr int kTxSymQueueDepth = 2048;
    struct TxSymSlot
    {
        alignas(32) float i[kTxSymFrame];
        alignas(32) float q[kTxSymFrame];
    };
    std::vector<TxSymSlot> txSymPool_;
    int txSymHead_ = 0;
    int txSymTail_ = 0;
    int txSymCount_ = 0;
    std::mutex mtxTxSymQ_;
    std::condition_variable cvTxSymReady_;
    std::condition_variable cvTxSymSpace_;
    bool pushTxSymFrameBlocking_(const float* i, const float* q);
#endif

    void GenerateData(void);
    void FilterData(void);
    void Split3(void);
    void interleave_iq_f32_to_f32_avx2(const float* xi, const float* xq, float* y, int n) ;
    void interleave_iq_f32_to_i16_avx2(const float* xi, const float* xq, short* y, int n) ;
public:
    Transmitter(/* args */);
    ~Transmitter();
    void StartThreads(TxModes TxMode, float RollOff, string FileName = "");
    void StopThreads(void);
	float *GetOutput(void);
    void AdvanceOut(void)
    {
        {
            std::lock_guard<std::mutex> lk(mtxFilterQ_);
            FilterQ.AdvanceRead();
        }
        CvOutFilter.notify_all();
    }
    // Atomic copy (internal mutex) from the TX output FIFO (FilterQ).
    // Returns false if stopAll becomes true while waiting.
    bool CopyOutputSamples(float* dst, int nFloats, std::atomic<bool>& stopAll);

    // TX symbol reference accessors.
    // When ENABLE_RAW_PREVITERBI_METRICS is OFF, these return false and do nothing (no backpressure / no queue).
    bool TryPopTxSymbolFrame(float* dstI, float* dstQ, int nSym);
    bool WaitPopTxSymbolFrame(float* dstI, float* dstQ, int nSym);
    void GetCVOut(condition_variable* &CVIn, condition_variable* &CVOut )
    {
        CVIn = &CvOutFilter;
        CVOut = &CvFilterOut;
    }
    double GetTxPower()
    {
        return oTxFilter.TxPower;
    }
};