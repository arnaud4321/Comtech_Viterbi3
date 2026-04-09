/**
 * @file AWGNChannel.cpp
 * @brief Implementation of the multithreaded AWGN channel (noise mix, Es/N0, profile, NCO, SCO, shorts out).
 *
 * @details **Pipeline (float domain until @ref ChannelSamplingClockOffset writes @c short IQ into
 * @c AWGNChannel::OutputBuffer):**
 *
 * 1. **GenerateNoise** — Fills a FIFO of AVX-aligned Gaussian buffers (`randn`) consumed by the output thread.
 * 2. **GenerateOutput** — For each batch: copies transmitter IQ (@ref Transmitter::CopyOutputSamples), combines with
 *    noise as `mkSig * tx + Stdn * noise` (AVX). **Es/N0** is enforced via `Stdn` and `mkSig` computed in
 *    @c AWGNChannel::StartThreads() from @ref Params::EsN0, TX average power, a reference constant `Pow0`, and a
 *    **Backoff** (dB) so that signal+noise scaled together stay below int16 headroom after downstream processing.
 * 3. Still in **GenerateOutput** — Applies a **time-varying amplitude** `10^(rangeDb/20)` where `rangeDb` follows
 *    the same four-segment piecewise schedule as frequency (stable → ramp → stable → ramp), driven by
 *    *simulation time* `t_rate = batchCount * TxOutputBatchDuration`. Pushes a @c FreqChunk (interleaved floats +
 *    metadata) onto `freqQ_`.
 * 4. **ApplyFrequencyOffset** — Dequeues each chunk: sets NCO frequency (@ref FrequencyOffset), updates total SCO ppm
 *    (`ClockMismatchPpm` + Doppler term from carrier offset / `CarrierToSymbolRateRatio` / `SymbolRate`),
 *    publishes @c GetCurrentApplied() atomics, rotates IQ in place, then forwards interleaved floats to
 *    @ref ChannelSamplingClockOffset::EnqueueNoisyInterleaved.
 * 5. **ChannelSamplingClockOffset** (started from @c AWGNChannel::StartThreads()) — Resamples to emulate Rx clock
 *    error and writes saturated **short** IQ into @c OutputBuffer; @c CopyOutputSamples() is the consumer API for
 *    @ref Sampler.
 *
 * **Stop condition:** after batch index reaches `TransitionCounter[4]` (end of the fourth segment), sets the
 * global @c Finish flag and stops. **Thread order on shutdown:** @c AWGNChannel::StopThreads() joins noise/output/freq
 * threads then stops the SCO worker.
 */

#include "AWGNChannel.h"
#include "Transmitter.h"
#include <cstring>
#include <deque>
#include <iomanip>
#include <iostream>
#include <atomic>
extern std::atomic<bool> Finish;
extern mutex mtxfilethr;
AWGNChannel::AWGNChannel(unsigned int Seed)
    : NoiseQ(LengthQueue)
    , OutputBuffer(kChannelOutShortRingSize, kChannelOutShortRingExtra)
{
    // NoiseBatchSize: slightly > one TX batch so slowed-time runs do not run out of noise; 16-aligned for AVX loops.
    oNoiseGen.set_seed(Seed);
    NoiseBatchSize = (int) ((double) TxOutputBatchSize * 1.001); //make sure there will be enough noise even if time is slowed down
    int Tmp = (NoiseBatchSize >> 4)<<4; //make a length of 16
    if(NoiseBatchSize > Tmp)
        Tmp += 16; 
    NoiseBatchSize = Tmp;
    Noise[0] = (float *) _mm_malloc(sizeof(float) * NoiseBatchSize * LengthQueue,32);
    for(int i = 1; i < LengthQueue; i++)
    {
        // Each slot must be spaced by NoiseBatchSize (randn writes NoiseBatchSize floats).
        // Using TxOutputBatchSize here would overlap slots and corrupt memory when NoiseBatchSize > TxOutputBatchSize.
        Noise[i] = Noise[i-1] + NoiseBatchSize;
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
    // Reference symbol energy convention: Pow0 ties Es/N0 to the noise variance used in combination with TX power.
    double Pow0 = 2.0;//constant
    EsN0db = pParams->EsN0;
    double N0 = Pow0*pow(10.0,-0.1*EsN0db);
    // Digital peak target for I^2+Q^2 after scaling (Backoff dB below full-scale int16 energy).
    double TargetPower = pow(2.0,2*NumBits-1)*pow(10.0,-0.1*Backoff);
    double TotalPower = TxPower + N0;
    double kAll = sqrt(TargetPower/TotalPower);
    double kSig = kAll/kTx;
    // Per real dimension variance N0/2, scaled same as signal so Es/N0 matches the chosen N0.
    Stdn = sqrt(N0*0.5)*kAll;
    mkSig = _mm256_set1_ps(kSig);
    // Piecewise schedule in *batches*: stable1 | acc1 | stable2 | acc2 | stop (see GenerateOutput segment switch).
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

    displayPeriodSec_ = (pParams) ? pParams->DisplayPeriodSec : 1.0;

    const double freqShift0Hz = pParams->InitialFrequencyShift;
    const double freqShiftDeltaHz = pParams->FrequencyShift;
    const double freqShiftStartHz = freqShift0Hz;
    const double freqShiftPeakHz = freqShift0Hz + freqShiftDeltaHz;

    const double rangeDb0 = pParams->InitialGainDb;
    const double rangeDbDelta = pParams->DynamicRangeDb;
    const double rangeDbStart = rangeDb0;
    const double rangeDbLow = rangeDb0 + rangeDbDelta;
    const double rangeDbMin = std::min(rangeDbStart, rangeDbLow);
    const double rangeDbMax = std::max(rangeDbStart, rangeDbLow);
    const double ppmDoppler0 =
        freqShiftStartHz * 1.0e6 / (pParams->CarrierToSymbolRateRatio * SymbolRate);
    const double ppmDopplerPeak =
        freqShiftPeakHz * 1.0e6 / (pParams->CarrierToSymbolRateRatio * SymbolRate);
    const double ppmDopplerMin = std::min(ppmDoppler0, ppmDopplerPeak);
    const double ppmDopplerMax = std::max(ppmDoppler0, ppmDopplerPeak);
    const double totalPpm0 = pParams->ClockMismatchPpm + ppmDoppler0;
    const double totalPpmPeak = pParams->ClockMismatchPpm + ppmDopplerPeak;
    const double totalPpmMin = std::min(totalPpm0, totalPpmPeak);
    const double totalPpmMax = std::max(totalPpm0, totalPpmPeak);
    const double eps = totalPpm0 * 1.0e-6;
    const double actualSymbolRate = SymbolRate * (1.0 + eps);

    if (displayPeriodSec_ > 0.0) {
        std::cout << std::fixed << std::setprecision(6) << "[AWGN] Channel parameters\n"
                  << "        Es/N0 (dB):                    " << pParams->EsN0 << '\n'
                  << "        InitialGainDb (dB):          " << pParams->InitialGainDb << '\n'
                  << "        DynamicRangeDb (dB):         " << pParams->DynamicRangeDb << '\n'
                  << "        RangeDb range (dB):          [" << rangeDbMin << " .. " << rangeDbMax << "]\n"
                  << "        InitialFrequencyShift (Hz):  " << pParams->InitialFrequencyShift << '\n'
                  << "        FrequencyShift (Hz):         " << pParams->FrequencyShift << '\n'
                  << "        FrequencyShift range (Hz):   [" << freqShiftStartHz << " .. " << freqShiftPeakHz << "]\n"
                  << "        ClockMismatchPpm:            " << pParams->ClockMismatchPpm << '\n'
                  << "        CarrierToSymbolRateRatio:    " << pParams->CarrierToSymbolRateRatio
                  << " (fc/Rs)\n"
                  << "        Symbol rate Rs (sym/s):      " << SymbolRate << '\n'
                  << "        Actual Rs from ppm (sym/s):  " << actualSymbolRate << '\n'
                  << "        Stable1 (s):                 " << StablePeriod << '\n'
                  << "        Acceleration1 (s):           " << AccelerationPeriod << '\n'
                  << "        Stable2 (s):                 " << StablePeriod << '\n'
                  << "        Acceleration2 (s):           " << AccelerationPeriod << '\n'
                  << "        TotalPeriod (s):             " << TotalPeriod << '\n'
                  << "        DisplayPeriodSec (s):        " << displayPeriodSec_
                  << (displayPeriodSec_ > 0.0 ? "" : " (periodic status off)") << '\n'
                  << "        SamplingClockOffset total ppm (t=0): " << totalPpm0 << '\n'
                  << "        SamplingClockOffset total ppm range: [" << totalPpmMin << " .. " << totalPpmMax << "]\n"
                  << "          clock mismatch component:  " << pParams->ClockMismatchPpm << '\n'
                  << "          Doppler component (t=0):   " << ppmDoppler0 << '\n'
                  << "          Doppler component range:   [" << ppmDopplerMin << " .. " << ppmDopplerMax << "]\n"
                  << "          (Doppler from FrequencyShift / (fc/Rs) / Rs)\n";
    }

    samplingClockOffset_.Configure(totalPpm0);
    samplingClockOffset_.Start(&OutputBuffer, &mtxOutputBuffer_, &CvUserOut, &CvOutUser, &StopAll);

    // NCO initial Hz (per-chunk updates happen in ApplyFrequencyOffset from FreqChunk::freqHz).
    objFreqOffset.SetSamplingFrequency(SamplingFrequency);
    objFreqOffset.SetFrequency(static_cast<Ipp64f>(freqShiftStartHz));

    NoiseThread = std::thread(&AWGNChannel::GenerateNoise, this);
    OutputThread = std::thread(&AWGNChannel::GenerateOutput, this);
    FreqOffsetThread = std::thread(&AWGNChannel::ApplyFrequencyOffset, this);
}



void AWGNChannel::StopThreads(void)
{
    StopAll = true;
    CvOutNoise.notify_all();
    CvOutUser.notify_all();
    CvNoiseOut.notify_all();
    CvUserOut.notify_all();
    cvFreqQData_.notify_all();
    cvFreqQSpace_.notify_all();
    if(NoiseThread.joinable())
        NoiseThread.join();
    if(OutputThread.joinable())
        OutputThread.join();
    if(FreqOffsetThread.joinable())
        FreqOffsetThread.join();
    samplingClockOffset_.StopJoin();
}


/**
 * @brief Producer thread: Gaussian noise buffers into @c NoiseQ for @ref GenerateOutput.
 */
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
        std::unique_lock<std::mutex> lk(mtxNoiseQ_);
        CvOutNoise.wait(lk, [&] { return StopAll || NoiseQ.AvailableWrite(); });
        if (StopAll)
            break;
        const unsigned int PtrWr = NoiseQ.GetPtrWr();
        oNoiseGen.randn(Noise[PtrWr], NoiseBatchSize);
        NoiseQ.AdvanceWrite();
        lk.unlock();
        CvNoiseOut.notify_one();
    }
}


/**
 * @brief Channel output thread: TX + AWGN mix, dynamic range, enqueue chunks for @ref ApplyFrequencyOffset.
 */
void AWGNChannel::GenerateOutput(void)
{
    #ifdef WRITE_LOG_THR
	mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"Channel Output Thread %d\n", gettid());
    fclose(fidthr);
    mtxfilethr.unlock();
    #endif
    int batchCount = 0;
    const auto wall_start = std::chrono::steady_clock::now();
    auto lastStatusDisplay = wall_start;
   
    while(!StopAll)
    {
        unsigned int PtrRdNoise = 0;
        {
            std::unique_lock<std::mutex> lkNoise(mtxNoiseQ_);
            CvNoiseOut.wait(lkNoise, [&] { return StopAll || NoiseQ.AvailableRead(); });
            if (StopAll)
                break;
            PtrRdNoise = NoiseQ.GetPtrRd();
        }

        alignas(32) float noisyChunk[TxOutputBatchSize];
        if (!pTx->CopyOutputSamples(noisyChunk, TxOutputBatchSize, StopAll))
            break;

        // Calculate current t_rate and rangeDb BEFORE adding noise so we can apply gain before noise if requested.
        const double t_rate = static_cast<double>(batchCount) * TxOutputBatchDuration;
        int seg = 0;
        double deltaHz = 0.0;
        double rangeDb = pParams->InitialGainDb;
        const double a = AccelerationPeriod;
        const double s = StablePeriod;
        if (t_rate < s)
        {
            seg = 0; // stable1
            deltaHz = 0.0;
            rangeDb = pParams->InitialGainDb;
        }
        else if (t_rate < s + a)
        {
            seg = 1; // acc1: 0 -> +FrequencyShift
            const double x = (t_rate - s) / std::max(1e-12, a);
            deltaHz = std::max(0.0, std::min(1.0, x)) * pParams->FrequencyShift;
            rangeDb = pParams->InitialGainDb +
                      std::max(0.0, std::min(1.0, x)) * pParams->DynamicRangeDb;
        }
        else if (t_rate < s + a + s)
        {
            seg = 2; // stable2
            deltaHz = pParams->FrequencyShift;
            rangeDb = pParams->InitialGainDb + pParams->DynamicRangeDb;
        }
        else
        {
            seg = 3; // acc2: +FrequencyShift -> 0
            const double x = (t_rate - (s + a + s)) / std::max(1e-12, a);
            deltaHz = (1.0 - std::max(0.0, std::min(1.0, x))) * pParams->FrequencyShift;
            rangeDb = pParams->InitialGainDb +
                      (1.0 - std::max(0.0, std::min(1.0, x))) * pParams->DynamicRangeDb;
        }

        const double freqHz = pParams->InitialFrequencyShift + deltaHz;
        const double ppmDoppler = freqHz * 1.0e6 / (pParams->CarrierToSymbolRateRatio * SymbolRate);
        const double totalPpm = pParams->ClockMismatchPpm + ppmDoppler;

        const float gain = static_cast<float>(std::pow(10.0, rangeDb / 20.0));

        // noisyChunk <- mkSig * tx + Stdn * noise (interleaved float IQ, one float per I or Q sample).
        {
            __m256 mStdn = _mm256_set1_ps(Stdn);
            __m256 mg = _mm256_set1_ps(gain);
            __m256 mSigScale = mkSig;
            if (pParams->ApplyGainBeforeNoise) {
                mSigScale = _mm256_mul_ps(mSigScale, mg);
            }

            for (int i = 0; i < TxOutputBatchSize;)
            {
                __m256 mInl = _mm256_loadu_ps(noisyChunk + i);
                __m256 mOutl = _mm256_loadu_ps(Noise[PtrRdNoise] + i);
                i += 8;
                __m256 mInh = _mm256_loadu_ps(noisyChunk + i);
                __m256 mOuth = _mm256_loadu_ps(Noise[PtrRdNoise] + i);
                i += 8;
                mOutl = _mm256_mul_ps(mStdn, mOutl);
                mOuth = _mm256_mul_ps(mStdn, mOuth);

                mOutl = _mm256_fmadd_ps(mInl, mSigScale, mOutl);
                mOuth = _mm256_fmadd_ps(mInh, mSigScale, mOuth);

                if (!pParams->ApplyGainBeforeNoise) {
                    mOutl = _mm256_mul_ps(mOutl, mg);
                    mOuth = _mm256_mul_ps(mOuth, mg);
                }

                _mm256_storeu_ps(noisyChunk + i - 16, mOutl);
                _mm256_storeu_ps(noisyChunk + i - 8, mOuth);
            }
            batchCount++;
        }

#ifdef DEBUG_AWGN
        {
            alignas(32) short dbgOut[TxOutputBatchSize];
            unsigned int PtrOut = 0;
            for (int i = 0; i < TxOutputBatchSize;)
            {
                __m256 mInl = _mm256_loadu_ps(noisyChunk + i);
                i += 8;
                __m256 mInh = _mm256_loadu_ps(noisyChunk + i);
                i += 8;
                __m256i mOuti = floats_to_shorts_sat_perm_avx2(mInl, mInh);
                _mm256_storeu_si256((__m256i*)(dbgOut + PtrOut), mOuti);
                PtrOut += 16;
            }
            std::copy(dbgOut, dbgOut + TxOutputBatchSize, OutAllI + PtrOutAll);
            PtrOutAll += TxOutputBatchSize;
            cout << "Collected " << PtrOutAll << endl;
            if (PtrOutAll >= 1100000)
            {
                FILE* fid = fopen("AWGNOut.bin", "wb");
                fwrite(OutAllI, sizeof(short), PtrOutAll, fid);
                fclose(fid);
                cout << "Saved AWGN " << PtrOutAll << endl;
                std::this_thread::sleep_for(10ms);

                Finish.store(true, std::memory_order_relaxed);
            }
        }
#endif

        // Decouple heavy NCO + SCO work from TX+noise so OutputThread keeps pulling transmitter samples.
        {
            std::unique_lock<std::mutex> lk(mtxFreqQ_);
            cvFreqQSpace_.wait(lk, [&] { return StopAll || static_cast<int>(freqQ_.size()) < kFreqQDepth; });
            if (StopAll)
                break;

            if (batchCount >= static_cast<int>(TransitionCounter[4]))
            {
                Finish.store(true, std::memory_order_relaxed);
                StopAll = true;
                break;
            }

            FreqChunk chunk;
            chunk.interleaved.resize(static_cast<size_t>(TxOutputBatchSize));
            std::copy(noisyChunk, noisyChunk + TxOutputBatchSize, chunk.interleaved.begin());
            chunk.freqHz = freqHz;
            chunk.totalPpm = totalPpm;
            chunk.gainDb = rangeDb;
            chunk.segment = seg;
            chunk.t_rate = t_rate;
            freqQ_.push_back(std::move(chunk));
        }
        cvFreqQData_.notify_one();

        {
            std::unique_lock<std::mutex> lkNoise(mtxNoiseQ_);
            NoiseQ.AdvanceRead();
        }
        CvOutNoise.notify_one();

        // Periodic channel status line.
        {
            const auto now = std::chrono::steady_clock::now();
        if (displayPeriodSec_ > 0.0 &&
            now - lastStatusDisplay >= std::chrono::duration<double>(displayPeriodSec_))
        {
            // const char* segName = "stable1";
            // // Best-effort read of last enqueued segment (no extra locking: we reuse the already pushed logic).
            // // For clarity, we just derive it again from batchCount/t_rate.
            // const double t_rate = static_cast<double>(batchCount) * TxOutputBatchDuration;
            // const double a = AccelerationPeriod;
            // const double s = StablePeriod;
            // double deltaHz = 0.0;
            // if (t_rate < s) { segName = "stable1"; deltaHz = 0.0; }
            // else if (t_rate < s + a) { segName = "acc1"; deltaHz = (t_rate - s) / std::max(1e-12, a) * pParams->FrequencyShift; }
            // else if (t_rate < s + a + s) { segName = "stable2"; deltaHz = pParams->FrequencyShift; }
            // else { segName = "acc2"; deltaHz = (1.0 - (t_rate - (s + a + s)) / std::max(1e-12, a)) * pParams->FrequencyShift; }
            // deltaHz = std::max(0.0, std::min(pParams->FrequencyShift, deltaHz));
            // const double freqHz = pParams->InitialFrequencyShift + deltaHz;
            // const double ppmDoppler =
            //     freqHz * 1.0e6 / (pParams->CarrierToSymbolRateRatio * SymbolRate);
            // const double totalPpm = pParams->ClockMismatchPpm + ppmDoppler;
            // double gainDb = pParams->InitialGainDb;
            // if (t_rate < s) { gainDb = pParams->InitialGainDb; }
            // else if (t_rate < s + a) { gainDb = pParams->InitialGainDb + (t_rate - s) / std::max(1e-12, a) * pParams->DynamicRangeDb; }
            // else if (t_rate < s + a + s) { gainDb = pParams->InitialGainDb + pParams->DynamicRangeDb; }
            // else { gainDb = pParams->InitialGainDb + (1.0 - (t_rate - (s + a + s)) / std::max(1e-12, a)) * pParams->DynamicRangeDb; }
            // const double t_sim = std::chrono::duration<double>(now - wall_start).count();
            // double snrDb = pParams->ApplyGainBeforeNoise ? (pParams->EsN0 + gainDb) : pParams->EsN0;
            // std::cout << std::fixed << std::setprecision(6)
            //           << "[AWGN] status seg=" << segName
            //           << " freqShiftHz=" << freqHz
            //           << " gainDb=" << gainDb
            //           << " snrDb=" << snrDb
            //           << " ppm=" << totalPpm
            //           << " t_sim=" << t_sim << " s"
            //           << " t_rate=" << t_rate << " s"
            //           << std::endl;
            lastStatusDisplay = now;
        }
        }
    }
}

/**
 * @brief Applies per-chunk carrier rotation and forwards floats to @ref ChannelSamplingClockOffset.
 */
void AWGNChannel::ApplyFrequencyOffset(void)
{
    const unsigned int nComplex = static_cast<unsigned int>(TxOutputBatchSize / 2);
    std::vector<float> iBuf(nComplex);
    std::vector<float> qBuf(nComplex);
    FreqChunk chunk;

    while (!StopAll)
    {
        {
            std::unique_lock<std::mutex> lk(mtxFreqQ_);
            cvFreqQData_.wait(lk, [&] { return StopAll || !freqQ_.empty(); });
            if (StopAll)
                break;
            chunk = std::move(freqQ_.front());
            freqQ_.pop_front();
            lk.unlock();
            cvFreqQSpace_.notify_one();
        }

        if (chunk.interleaved.size() != static_cast<size_t>(TxOutputBatchSize))
            continue;

        // Update time-varying impairments for this chunk.
        objFreqOffset.SetFrequency(static_cast<Ipp64f>(chunk.freqHz));
        samplingClockOffset_.UpdateTotalPpm(chunk.totalPpm);

        // Publish current applied values (for GUI/status overlay).
        currFreqHz_.store(chunk.freqHz, std::memory_order_relaxed);
        currTotalPpm_.store(chunk.totalPpm, std::memory_order_relaxed);
        currGainDb_.store(chunk.gainDb, std::memory_order_relaxed);
        currSegment_.store(chunk.segment, std::memory_order_relaxed);
        currTRate_.store(chunk.t_rate, std::memory_order_relaxed);

        for (unsigned int k = 0; k < nComplex; ++k)
        {
            iBuf[k] = chunk.interleaved[2 * k];
            qBuf[k] = chunk.interleaved[2 * k + 1];
        }
        objFreqOffset.CreateOutputs(iBuf.data(), qBuf.data(), nComplex);
        for (unsigned int k = 0; k < nComplex; ++k)
        {
            chunk.interleaved[2 * k] = iBuf[k];
            chunk.interleaved[2 * k + 1] = qBuf[k];
        }

        samplingClockOffset_.EnqueueNoisyInterleaved(chunk.interleaved.data(), TxOutputBatchSize);
    }
}

/**
 * @brief Blocks until @c nShorts IQ shorts are available from @ref OutputBuffer (post-SCO).
 */
bool AWGNChannel::CopyOutputSamples(short* dst, int nShorts, std::atomic<bool>& stopAll)
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
