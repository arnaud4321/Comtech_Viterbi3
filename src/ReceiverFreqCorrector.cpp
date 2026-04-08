/**
 * @file ReceiverFreqCorrector.cpp
 * @brief Single thread: pull 2× oversampled IQ, optional VIVA blocks, AGC, NCO, push to Gardner input ring.
 *
 * @details **ThreadMain** — Reads @c SPB-sized chunks from @c inRing_. Accumulates samples for
 * @ref CentralFreqEstimatorViva::Run when enabled; on detection smooths Hz with EMA and slew limit into
 * @ref FrequencyOffset. **AGC:** chunk mean power → EMA → linear gain @f$\sqrt{P_{\mathrm{ref}}/P_{\mathrm{ema}}}@f$.
 * Re-enables estimation on @ref ReceiverPhaseTrackingDD delock. Forwards to @c outRing_ only when
 * @c freqReady_ (or estimator disabled). Drops chunks while waiting for first lock if gating is on.
 */

#include "ReceiverFreqCorrector.h"
#include "ConsoleAlert.h"

#include <algorithm>
#include <cmath>
#include <chrono>
#include <iostream>
extern mutex mtxfilethr;

ReceiverFreqCorrector::ReceiverFreqCorrector()
{
    nco_.SetSamplingFrequency(SamplingFrequency);
}

ReceiverFreqCorrector::~ReceiverFreqCorrector()
{
    StopJoin();
}

void ReceiverFreqCorrector::Configure(const ReceiverFreqCorrectorConfig& cfg, double displayPeriodSec)
{
    cfg_ = cfg;
    if (cfg_.EstimationBlockSamples < 4096)
        cfg_.EstimationBlockSamples = 4096;
    if (cfg_.EstimatePeriodSec <= 0.0)
        cfg_.EstimatePeriodSec = 1.0;
    if (cfg_.NcoHzEmaAlpha <= 0.0 || cfg_.NcoHzEmaAlpha > 1.0)
        cfg_.NcoHzEmaAlpha = 0.2;
    if (cfg_.PowerEmaAlpha <= 0.0 || cfg_.PowerEmaAlpha > 1.0)
        cfg_.PowerEmaAlpha = 0.01;
    displayPeriodSec_ = displayPeriodSec;

    accI_.assign(static_cast<size_t>(cfg_.EstimationBlockSamples), 0.0f);
    accQ_.assign(static_cast<size_t>(cfg_.EstimationBlockSamples), 0.0f);
    accCount_ = 0;
    allowEstimation_ = true;
}

void ReceiverFreqCorrector::Start(BufferFloat* inRing, std::mutex* inMtx, std::condition_variable* inCvData,
                                  BufferFloat* outRing, std::mutex* outMtx, std::condition_variable* outCvData,
                                  ReceiverPhaseTrackingDD* phaseDD, std::atomic<bool>* stopAll)
{
    StopJoin();
    inRing_ = inRing;
    inMtx_ = inMtx;
    inCvData_ = inCvData;
    outRing_ = outRing;
    outMtx_ = outMtx;
    outCvData_ = outCvData;
    phaseDD_ = phaseDD;
    stopAll_ = stopAll;

    pwrEma_ = 0.0;
    pwrRef_ = 0.0;
    gainDb_.store(0.0, std::memory_order_relaxed);
    prevPhaseLocked_ = false;
    allowEstimation_ = true;
    freqReady_.store(!cfg_.Enable, std::memory_order_relaxed);
    accCount_ = 0;
    samplesTotal_ = 0;
    droppedChunks_ = 0;
    targetHz_.store(0.0, std::memory_order_relaxed);
    hasTargetHz_ = false;

    est_.Configure(SamplingFrequency, cfg_.Estimator);

    if (cfg_.TargetAvgPower > 0.0)
        pwrRef_ = cfg_.TargetAvgPower;

    running_ = true;
    thread_ = std::thread(&ReceiverFreqCorrector::ThreadMain, this);
}

void ReceiverFreqCorrector::StopJoin()
{
    running_ = false;
    if (thread_.joinable())
        thread_.join();
}

/**
 * @brief NCO + AGC + optional VIVA processing loop.
 */
void ReceiverFreqCorrector::ThreadMain()
{


    
    #ifdef WRITE_LOG_THR
	mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"Receiver Freq Corrector  Thread %d\n", gettid());
    fclose(fidthr);
    mtxfilethr.unlock();
    #endif


    constexpr int kChunk = SPB; // process in SPB-sized blocks (2 sps complex)
    std::vector<float> tmpI(static_cast<size_t>(kChunk));
    std::vector<float> tmpQ(static_cast<size_t>(kChunk));

    const auto wall_start = std::chrono::steady_clock::now();
    auto lastStatus = std::chrono::steady_clock::now();
    auto lastEstLog = lastStatus;
    auto nextEstimateStart = std::chrono::steady_clock::now();
    bool estimateCollecting = true;
    const double chunkDt = static_cast<double>(kChunk) / SamplingFrequency;
    auto lastOutRingSatLog = std::chrono::steady_clock::now();

    while (running_ && (stopAll_ == nullptr || !*stopAll_))
    {
        // Pull from input ring
        {
            std::unique_lock<std::mutex> lk(*inMtx_);
            inCvData_->wait(lk, [&] {
                return !running_ || (stopAll_ && *stopAll_) || inRing_->GetSizeInBuffer() >= kChunk;
            });
            if (!running_ || (stopAll_ && *stopAll_))
                break;
            float* inI = nullptr;
            float* inQ = nullptr;
            inRing_->GetReadBuffer(inI, inQ, kChunk);
            if (!inI)
                continue;
            std::copy(inI, inI + kChunk, tmpI.begin());
            std::copy(inQ, inQ + kChunk, tmpQ.begin());
            inRing_->AdvancePtrRd(kChunk);
        }
        inCvData_->notify_one();
        samplesTotal_ += static_cast<uint64_t>(kChunk);

        // Detect PhaseDD delock to re-enable estimation
        const bool phaseLocked = phaseDD_ ? phaseDD_->IsLocked() : false;
        if (prevPhaseLocked_ && !phaseLocked)
        {
            freqReady_.store(!cfg_.Enable, std::memory_order_relaxed);
            accCount_ = 0;
            estimateCollecting = true;
            nextEstimateStart = std::chrono::steady_clock::now();
            hasTargetHz_ = false;
            std::cout << "[CentralFreq] \033[31mDELOCK\033[0m"
                      << " PhaseDD -> re-enable estimation" << std::endl;
        }
        prevPhaseLocked_ = phaseLocked;

        // Periodic estimation (every EstimatePeriodSec). We always collect samples; we only apply
        // the update smoothly to avoid delocking PhaseDD.
        const auto now_tp = std::chrono::steady_clock::now();
        if (now_tp >= nextEstimateStart)
        {
            estimateCollecting = true;
            accCount_ = 0;
            nextEstimateStart += std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                std::chrono::duration<double>(cfg_.EstimatePeriodSec));
            while (nextEstimateStart <= now_tp)
            {
                nextEstimateStart += std::chrono::duration_cast<std::chrono::steady_clock::duration>(
                    std::chrono::duration<double>(cfg_.EstimatePeriodSec));
            }
        }

        // Collect enough samples for the next estimate window.
        if (cfg_.Enable && estimateCollecting)
        {
            const int canCopy = std::min(kChunk, cfg_.EstimationBlockSamples - accCount_);
            if (canCopy > 0)
            {
                std::copy(tmpI.begin(), tmpI.begin() + canCopy, accI_.begin() + accCount_);
                std::copy(tmpQ.begin(), tmpQ.begin() + canCopy, accQ_.begin() + accCount_);
                accCount_ += canCopy;
            }
            if (accCount_ >= cfg_.EstimationBlockSamples)
            {
                estimateCollecting = false;
                const CentralFreqVivaResult r = est_.Run(accI_.data(), accQ_.data(), accCount_);

                // Use VIVA power as a robust reference for normalization (track slowly) when not using an absolute target.
                if (cfg_.TargetAvgPower <= 0.0 && r.AvgPower > 1e-20)
                {
                    if (pwrRef_ <= 0.0)
                        pwrRef_ = r.AvgPower;
                    else
                        pwrRef_ = 0.01 * r.AvgPower + 0.99 * pwrRef_;
                }

                if (r.Detected)
                {
                    const bool wasReady = freqReady_.load(std::memory_order_relaxed);
                    const double targetHz = r.OffsetHz;
                    targetHz_.store(targetHz, std::memory_order_relaxed);
                    hasTargetHz_ = true;

                    // On first acquisition, jump directly to target to remove large offset quickly.
                    if (!wasReady)
                    {
                        ncoHz_.store(targetHz, std::memory_order_relaxed);
                        nco_.SetFrequency(static_cast<Ipp64f>(-targetHz));
                    }
                    freqReady_.store(true, std::memory_order_relaxed);

                    const auto now = std::chrono::steady_clock::now();
                    if (!wasReady || (displayPeriodSec_ > 0.0 &&
                                      now - lastEstLog >= std::chrono::duration<double>(displayPeriodSec_)))
                    {
                        const double t_rate = static_cast<double>(samplesTotal_) / SamplingFrequency;
                        const double t_sim = std::chrono::duration<double>(now - wall_start).count();
                        std::cout << "[CentralFreq] \033[32mDetected\033[0m offset=" << targetHz << " Hz"
                                  << " [window " << cfg_.EstimationBlockSamples << " samples]"
                                  << " t_sim=" << t_sim << " s"
                                  << " t_rate=" << t_rate << " s"
                                  << " peak/median=" << r.PeakToMedian
                                  << " (threshold " << cfg_.Estimator.PeakToMedianThreshold << ")"
                                  << " dF=" << r.FftResolutionHz << " Hz"
                                  << " ncoHz=" << ncoHz_.load(std::memory_order_relaxed) << " Hz"
                                  << " slew=" << cfg_.MaxHzSlewRate << " Hz/s"
                                  << " alpha=" << cfg_.NcoHzEmaAlpha
                                  << " gainDb=" << gainDb_.load(std::memory_order_relaxed)
                                  << " pwrEma=" << pwrEma_
                                  << " pwrRef=" << pwrRef_
                                  << " vivaAvgPwr=" << r.AvgPower
                                  << std::endl;
                        lastEstLog = now;
                    }
                }
                else
                {
                    const auto now = std::chrono::steady_clock::now();
                    if (displayPeriodSec_ > 0.0 &&
                        now - lastEstLog >= std::chrono::duration<double>(displayPeriodSec_))
                    {
                        std::cout << "[CentralFreq] \033[31mNot detected\033[0m peak/median=" << r.PeakToMedian
                                  << " (thresh=" << cfg_.Estimator.PeakToMedianThreshold << ")"
                                  << std::endl;
                        lastEstLog = now;
                    }
                }
                accCount_ = 0;
            }
        }

        // Power estimate (for gain normalization)
        double p = 0.0;
        for (int i = 0; i < kChunk; ++i)
            p += static_cast<double>(tmpI[i]) * tmpI[i] + static_cast<double>(tmpQ[i]) * tmpQ[i];
        p /= static_cast<double>(kChunk);
        if (pwrEma_ <= 0.0)
        {
            pwrEma_ = p;
            if (pwrRef_ <= 0.0)
                pwrRef_ = p;
        }
        else
        {
            pwrEma_ = cfg_.PowerEmaAlpha * p + (1.0 - cfg_.PowerEmaAlpha) * pwrEma_;
        }
        const double gain = (pwrEma_ > 1e-20 && pwrRef_ > 1e-20) ? std::sqrt(pwrRef_ / pwrEma_) : 1.0;
        gainDb_.store(20.0 * std::log10(std::max(1e-20, gain)), std::memory_order_relaxed);

        // Apply gain + NCO
        const float g = static_cast<float>(gain);
        for (int i = 0; i < kChunk; ++i)
        {
            tmpI[i] *= g;
            tmpQ[i] *= g;
        }
        if (freqReady_.load(std::memory_order_relaxed))
        {
            // Smooth NCO tracking continuously (per chunk) toward the latest detected target.
            if (hasTargetHz_ && cfg_.MaxHzSlewRate > 0.0)
            {
                const double prevHz = ncoHz_.load(std::memory_order_relaxed);
                const double maxStep = cfg_.MaxHzSlewRate * chunkDt;
                const double err = targetHz_.load(std::memory_order_relaxed) - prevHz;
                const double step = std::clamp(err, -maxStep, maxStep);
                const double steppedHz = prevHz + step;
                const double appliedHz = (1.0 - cfg_.NcoHzEmaAlpha) * prevHz + cfg_.NcoHzEmaAlpha * steppedHz;
                if (std::abs(appliedHz - prevHz) > 1e-9)
                {
                    ncoHz_.store(appliedHz, std::memory_order_relaxed);
                    nco_.SetFrequency(static_cast<Ipp64f>(-appliedHz));
                }
            }
            nco_.CreateOutputs(tmpI.data(), tmpQ.data(), static_cast<unsigned int>(kChunk));
        }

        // Push to output ring only once central frequency is ready (like symbol-rate gating).
        if (freqReady_.load(std::memory_order_relaxed))
        {
            std::unique_lock<std::mutex> lk(*outMtx_);
            while (running_ && (stopAll_ == nullptr || !*stopAll_) && outRing_->AlmostFull())
            {
                const auto nowSat = std::chrono::steady_clock::now();
                if (nowSat - lastOutRingSatLog >= std::chrono::seconds(1))
                {
                    lastOutRingSatLog = nowSat;
                    CONSOLE_ALERT_STMT(std::cout << ConsoleAlert::kRedOpen
                                                 << "[CentralFreq] output float ring almost full; fill="
                                                 << outRing_->GetSizeInBuffer() << "/" << outRing_->GetBufferSize() - 1
                                                 << ConsoleAlert::kReset << std::endl;);
                }
                outCvData_->wait_for(lk, std::chrono::milliseconds(1));
            }
            if (!running_ || (stopAll_ && *stopAll_))
                break;
            float* outI = nullptr;
            float* outQ = nullptr;
            outRing_->GetWriteBuffer(outI, outQ, kChunk);
            if (!outI)
                continue;
            std::copy(tmpI.begin(), tmpI.end(), outI);
            std::copy(tmpQ.begin(), tmpQ.end(), outQ);
            outRing_->AdvancePtrWr(kChunk);
            outCvData_->notify_one();
        }
        else
        {
            droppedChunks_++;
        }

        // Periodic status
        const auto now = std::chrono::steady_clock::now();
        if (displayPeriodSec_ > 0.0 &&
            now - lastStatus >= std::chrono::duration<double>(displayPeriodSec_))
        {
            std::cout << "[CentralFreq] status"
                      << " phaseLocked=" << (phaseLocked ? 1 : 0)
                      << " ready=" << (freqReady_.load(std::memory_order_relaxed) ? 1 : 0)
                      << " targetHz=" << (hasTargetHz_ ? targetHz_.load(std::memory_order_relaxed) : 0.0)
                      << " ncoHz=" << ncoHz_.load(std::memory_order_relaxed)
                      << " gainDb=" << gainDb_.load(std::memory_order_relaxed)
                      << " pwrEma=" << pwrEma_
                      << " pwrRef=" << pwrRef_
                      << " droppedChunks=" << droppedChunks_
                      << std::endl;
            lastStatus = now;
        }
    }
}

