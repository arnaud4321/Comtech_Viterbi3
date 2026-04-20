/**
 * @file ChannelSamplingClockOffset.h
 * @brief Optional Rx sampling-clock error: resample Tx/channel stream by @f$(1+\mathrm{ppm}\cdot10^{-6})@f$.
 */
#pragma once

#include "BufferShort.h"
#include "definitions.h"
#include <condition_variable>
#include <deque>
#include <mutex>
#include <atomic>
#include <thread>
#include <vector>

/**
 * @brief Emulates receiver clock error by asynchronous resampling of the complex baseband stream.
 *
 * @details **Input path:** @ref EnqueueNoisyInterleaved appends one @c TxOutputBatchSize interleaved-float
 * batch to a bounded queue (@c kInputQueueDepth); the worker dequeues, de-interleaves to I/Q, and
 * @c PushBlock into a large ring (@c ScoRingBuffer, power-of-two length).
 *
 * **Time step:** @c UpdateTotalPpm sets @f$\varepsilon=\mathrm{ppm}\cdot10^{-6}@f$ and stores
 * @f$\texttt{step}=1/(1+\varepsilon)@f$ in an atomic (read each batch in @c ThreadMain). The fractional
 * read cursor @c tScoAbs advances by @c step per **output** sample; average output rate is
 * @f$F_s(1+\varepsilon)@f$ vs input @f$F_s@f$ (positive ppm ⇒ more output samples per wall-clock time).
 *
 * **Interpolation / startup:** After enough input has accumulated, @c tScoAbs is initialized with a fixed
 * lag (@c latencySamples = @c kComplexPerBatch·64) behind the newest ring index so Lagrange-4 has margin.
 * Each output uses Lagrange-4 at @c tScoAbs (@ref Lagrange4Simd, AVX in the hot loop); if the cursor leaves a safe band inside the ring,
 * the lag state is dropped until warmup completes again.
 *
 * **Long runs (weeks / months):** @c absWrite and @c tScoAbs are periodically reduced by one full ring length
 * (@c kSize, power of two) so @c double fractional precision for Lagrange @c mu stays adequate and ring
 * indices stay well inside @c int range for masking — same idea as @c GardnerTiming wrap.
 *
 * **Output:** @c FloatBatchToShortsInterleaved applies @c lrintf and clamps to int16; results go to the shared
 * @c BufferShort like the rest of the channel. Distinct from Gardner recovery, which acts on **symbol**
 * timing after matched filtering in @ref Receiver.
 */
class ChannelSamplingClockOffset
{
public:
    ChannelSamplingClockOffset();
    ~ChannelSamplingClockOffset();

    void Configure(double totalPpm);
    void UpdateTotalPpm(double totalPpm);

    /// Output buffer shared with CopyOutputSamples (same mutex and condition variables as the channel).
    void Start(BufferShort* outputBuffer, std::mutex* outputMutex,
               std::condition_variable* cvSpace, std::condition_variable* cvData,
               std::atomic<bool>* stopAll);

    void StopJoin();

    /// Called from the noise/output thread: enqueues one batch of TxOutputBatchSize interleaved IQ floats.
    void EnqueueNoisyInterleaved(const float* interleaved, int nFloats);

    void ThreadMain();

private:
    static constexpr int kInputQueueDepth = kScoInputQueueDepth;
    static constexpr int kComplexPerBatch = TxOutputBatchSize / 2;
    static constexpr int kScoOutMax = kComplexPerBatch + 512;

    struct ScoRingBuffer
    {
        static constexpr int kSize = 1 << 19;
        static constexpr int kMask = kSize - 1;

        alignas(32) float ringI[kSize] = {};
        alignas(32) float ringQ[kSize] = {};
        long long absWrite = 0;

        void PushBlock(const float* inI, const float* inQ, int n);
        long long OldestAbs() const;
        float GetI(long long absIdx) const;
        float GetQ(long long absIdx) const;
        float InterpLagrange4I(double t) const;
        float InterpLagrange4Q(double t) const;

        /// Periodically subtract @c kSize from @c absWrite and, when @a tCursor is non-null, from @a *tCursor
        /// so absolute sample time stays bounded during multi-week / multi-month runs (double ULP, int mask).
        void RewrapForLongRun(double* tCursor);
    };

    double totalPpm_ = 0.0;
    std::atomic<double> step_{1.0};

    BufferShort* pOutBuf_ = nullptr;
    std::mutex* pOutMtx_ = nullptr;
    std::condition_variable* pCvSpace_ = nullptr;
    std::condition_variable* pCvData_ = nullptr;
    std::atomic<bool>* pStopAll_ = nullptr;

    std::mutex mtxIn_;
    std::condition_variable cvInSpace_;
    std::condition_variable cvInData_;
    std::deque<std::vector<float>> inQueue_;

    std::thread thread_;
    bool threadRunning_ = false;

    static void FloatBatchToShortsInterleaved(const float* iPtr, const float* qPtr, int nComplex,
                                              short* dst);
};
