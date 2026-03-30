#pragma once

#include "BufferShort.h"
#include "definitions.h"
#include <condition_variable>
#include <deque>
#include <mutex>
#include <thread>
#include <vector>

/// Resamples the 2×Fs complex stream to emulate a sampling-clock offset on the receive path:
/// output at effective Fs·(1 + ppm·10⁻⁶) relative to the input.
/// +ppm means a faster Rx clock (same convention as the former Gardner-input SCO).
class ChannelSamplingClockOffset
{
public:
    ChannelSamplingClockOffset();
    ~ChannelSamplingClockOffset();

    void Configure(double totalPpm);

    /// Output buffer shared with CopyOutputSamples (same mutex and condition variables as the channel).
    void Start(BufferShort* outputBuffer, std::mutex* outputMutex,
               std::condition_variable* cvSpace, std::condition_variable* cvData,
               bool* stopAll);

    void StopJoin();

    /// Called from the noise/output thread: enqueues one batch of TxOutputBatchSize interleaved IQ floats.
    void EnqueueNoisyInterleaved(const float* interleaved, int nFloats);

    void ThreadMain();

private:
    static constexpr int kInputQueueDepth = 4;
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
    };

    double totalPpm_ = 0.0;
    double step_ = 1.0;

    BufferShort* pOutBuf_ = nullptr;
    std::mutex* pOutMtx_ = nullptr;
    std::condition_variable* pCvSpace_ = nullptr;
    std::condition_variable* pCvData_ = nullptr;
    bool* pStopAll_ = nullptr;

    std::mutex mtxIn_;
    std::condition_variable cvInSpace_;
    std::condition_variable cvInData_;
    std::deque<std::vector<float>> inQueue_;

    std::thread thread_;
    bool threadRunning_ = false;

    static void FloatBatchToShortsInterleaved(const float* iPtr, const float* qPtr, int nComplex,
                                              short* dst);
};
