/**
 * @file definitions.h
 * @brief Global simulation rates, batch sizes, and transmitter/receiver mode tags.
 */
#pragma once

/** @brief Nominal symbol rate (sym/s). */
inline double SymbolRate = 10.71e6;
/** @brief Complex baseband sampling rate at 2 samples per symbol (Hz). */
inline double SamplingRate = 2.0*SymbolRate;
const int BatchSize1 = 1024;
const int BatchSize3 = 3 * BatchSize1;
const int BatchSize1Bytes = BatchSize1 >> 3;
const int BatchConvSize1 = BatchSize1 * 2;
const int TxOutputBatchSize = BatchSize3 * 4;
// TxOutputBatchSize is the number of float values in interleaved IQ (I,Q,I,Q,...),
// so the number of complex samples per batch is TxOutputBatchSize/2.
inline double TxOutputBatchDuration = (static_cast<double>(TxOutputBatchSize) / 2.0) / SamplingRate;
const int ReceiverInputBatchIQSamples = BatchSize3 * 2; //BatchSize3 x2 (conv) /2 (IQ) x2 (2x)
const int ReceiverInputBatchIQSymbols = BatchSize3;
const int SPB = 8192;
inline double SamplingFrequency = 21.42e6;

/** @brief Float ring length (per I/Q) for RX @c BufferFloat stages (filter, resampler, freq corrector). Larger ⇒ more headroom under backpressure. */
inline constexpr int kRxRingFloatLen = SPB * 1024;
/** @brief Extra margin for wrapped access in RX @c BufferFloat rings. */
inline constexpr int kRxRingFloatExtra = SPB * 64;
/** @brief Interleaved-short ring length for AWGN channel output and @ref Sampler (shorts). */
inline constexpr int kSamplerShortRingSize = 1024 * SPB;
inline constexpr int kSamplerShortRingExtra = 4 * SPB;
/** @brief @c BufferShort capacity (shorts) on AWGN output path (ties to @ref TxOutputBatchSize). */
inline constexpr int kChannelOutShortRingSize = 256 * TxOutputBatchSize;
inline constexpr int kChannelOutShortRingExtra = 4 * TxOutputBatchSize;
/** @brief Resampler internal FIFO cap (floats per I/Q); bounds drain from matched-filter ring. */
inline constexpr int kResamplerFifoMaxSamples = SPB * 1024;
/** @brief SCO input queue depth (float batches) before @ref ChannelSamplingClockOffset::EnqueueNoisyInterleaved blocks. */
inline constexpr int kScoInputQueueDepth = 36;

/** @brief Transmitter source: internal PRBS or IQ file. */
enum TxModes
{ PRBS_TX,FILE_TX};

/** @brief High-level simulation mode (acquisition vs continuous). */
enum SimModes
{
    NOT_SIM,ACQ_SIM,CONT_SIM
};
enum OpModes
{
    NOT_OP,TX_ONLY,RX_ONLY,TX_RX
};
const float ViterbiThreshold1 = 80; //difference between best metric and other 3 options
const int PRBSThreshold = 12;

/**
 * @brief PRBS test injection: deterministic bit flips applied before PRBS sync search.
 *
 * When @c PRBSInjectStride > 0, bit @c in[i] is flipped for indices where
 * (i - PRBSInjectStart) % PRBSInjectStride == 0, for all i in [PRBSInjectStart, Length).
 * Use this to verify that @c PRBSSynchronized remains false when the descrambled PRBS
 * stream is corrupted.
 *
 * Set @c PRBSInjectStride = 0 to disable injection.
 */
const int PRBSInjectStride = -1;   ///< 0 disables; otherwise flips ~1/stride bits.
const int PRBSInjectStart = 0;    ///< Start index (in bits / bytes) for injection.

#define WRITE_LOG_THR
#define RX_INITIAL_GAIN 30
#define AGCBACKOFF 18
