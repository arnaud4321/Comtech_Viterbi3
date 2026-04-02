/**
 * @file definitions.h
 * @brief Global simulation rates, batch sizes, and transmitter/receiver mode tags.
 */
#pragma once

/** @brief Nominal symbol rate (sym/s). */
const double SymbolRate = 10.71e6;
/** @brief Complex baseband sampling rate at 2 samples per symbol (Hz). */
const double SamplingRate = 2.0*SymbolRate;
const int BatchSize1 = 1024;
const int BatchSize3 = 3 * BatchSize1;
const int BatchSize1Bytes = BatchSize1 >> 3;
const int BatchConvSize1 = BatchSize1 * 2;
const int TxOutputBatchSize = BatchSize3 * 4;
// TxOutputBatchSize is the number of float values in interleaved IQ (I,Q,I,Q,...),
// so the number of complex samples per batch is TxOutputBatchSize/2.
const double TxOutputBatchDuration = (static_cast<double>(TxOutputBatchSize) / 2.0) / SamplingRate;
const int ReceiverInputBatchIQSamples = BatchSize3 * 2; //BatchSize3 x2 (conv) /2 (IQ) x2 (2x)
const int ReceiverInputBatchIQSymbols = BatchSize3;
const int SPB = 8192;
const double SamplingFrequency = 21.42e6;

/** @brief Transmitter source: internal PRBS or IQ file. */
enum TxModes
{ PRBS_TX,FILE_TX};

/** @brief High-level simulation mode (acquisition vs continuous). */
enum SimModes
{
    ACQ_SIM,CONT_SIM
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
