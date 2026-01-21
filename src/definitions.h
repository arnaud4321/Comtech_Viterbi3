#pragma once
const int BatchSize1 = 1024;
const int BatchSize3 = 3 * BatchSize1;
const int BatchSize1Bytes = BatchSize1 >> 3;
const int BatchConvSize1 = BatchSize1 * 2;
const int TxOutputBatchSize = BatchSize3 * 4;
const int ReceiverInputBatchIQSamples = BatchSize3 * 2; //BatchSize3 x2 (conv) /2 (IQ) x2 (2x)
const int ReceiverInputBatchIQSymbols = BatchSize3;
const int SPB = 8192;
const double SamplingFrequency = 21.42e6;

enum TxModes
{ PRBS_TX,FILE_TX};

enum SimModes
{
    ACQ_SIM,CONT_SIM
};
const float ViterbiThreshold1 = 150; //difference between best metric and other 3 options
const int PRBSThreshold = 12;

#define WRITE_LOG_THR
