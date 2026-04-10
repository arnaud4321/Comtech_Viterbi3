/**
 * @file TransmitterManager.h
 * @brief USRP hardware transmitter integration using the UHD API.
 */
#pragma once
#include "Transmitter.h"
#include "random_generator_new.h"
#include <thread>
#include <uhd/usrp/multi_usrp.hpp>
#include <condition_variable> // std::condition_variable
#include <mutex>
using namespace std;
class AWGNChannel;

/**
 * @brief Manages the continuous transmission of IQ samples to a USRP device.
 * 
 * @details This class runs a dedicated thread that constantly polls the @ref Transmitter
 * for new batches of floating-point IQ samples via `CopyOutputSamples`. It then scales
 * the samples into 16-bit integers and feeds them into the USRP using the `uhd::tx_streamer` API.
 */
class TransmitterManager
{
    uhd::usrp::multi_usrp::sptr usrp;
    std::atomic<bool> StopAll{false};
    Transmitter *pTx;
    uhd::tx_streamer::sptr tx_stream = 0;
    void OperateThread(void);
    std::thread TxThread;
    public:
    uint64_t NumBatches=0;
    TransmitterManager(uhd::usrp::multi_usrp::sptr usrpx);
    ~TransmitterManager();
    void StartThreads();
    void StopThreads();
    void SetTransmitter(Transmitter *p)
    {
        pTx = p;
    }
};