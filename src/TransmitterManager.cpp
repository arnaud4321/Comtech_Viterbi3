#include "TransmitterManager.h"
#include "AWGNChannel.h"
#include <thread>
extern mutex mtxfilethr;

TransmitterManager::TransmitterManager(uhd::usrp::multi_usrp::sptr usrpx):usrp(usrpx)
{

}
TransmitterManager::~TransmitterManager()
{

}

void TransmitterManager::OperateThread()
{
    NumBatches = 0;
    mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"UHD Tx %d\n",gettid());
    fclose(fidthr);
    mtxfilethr.unlock();

    std::string cpu_format ="sc16";
    std::string wirefmt="sc16";
    uhd::stream_args_t stream_args(cpu_format, wirefmt);
    std::vector<size_t> channel_nums;
    channel_nums.push_back(boost::lexical_cast<size_t>("0"));
    stream_args.channels             = channel_nums;
    tx_stream = usrp->get_tx_stream(stream_args);
    uhd::tx_metadata_t md;
    md.start_of_burst = true;
    md.end_of_burst   = false;

    alignas(32) float floatsOut[SPB*2];
    alignas(32) short ChOut[SPB*2];

    bool first = true;

    while(!StopAll)
    {
        if(!pTx->CopyOutputSamples(floatsOut, SPB*2, StopAll))
            break;

        for(int i = 0; i < SPB*2; i++) {
            float v = floatsOut[i] * 8192.0f; // Scale appropriately, avoids clipping.
            if (v > 32767.0f) v = 32767.0f;
            if (v < -32768.0f) v = -32768.0f;
            ChOut[i] = (short)v;
        }
        
        int ns = tx_stream->send(ChOut, SPB, md, 0.1);
        NumBatches++;

        if (first) {
            md.start_of_burst = false;
            first = false;
        }
    }
    std::cout << "[TransmitterManager] OperateThread exiting..." << std::endl;
}

void TransmitterManager::StartThreads()
{
    TxThread = std::thread(&TransmitterManager::OperateThread, this);
}
void TransmitterManager::StopThreads()
{
    StopAll = true;
    if(TxThread.joinable()) {
        TxThread.detach();
    }
}