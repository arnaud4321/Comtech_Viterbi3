// Comtech_Viterbi_3_Utilities.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "definitions.h"
#include <iostream>
#include <thread>
#include <chrono>
#include <condition_variable>
#include <mutex>
using namespace std;
#include "Params.h"
#include "Transmitter.h"
#include "AWGNChannel.h"
#include "Receiver.h"
#include "Sampler.h"
using namespace std;

bool Finish = false;

mutex mtxfilethr;
int main(int argc, char* argv[])
{
	if(argc != 2)
	{
		cout<<"File configuration is needed. Exiting ..."<<endl;
		exit(-1);
	}
	#ifdef WRITE_LOG_THR
	FILE *fidthr = fopen("LogThreadsInfo.txt","wt");
    fprintf(fidthr,"Main Thread %d\n",gettid());
    fclose(fidthr);
	#endif
	cout<<"Starting program"<<endl;
	

 	Params objParams;
    objParams.ReadParams(argv[1]);
	Transmitter oTx;
	oTx.StartThreads(PRBS_TX,objParams.RollOff);
	double TxPower = oTx.GetTxPower();
	AWGNChannel oAWGN(objParams.Seed);
	oAWGN.SetTransmitter(&oTx);
	oAWGN.SetParameters(&objParams);
	oAWGN.StartThreads();
	Sampler objSampler(objParams.Debug);
	objSampler.SetChannel(&oAWGN);
	objSampler.StartThread();
	Receiver oRx;
	oRx.SetSampler(&objSampler);
	oRx.StartThreads(objParams.RollOff, objParams.TxMode, objParams.SymRateCfg);
	
	auto Start = std::chrono::high_resolution_clock::now();

	condition_variable *pcv_rx_out, *pcv_out_rx;
	if(objParams.TxMode == FILE_TX)
	{
		unsigned char *Out;
		
		while((Out == oRx.GetOutput())== 0)
		{
			std::mutex mtx;
            std::unique_lock<std::mutex> lck(mtx);
        	pcv_rx_out->wait(lck);
		}
		oRx.AdvanceOutput();
		pcv_out_rx->notify_one();
	}
	else
	{
		int n = 0;
		while(1)
		{
			std::this_thread::sleep_for((std::chrono::duration<double>(1)));
  			auto Now = std::chrono::high_resolution_clock::now();
        	std::chrono::duration<double> elapsed = Now - Start;
			cout<<"Elapsed Time "<<elapsed.count()<<endl;
			cout<<"EsN0 "<<objParams.EsN0<<" Number of Decoded Bits "<<oRx.NumBitsAll<<" Number of Errors "<<oRx.NumErrorsAll<<endl;
			cout<<"BER "<<double(oRx.NumErrorsAll)/double(oRx.NumBitsAll)<<endl;
			#ifdef DEBUG_STATISTICS
			cout<<"Average Metrics Growth "<<oRx.CurrDebugStatistics.MeanMetricsGrowth<<endl;
			cout<<"Max Metrics Growth "<<oRx.CurrDebugStatistics.MaxMetricsGrowth<<endl;

			#endif
		}
	}

	oTx.StopThreads();
	oAWGN.StopThreads();
	objSampler.StopThread();
	oRx.StopThreads();
	return 0;
}

