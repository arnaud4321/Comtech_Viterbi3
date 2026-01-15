// Comtech_Viterbi_3_Utilities.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <thread>
#include <chrono>
using namespace std;
#include "Params.h"
#include "Transmitter.h"
#include "AWGNChannel.h"
#include "Receiver.h"
#include "Sampler.h"
using namespace std;

bool Finish = false;


int main(int argc, char* argv[])
{
	if(argc != 2)
	{
		cout<<"File configuration is needed. Exiting ..."<<endl;
		exit(-1);
	}

	
 	Params objParams;
    objParams.ReadParams(argv[1]);
	Transmitter oTx;
	oTx.StartThreads(PRBS_TX,objParams.RollOff);
	double TxPower = oTx.GetTxPower();
	AWGNChannel oAWGN(objParams.Seed);
	oAWGN.SetTransmitter(&oTx);
	oAWGN.StartThreads(objParams.EsN0,objParams.TimeDrift,objParams.FrequencyShift, objParams.RandomFrequency);
	Sampler objSampler(objParams.Debug);
	objSampler.SetChannel(&oAWGN);
	objSampler.StartThread();
	Receiver oRx;
	oRx.SetSampler(&objSampler);
	oRx.StartThreads(objParams.RollOff);
	std::this_thread::sleep_for(std::chrono::seconds(1000));
	oTx.StopThreads();
	return 0;
}

