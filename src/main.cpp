/**
 * @file main.cpp
 * @brief Program entry: JSON config, start TX / channel / sampler / RX, wait for @c Finish, join threads.
 *
 * @details Typical flow: @ref Params::ReadParams, construct @ref Transmitter / @ref AWGNChannel / @ref Sampler /
 * @ref Receiver, wire pointers (@c SetTransmitter, @c SetChannel, @c SetSampler), @c StartThreads on each block
 * with rolloff and RX options from JSON. Main thread loop advances TX output, optionally drains @c Receiver
 * output for FILE mode, prints periodic status until @c Finish. Stops child threads in reverse dependency order
 * (RX → sampler → channel → TX).
 */

#include "definitions.h"
#include <iomanip>
#include <iostream>
#include <thread>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <atomic>
using namespace std;
#include "ConsoleAlert.h"
#include "Params.h"
#include "Transmitter.h"
#include "AWGNChannel.h"
#include "Receiver.h"
#include "Sampler.h"
using namespace std;

std::atomic<bool> Finish{false};

mutex mtxfilethr;

/**
 * @brief Loads config, starts TX/channel/sampler/RX threads, runs until @c Finish or FILE_TX drain.
 * @param argc Argument count (expects @c 2: program name + path to JSON config).
 * @param argv @a argv[1] = configuration file path.
 * @return 0 on normal shutdown, -1 if config path missing.
 */
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
	{
		// Match RX filter: 1 s rollup when DisplayPeriodSec<=0 (avoid coupling to small Constellation PeriodSec).
		const double tpm =
		    (objParams.DisplayPeriodSec > 0.0) ? objParams.DisplayPeriodSec : 1.0;
		objSampler.SetThroughputMeasurePeriodSec(tpm);
	}
	objSampler.SetDisplayPeriodSec(objParams.DisplayPeriodSec);
	objSampler.SetChannel(&oAWGN);
	objSampler.StartThread();
	Receiver oRx;
	oRx.SetSampler(&objSampler);
	oRx.SetChannel(&oAWGN);
	oRx.SetTransmitter(&oTx);
	oRx.StartThreads(objParams.RollOff, objParams.TxMode, objParams.SymRateCfg, objParams.DisplayPeriodSec, objParams.CentralFreqCfg, objParams.ConstellationCfg);
	
	auto Start = std::chrono::high_resolution_clock::now();

	condition_variable *pcv_rx_out, *pcv_out_rx;
	if(objParams.TxMode == FILE_TX)
	{
		unsigned char *Out = nullptr;
		
		while((Out = oRx.GetOutput()) == nullptr)
		{
			std::this_thread::sleep_for(std::chrono::milliseconds(1));
		}
		oRx.AdvanceOutput();
	}
	else
	{
		while(!Finish.load(std::memory_order_relaxed))
		{
			// DisplayPeriodSec <= 0: no periodic [Main] line; sleep to avoid busy-wait.
			const double dps = objParams.DisplayPeriodSec;
			const double sleepSec = (dps > 0.0) ? dps : 1.0;
			std::this_thread::sleep_for(std::chrono::duration<double>(sleepSec));
			if (dps <= 0.0)
				continue;
			auto Now = std::chrono::high_resolution_clock::now();
			std::chrono::duration<double> elapsed = Now - Start;
			const double t_rate = static_cast<double>(oRx.NumBitsAll) / SymbolRate;
			auto awgnState = oAWGN.GetCurrentApplied();
			double snrDb = awgnState.SnrDb;
			double evmRms = oRx.GetPhaseTrackingEvmRms();
			double snrEvmDb = (evmRms > 1e-12) ? (-20.0 * std::log10(evmRms)) : std::numeric_limits<double>::quiet_NaN();
#ifdef ENABLE_RAW_PREVITERBI_METRICS
			const bool rawLocked = oRx.RawSyncLocked.load(std::memory_order_relaxed);
			cout << std::fixed << std::setprecision(6)
			     << "[Main][RAW] t_sim=" << elapsed.count() << " s"
			     << " t_rate=" << t_rate << " s"
			     << " snrDb=" << snrDb
			     << " snrEvmDb=" << snrEvmDb
			     << " RAW_LOCK=" << (rawLocked ? (std::string(ConsoleAlert::kGreenOpen) + "1" + ConsoleAlert::kReset) : "0")
			     << " syms=" << oRx.RawNumSymsAll
			     << " symErr=" << oRx.RawNumSymErrorsAll
			     << " SER=" << std::scientific << std::setprecision(6)
			     << (static_cast<double>(oRx.RawNumSymErrorsAll) / std::max(1.0, static_cast<double>(oRx.RawNumSymsAll)))
			     << std::defaultfloat
			     << " syncPeakAbs=" << std::scientific << std::setprecision(6) << oRx.RawSyncPeakAbs.load(std::memory_order_relaxed)
			     << std::defaultfloat
			     << " syncPeakToMean=" << std::scientific << std::setprecision(6) << oRx.RawSyncBestAbs.load(std::memory_order_relaxed)
			     << std::defaultfloat
			     << " syncThrPeakToMean=" << std::scientific << std::setprecision(6) << oRx.RawSyncThrAbs.load(std::memory_order_relaxed)
			     << std::defaultfloat
			     << " syncPhaseRad=" << std::fixed << std::setprecision(6) << oRx.RawSyncPeakPhaseRad.load(std::memory_order_relaxed)
			     << " syncAppliedPhaseRad=" << std::fixed << std::setprecision(6) << oRx.RawSyncAppliedPhaseRad.load(std::memory_order_relaxed)
			     << " syncAppliedQuad=" << oRx.RawSyncAppliedQuad.load(std::memory_order_relaxed)
			     << " syncLagSym=" << oRx.RawSyncLagSym.load(std::memory_order_relaxed)
			     << " bits=" << oRx.RawNumBitsAll
			     << " bitErr=" << oRx.RawNumErrorsAll
			     << " rawBER=" << std::scientific << std::setprecision(6)
			     << (static_cast<double>(oRx.RawNumErrorsAll) / std::max(1.0, static_cast<double>(oRx.RawNumBitsAll)))
			     << std::defaultfloat
			     << endl;
#endif
			cout << std::fixed << std::setprecision(6)
			     << "[Main] t_sim=" << elapsed.count() << " s"
			     << " t_rate=" << t_rate << " s"
			     << " snrDb=" << snrDb
			     << " snrEvmDb=" << snrEvmDb
			     << " bits=" << oRx.NumBitsAll
			     << " errors=" << oRx.NumErrorsAll
			     << " BER=" << std::scientific << std::setprecision(6)
			     << (static_cast<double>(oRx.NumErrorsAll) / std::max(1.0, static_cast<double>(oRx.NumBitsAll)))
			     << std::defaultfloat
			     << endl;
			#ifdef DEBUG_STATISTICS
			cout << std::fixed << std::setprecision(6)
			     << "[Viterbi] mean_metrics_growth=" << oRx.CurrDebugStatistics.MeanMetricsGrowth
			     << " max_metrics_growth=" << oRx.CurrDebugStatistics.MaxMetricsGrowth
			     << " t_sim=" << elapsed.count() << " s"
			     << " t_rate=" << t_rate << " s"
			     << endl;

			#endif
		}
	}

	oTx.StopThreads();
	oAWGN.StopThreads();
	oRx.StopThreads();
	objSampler.StopThread();
	return 0;
}

