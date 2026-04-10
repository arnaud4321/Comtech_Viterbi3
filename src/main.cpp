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
#include <csignal>
using namespace std;
#include "ConsoleAlert.h"
#include "Params.h"
#include "Transmitter.h"
#include "AWGNChannel.h"
#include "Receiver.h"
#include "Sampler.h"
#include "UHDSampler.h"
#include "TransmitterManager.h"
#include "USRPInit.h"
using namespace std;

std::atomic<bool> Finish{false};

mutex mtxfilethr;

void signal_handler(int signal)
{
    if (signal == SIGINT || signal == SIGTERM) {
        // Do not use std::cout here, as it is not async-signal-safe and can cause deadlocks!
        Finish.store(true, std::memory_order_relaxed);
    }
}

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
	
	std::signal(SIGINT, signal_handler);
	std::signal(SIGTERM, signal_handler);

 	Params objParams;
    objParams.ReadParams(argv[1]);

	Transmitter* oTx = nullptr;
	if (objParams.OpMode != RX_ONLY) {
		oTx = new Transmitter();
		oTx->StartThreads(PRBS_TX,objParams.RollOff);
	}
	
	AWGNChannel* oAWGN = nullptr;
	if (objParams.OpMode == NOT_OP) {
		oAWGN = new AWGNChannel(objParams.Seed);
		oAWGN->SetTransmitter(oTx);
		oAWGN->SetParameters(&objParams);
		oAWGN->StartThreads();
	}

	USRPInit* pUSRPInit = nullptr;
	TransmitterManager* pTxManager = nullptr;
	Sampler* objSampler = nullptr;

	if (objParams.OpMode != NOT_OP) {
		pUSRPInit = new USRPInit(objParams.RxFreq, objParams.TxFreq, 0, objParams.TxGaindb, 0, objParams.ref, objParams.RxSampleRate, objParams.TxSampleRate);
	}

	if (objParams.OpMode == NOT_OP) {
		objSampler = new Sampler(objParams.Debug);
		objSampler->SetChannel(oAWGN);
	} else if (objParams.OpMode != TX_ONLY) {
		objSampler = new UHDSampler(pUSRPInit->usrp);
	}

	if (objParams.OpMode == TX_RX || objParams.OpMode == TX_ONLY) {
		pTxManager = new TransmitterManager(pUSRPInit->usrp);
		pTxManager->SetTransmitter(oTx);
		std::this_thread::sleep_for(100ms);
		pTxManager->StartThreads();
	}

	if (objSampler) {
		const double tpm = (objParams.DisplayPeriodSec > 0.0) ? objParams.DisplayPeriodSec : 1.0;
		objSampler->SetThroughputMeasurePeriodSec(tpm);
		objSampler->SetDisplayPeriodSec(objParams.DisplayPeriodSec);
		objSampler->StartThread();
	}

	Receiver* oRx = nullptr;
	if (objParams.OpMode != TX_ONLY) {
		oRx = new Receiver();
		oRx->SetSampler(objSampler);
		oRx->SetChannel(oAWGN);
		oRx->SetTransmitter(oTx);
		oRx->StartThreads(objParams.RollOff, objParams.TxMode, objParams.SymRateCfg, objParams.TimingTrackingCfg,
		                  objParams.DisplayPeriodSec, objParams.CentralFreqCfg, objParams.PhaseTrackingCfg, objParams.ConstellationCfg);
	}
	
	auto Start = std::chrono::high_resolution_clock::now();

	condition_variable *pcv_rx_out, *pcv_out_rx;
	if(oRx != nullptr && objParams.TxMode == FILE_TX)
	{
		unsigned char *Out = nullptr;
		
		while((Out = oRx->GetOutput()) == nullptr)
		{
			std::this_thread::sleep_for(std::chrono::milliseconds(1));
		}
		oRx->AdvanceOutput();
	}
	else if (oRx != nullptr)
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
			
            double snrDb = std::numeric_limits<double>::quiet_NaN();
            if (oAWGN) {
                auto awgnState = oAWGN->GetCurrentApplied();
                snrDb = awgnState.SnrDb;
            }
			
			RxStatistics stats = oRx->GetRxStatistics(elapsed.count());
			
#ifdef ENABLE_RAW_PREVITERBI_METRICS
			const bool rawLocked = oRx->RawSyncLocked.load(std::memory_order_relaxed);
			cout << std::fixed << std::setprecision(6)
			     << "[Main][RAW] t_sim=" << elapsed.count() << " s"
			     << " t_rate=" << stats.t_rate << " s"
			     << " snrDb=" << snrDb
			     << " snrEvmDb=" << stats.snrFromEvmDb
			     << " RAW_LOCK=" << (rawLocked ? "\033[32mLOCKED\033[0m" : "\033[31mUNLOCKED\033[0m")
			     << " syms=" << oRx->RawNumSymsAll
			     << " symErr=" << oRx->RawNumSymErrorsAll
			     << " SER=" << std::scientific << std::setprecision(6)
			     << (static_cast<double>(oRx->RawNumSymErrorsAll) / std::max(1.0, static_cast<double>(oRx->RawNumSymsAll)))
			     << std::defaultfloat
			     << " syncPeakAbs=" << std::scientific << std::setprecision(6) << oRx->RawSyncPeakAbs.load(std::memory_order_relaxed)
			     << std::defaultfloat
			     << " syncPeakToMean=" << std::scientific << std::setprecision(6) << oRx->RawSyncBestAbs.load(std::memory_order_relaxed)
			     << std::defaultfloat
			     << " syncThrPeakToMean=" << std::scientific << std::setprecision(6) << oRx->RawSyncThrAbs.load(std::memory_order_relaxed)
			     << std::defaultfloat
			     << " syncPhaseRad=" << std::fixed << std::setprecision(6) << oRx->RawSyncPeakPhaseRad.load(std::memory_order_relaxed)
			     << " syncAppliedPhaseRad=" << std::fixed << std::setprecision(6) << oRx->RawSyncAppliedPhaseRad.load(std::memory_order_relaxed)
			     << " syncAppliedQuad=" << oRx->RawSyncAppliedQuad.load(std::memory_order_relaxed)
			     << " syncLagSym=" << oRx->RawSyncLagSym.load(std::memory_order_relaxed)
			     << " bits=" << oRx->RawNumBitsAll
			     << " bitErr=" << oRx->RawNumErrorsAll
			     << " rawBER=" << std::scientific << std::setprecision(6)
			     << (static_cast<double>(oRx->RawNumErrorsAll) / std::max(1.0, static_cast<double>(oRx->RawNumBitsAll)))
			     << std::defaultfloat
			     << endl;
#endif
		cout << std::fixed << std::setprecision(6)
			     << "[Main] t_sim=" << stats.t_sim << " s"
			     << " t_rate=" << stats.t_rate << " s"
			     << " snrDb=" << snrDb
			     << " snrEvmDb=" << stats.snrFromEvmDb
			     << " evmRms=" << stats.evmRms
			     << " samplerMsps=" << stats.samplerMsps
			     << " rxFilterMsps=" << stats.rxFilterMsps
			     << " symRateLk=" << stats.isSymRateLocked
			     << " symRateMsps=" << stats.symRateMsps
			     << " centralRdy=" << stats.isCentralFreqReady
			     << " ncoHz=" << stats.centralNcoHz
			     << " gainDb=" << stats.gainDb
			     << " gardnerLk=" << stats.isGardnerLocked
			     << " gardnerPpm=" << stats.gardnerPpm
			     << " phaseLk=" << stats.isPhaseLocked
			     << " phaseFreqHz=" << stats.phaseFreqHz
			     << " vitLk=" << stats.isViterbiLocked
			     << " prbsLk=" << stats.isPrbsLocked
			     << " bits=" << stats.numBits
			     << " errors=" << stats.numErrors
			     << " BER=" << std::scientific << std::setprecision(6)
			     << stats.ber
			     << std::defaultfloat;
		if (stats.uhdOverflows > 0)
			cout << " uhdOverflows=" << stats.uhdOverflows;
		cout << endl;
			#ifdef DEBUG_STATISTICS
			cout << std::fixed << std::setprecision(6)
			     << "[Viterbi] mean_metrics_growth=" << oRx->CurrDebugStatistics.MeanMetricsGrowth
			     << " max_metrics_growth=" << oRx->CurrDebugStatistics.MaxMetricsGrowth
			     << " t_sim=" << stats.t_sim << " s"
			     << " t_rate=" << stats.t_rate << " s"
			     << endl;
			#endif
		}
	}
    else {
        while(!Finish.load(std::memory_order_relaxed)) {
            std::this_thread::sleep_for(std::chrono::duration<double>(1.0));
        }
    }

    std::cout << "[Main] Stopping Receiver..." << std::endl;
    if (oRx) oRx->StopThreads();
    std::cout << "[Main] Stopping TransmitterManager..." << std::endl;
    if (pTxManager) pTxManager->StopThreads();
    std::cout << "[Main] Stopping Sampler..." << std::endl;
    if (objSampler) objSampler->StopThread();
    std::cout << "[Main] Stopping AWGNChannel..." << std::endl;
    if (oAWGN) oAWGN->StopThreads();
    std::cout << "[Main] Stopping Transmitter..." << std::endl;
    if (oTx) oTx->StopThreads();

    // Give detached UHD threads a moment to finish cleanly before we yank the memory
    std::this_thread::sleep_for(std::chrono::milliseconds(500));

    std::cout << "[Main] Shutdown complete." << std::endl;
    // We intentionally leak the top-level objects (oRx, pTxManager, etc.) 
    // because detaching UHD threads means they might still be executing.
    // The OS will reclaim all memory and USB handles immediately upon exit.
    std::exit(0);

	return 0;
}

