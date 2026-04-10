/** @file Params.cpp
 *  @brief JSON configuration loader for @ref Params.
 */
#include "Params.h"
#include "json.hpp"
#include <fstream>

Params::Params()
{

}

void Params::ReadParams(string FileName)
{
	using json = nlohmann::json;
	// read a JSON file
	std::ifstream i(FileName);
	json j;
	i >> j;

    const json opLegacy = j.value("Operation", json::object());
    const json recv = j.value("Receiver", json::object());
    const json sim = j.value("Simulation", json::object());

    unsigned int Tmp;
    const json& tx = j.at("Transmitter");
    Tmp = tx["Method"];
    if(Tmp == 0)
        TxMode = PRBS_TX;
    else
        TxMode = FILE_TX;
    RollOff = tx["RollOff"];

    TxFileName = tx["FileName"];
    TxFreq = tx.value("FreqHz", opLegacy.value("TxFreq", 2000e6));
    TxGaindb = tx.value("GainDb", opLegacy.value("TxGaindb", 20.0));
    TxSampleRate = tx.value("SampleRate", opLegacy.value("TxSampleRate", 21.42e6));

	EsN0 = j["Channel"]["Esn0"];
    // Channel gain profile (dB). Backward-compatible with older key InitialRangeDb.
    InitialGainDb = j["Channel"].value("InitialGainDb", j["Channel"].value("InitialRangeDb", 0.0));
    // Keep name "DynamicRangeDb": it is a dB span applied over the ramp segments.
    DynamicRangeDb = j["Channel"].value("DynamicRangeDb", 0.0);
    ApplyGainBeforeNoise = j["Channel"].value("ApplyGainBeforeNoise", false);
    InitialFrequencyShift = j["Channel"].value("InitialFrequencyShift", 0.0);
    FrequencyShift = j["Channel"]["FrequencyShift"];
    ClockMismatchPpm = j["Channel"].value("ClockMismatchPpm", 0.0);
    CarrierToSymbolRateRatio = j["Channel"].value("CarrierToSymbolRateRatio", 935.0);
    if (CarrierToSymbolRateRatio <= 0.0)
        CarrierToSymbolRateRatio = 935.0;
    AccelerationPeriod = j["Channel"]["AccelerationPeriod"];
    StablePeriod = j["Channel"]["StablePeriod"];
    TotalPeriod = 2*(AccelerationPeriod + StablePeriod);
    
    

    Tmp = j["Simulation"]["Method"];
    if(Tmp == 0)
        SimMode = ACQ_SIM;
    else
        SimMode = CONT_SIM;

    // Prefer Operation.Mode / ref; fallback Simulation.OperationMode / ClockReference for older JSON.
    Tmp = opLegacy.value("Mode", sim.value("OperationMode", 0));
    switch (Tmp) {
        case 0: OpMode = NOT_OP; break;
        case 1: OpMode = TX_ONLY; break;
        case 2: OpMode = RX_ONLY; break;
        case 3: OpMode = TX_RX; break;
        default: OpMode = NOT_OP; break;
    }
    RxFreq = recv.value("FreqHz", opLegacy.value("RxFreq", 2000e6));
    RxSampleRate = recv.value("SampleRate", opLegacy.value("RxSampleRate", 21.42e6));
    ref = opLegacy.value("ref", sim.value("ClockReference", std::string("internal")));

    SamplingFrequency = RxSampleRate;
    SymbolRate = SamplingFrequency / 2.0;
    SamplingRate = SamplingFrequency;
    TxOutputBatchDuration = (static_cast<double>(TxOutputBatchSize) / 2.0) / SamplingRate;
    
	NumErrors = j["Simulation"]["NumErrors"];
	Seed = j["Simulation"]["Seed"];

	if (Seed == 0)
	{
		_rdrand32_step(&Seed);
	}
    Tmp = j["Simulation"]["Debug"];
    if(Tmp == 0)
        Debug = false;
    else
        Debug = true;

    DisplayPeriodSec = j["Simulation"].value("DisplayPeriodSec", 1.0);
    // > 0: periodic status interval (s). <= 0: no periodic logs; lock/unlock (and similar) events only.

    const json* tTiming = nullptr;
    if (recv.contains("TimingTracking"))
        tTiming = &recv["TimingTracking"];
    else if (j.contains("TimingTracking"))
        tTiming = &j["TimingTracking"];
    if (tTiming != nullptr)
    {
        const auto& t = *tTiming;
        TimingTrackingCfg.NominalOmega = t.value("NominalOmega", TimingTrackingCfg.NominalOmega);
        TimingTrackingCfg.Kp = t.value("Kp", TimingTrackingCfg.Kp);
        TimingTrackingCfg.Ki = t.value("Ki", TimingTrackingCfg.Ki);
        int upd = t.value("UpdatePeriodSymbols", TimingTrackingCfg.UpdatePeriodSymbols);
        TimingTrackingCfg.UpdatePeriodSymbols = (upd < 1) ? 1 : upd;
    }

    const json* pPhase = nullptr;
    if (recv.contains("PhaseTracking"))
        pPhase = &recv["PhaseTracking"];
    else if (j.contains("PhaseTracking"))
        pPhase = &j["PhaseTracking"];
    if (pPhase != nullptr)
    {
        const auto& p = *pPhase;
        PhaseTrackingCfg.Kp = p.value("Kp", PhaseTrackingCfg.Kp);
        PhaseTrackingCfg.Ki = p.value("Ki", PhaseTrackingCfg.Ki);
    }

    const json* cCentral = nullptr;
    if (recv.contains("CentralFrequency"))
        cCentral = &recv["CentralFrequency"];
    else if (j.contains("CentralFrequency"))
        cCentral = &j["CentralFrequency"];
    if (cCentral != nullptr)
    {
        const auto& c = *cCentral;
        CentralFreqCfg.Enable = c.value("Enable", 1) != 0;
        CentralFreqCfg.EstimationBlockSamples = c.value("EstimationBlockSamples", 65536);
        CentralFreqCfg.EstimatePeriodSec = c.value("EstimatePeriodSec", 1.0);
        CentralFreqCfg.Estimator.MaxOffsetHz = c.value("MaxOffsetHz", 200000.0);
        CentralFreqCfg.Estimator.PeakToMedianThreshold = c.value("PeakToMedianThreshold", 8.0);
        CentralFreqCfg.Estimator.FftSizeMultiplier = c.value("FftSizeMultiplier", 4);
        CentralFreqCfg.NcoHzEmaAlpha = c.value("NcoHzEmaAlpha", 0.2);
        CentralFreqCfg.MaxHzSlewRate = c.value("MaxHzSlewRate", 8000.0);
        CentralFreqCfg.PowerEmaAlpha = c.value("PowerEmaAlpha", 0.01);
        CentralFreqCfg.TargetAvgPower = c.value("TargetAvgPower", -1.0);
    }

    const json* sSym = nullptr;
    if (recv.contains("SymbolRateEstimator"))
        sSym = &recv["SymbolRateEstimator"];
    else if (j.contains("SymbolRateEstimator"))
        sSym = &j["SymbolRateEstimator"];
    if (sSym != nullptr)
    {
        const auto& s = *sSym;
        if (s.contains("FftSize")) SymRateCfg.FftSize = s["FftSize"];
        if (s.contains("MaxOffsetHz")) SymRateCfg.MaxOffsetHz = s["MaxOffsetHz"];
        if (s.contains("EstimatePeriodSec")) SymRateCfg.EstimatePeriodSec = s["EstimatePeriodSec"];
        if (s.contains("PeakToMedianThreshold")) SymRateCfg.PeakToMedianThreshold = s["PeakToMedianThreshold"];
        if (s.contains("MaxRelativeJump")) SymRateCfg.MaxRelativeJump = s["MaxRelativeJump"];
        if (s.contains("AllowRefinement")) SymRateCfg.AllowRefinement = s["AllowRefinement"].get<bool>();

        // Force refinement OFF if we are in loopback modes or simulation where TX/RX clocks are identical.
        if (OpMode == FILE_TX || OpMode == TX_RX) {
            SymRateCfg.AllowRefinement = false;
        }
    }

    const json* cConst = nullptr;
    if (recv.contains("ConstellationDisplay"))
        cConst = &recv["ConstellationDisplay"];
    else if (j.contains("ConstellationDisplay"))
        cConst = &j["ConstellationDisplay"];
    if (cConst != nullptr)
    {
        const auto& c = *cConst;
        ConstellationCfg.PeriodSec = c.value("PeriodSec", 0.0);
        ConstellationCfg.DrawPeriodSec = c.value("DrawPeriodSec", 0.0);
        ConstellationCfg.NumSymbols = c.value("NumSymbols", 2048);
        ConstellationCfg.MaxAbs = c.value("MaxAbs", 2.0);
        ConstellationCfg.Backend = c.value("Backend", std::string("matplotlib"));
        ConstellationCfg.PythonExe = c.value("PythonExe", std::string("python3"));
        ConstellationCfg.PythonScript = c.value("PythonScript", std::string("src/constellation_display.py"));
        ConstellationCfg.XDisplay = c.value("XDisplay", std::string(""));
        ConstellationCfg.Width = c.value("Width", 61);
        ConstellationCfg.Height = c.value("Height", 31);
        ConstellationCfg.ClearScreen = c.value("ClearScreen", 0) != 0;
        if (ConstellationCfg.NumSymbols <= 0)
            ConstellationCfg.NumSymbols = 2048;
    }

}