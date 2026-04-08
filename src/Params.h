/**
 * @file Params.h
 * @brief JSON-driven simulation and receiver configuration (@ref Params::ReadParams).
 */
#pragma once
#include <string>
using namespace std;
#include "definitions.h"
#include <immintrin.h>
#include "SymbolRateEstimator.h"
#include "ReceiverFreqCorrector.h"
#include "ConstellationDisplay.h"

/**
 * @brief Runtime parameters loaded from a JSON configuration file.
 *
 * Paths and numeric fields map to `config.json` sections (Channel, Simulation, SymbolRate, etc.).
 */
class Params
{
public:
	Params();
	/** @brief Populate fields from JSON at @a FileName. */
	void ReadParams(string FileName);
    TxModes TxMode;
    string TxFileName;
    double EsN0;
    double RollOff;
    double InitialGainDb = 0.0;
    /// Channel gain span (dB) over ramp segments (added to InitialGainDb).
    double DynamicRangeDb = 0.0;
    double InitialFrequencyShift = 0.0;
    double FrequencyShift;
    /// Fixed ppm: Tx/Rx oscillator mismatch on the symbol clock.
    double ClockMismatchPpm;
    /// Carrier frequency / symbol rate R_s (dimensionless). Maps RF shift (FrequencyShift, Hz)
    /// to relative error on R_s: ppm_doppler ≈ FrequencyShift·10⁶ / (CarrierToSymbolRateRatio · R_s).
    double CarrierToSymbolRateRatio;
    double AccelerationPeriod;
    double TotalPeriod, StablePeriod;
    /// If true, applies channel gain directly to the signal before adding noise (varies SNR).
    bool ApplyGainBeforeNoise = false;
    SimModes SimMode;
    unsigned int NumErrors;
    unsigned int Seed;
    bool Debug;
    /// Console status: >0 = periodic interval (s); <=0 = events only (lock/unlock, etc.), no periodic lines.
    double DisplayPeriodSec = 1.0;
    SymbolRateEstimatorConfig SymRateCfg;
    ReceiverFreqCorrectorConfig CentralFreqCfg;
    ConstellationDisplayConfig ConstellationCfg;

};
