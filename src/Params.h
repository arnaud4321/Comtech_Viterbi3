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
#include "ReceiverTimingTracking.h"
#include "ReceiverFreqCorrector.h"
#include "ConstellationDisplay.h"

/**
 * @brief Runtime parameters loaded from a JSON configuration file.
 *
 * JSON layout: @c Transmitter, @c Receiver, @c Channel, @c Operation (@c Mode, @c ref), @c Simulation.
 * Sous-objets RX peuvent être sous @c Receiver ou à la racine (compat.). @c Simulation.OperationMode /
 * @c ClockReference restent des repli si @c Operation est absent.
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
    OpModes OpMode;
    double RxFreq, TxFreq, TxGaindb;
    double RxSampleRate, TxSampleRate;
    string ref;
    unsigned int NumErrors;
    unsigned int Seed;
    bool Debug;
    /// When @c OpMode is @c NOT_OP, replay int16 interleaved IQ from disk instead of AWGN+TX (see @ref FileSampler).
    bool IqFileReplayEnable = false;
    /// Path to binary IQ file (short I,Q,...), same format as @c InputRecording.bin from @ref UHDSampler.
    string IqFileReplayPath;
    /// If true, rewind at EOF; if false, stop the file producer at EOF.
    bool IqFileReplayLoop = true;
    /// Console status: >0 = periodic interval (s); <=0 = events only (lock/unlock, etc.), no periodic lines.
    double DisplayPeriodSec = 1.0;
    SymbolRateEstimatorConfig SymRateCfg;
    TimingTrackingConfig TimingTrackingCfg;
    ReceiverFreqCorrectorConfig CentralFreqCfg;
    PhaseTrackingConfig PhaseTrackingCfg;
    ConstellationDisplayConfig ConstellationCfg;

};
