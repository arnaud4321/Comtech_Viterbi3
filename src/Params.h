#pragma once
#include <string>
using namespace std;
#include "definitions.h"
#include <immintrin.h>
#include "SymbolRateEstimator.h"

class Params
{
public:
	Params();
	void ReadParams(string FileName);
    TxModes TxMode;
    string TxFileName;
    double EsN0;
    double RollOff;
    double FrequencyShift;
    /// Fixed ppm: Tx/Rx oscillator mismatch on the symbol clock.
    double ClockMismatchPpm;
    /// Carrier frequency / symbol rate R_s (dimensionless). Maps RF shift (FrequencyShift, Hz)
    /// to relative error on R_s: ppm_doppler ≈ FrequencyShift·10⁶ / (CarrierToSymbolRateRatio · R_s).
    double CarrierToSymbolRateRatio;
    double AccelerationPeriod;
    double TotalPeriod, StablePeriod;
    SimModes SimMode;
    unsigned int NumErrors;
    unsigned int Seed;
    bool Debug;
    SymbolRateEstimatorConfig SymRateCfg;

};
