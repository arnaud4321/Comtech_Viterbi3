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
    double TimeDrift;
    double FrequencyShift;
    double AccelerationPeriod;
    double TotalPeriod, StablePeriod;
    SimModes SimMode;
    unsigned int NumErrors;
    unsigned int Seed;
    bool Debug;
    SymbolRateEstimatorConfig SymRateCfg;

};
