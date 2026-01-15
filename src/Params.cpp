#include "Params.h"
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

    unsigned int Tmp;
    Tmp = j["Transmitter"]["Method"];
    if(Tmp == 0)
        TxMode = PRBS_TX;
    else
        TxMode = FILE_TX;
    RollOff = j["Transmitter"]["RollOff"];

    TxFileName = j["Transmitter"]["FileName"];
	EsN0 = j["Channel"]["Esn0"];
    TimeDrift =  j["Channel"]["Drift"];
    FrequencyShift = j["Channel"]["FrequencyShift"];
    Tmp = j["Channel"]["RandomFrequency"];
    if(Tmp == 0)
        RandomFrequency = false;
    else
        RandomFrequency = true;

    Tmp = j["Simulation"]["Method"];
    if(Tmp == 0)
        SimMode = ACQ_SIM;
    else
        SimMode = CONT_SIM;
    
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

}