#pragma once
#include <cstring>
#include <cstdint>

class Prbs23
{
private:
	unsigned char* Outputs;
	unsigned int* NextStates;
	void InitializeTables();
	void CreateOutputs0(unsigned int &state, unsigned int Length, unsigned char *Out);

public:
	Prbs23();
	~Prbs23();
	void CreateOutputs(unsigned int& state, unsigned int Length, unsigned char* Out);

	

};

