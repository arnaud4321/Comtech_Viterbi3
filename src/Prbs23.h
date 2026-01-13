#pragma once
#include <cstring>
#include <cstdint>

class Prbs23
{
private:
	unsigned char* Outputs;
	unsigned int* NextStates;
	void InitializeTables();
public:
	Prbs23();
	~Prbs23();
	void CreateOutputs0(unsigned int &state, unsigned int Length, unsigned char *Out);
	void CreateOutputs(unsigned int& state, unsigned int Length, unsigned char* Out);

	

};

