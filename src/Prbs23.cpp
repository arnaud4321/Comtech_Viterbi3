#include "Prbs23.h"
Prbs23::Prbs23()
{
	NextStates = new unsigned int[8192];
	Outputs = new unsigned char[8192*8];
	InitializeTables();
}

Prbs23::~Prbs23()
{
	delete[] NextStates;
	delete[] Outputs;
}


void Prbs23::CreateOutputs0(unsigned int &state, unsigned int Length, unsigned char *Out)
{
	for (unsigned int i = 0; i < Length; i++)
	{
		unsigned char bit = ((state >> 22) ^ (state >> 17)) & 0x01;
		Out[i] = bit;
		state = (state << 1) | bit;
		state &= 0x7FFFFF; // Keep only the lower 23 bits
	}
}

void Prbs23::InitializeTables()
{
	// This function can be used to precompute the NextStates and Outputs tables if needed.
	// Currently not implemented.
	unsigned int PtrOut = 0;
	for(unsigned int State = 0; State < 8192; State++)
	{
		unsigned int currentState = State << 10;
		unsigned int NextState = 0;
		for (unsigned int i = 0; i < 8; i++)
		{
			unsigned char bit = ((currentState >> 22) ^ (currentState >> 17)) & 0x01;
			Outputs[PtrOut++] = bit;
			currentState = (currentState << 1) | bit;
			currentState &= 0x7FFFFF; // Keep only the lower 23 bits
			NextState = (NextState << 1) | bit;

		}
	
		NextStates[State] = NextState; 
	}
}

void Prbs23::CreateOutputs(unsigned int& state, unsigned int Length, unsigned char* Out)
{
	for(unsigned int i = 0; i < Length; i+=8)
	{
		unsigned int ReducedState = ((state >> 10) & 0x1FFF); // Ensure state is within 23 bits
		unsigned int index = ReducedState <<3; // Get the upper 13 bits for indexing
		std::memcpy(Out+i, Outputs + index, 8);
		
		state = ((state<<8)& 0x7FFFFF)| NextStates[ReducedState];
		
	}
}

