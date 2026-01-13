#pragma once
#include <svrng.h>
#include <random>
#include <cstring>

class NoiseGen
{
	svrng_engine_t   *engine;
	svrng_distribution_t *distr1;

public:
	NoiseGen();
	~NoiseGen();
	void CreateOutputs_svrng(float *In, int Length);
	void UpdateStd(float Std)
	{
		svrng_update_normal_distribution_float(*distr1, 0, Std);
	}


	
	


};
