#include "NoiseGen.h"

NoiseGen::NoiseGen()
{
	std::random_device rd;  // uses /dev/urandom on Linux
    unsigned int tt = rd();

	/* Create mt19937 engine */
	engine = new svrng_engine_t(svrng_new_mt19937_engine(tt));

	distr1 = new svrng_distribution_t(svrng_new_normal_distribution_float(0.0, 1.0));


}


NoiseGen::~NoiseGen()
{
	svrng_delete_distribution(*distr1);
	delete engine;
	delete distr1;

}

void NoiseGen::CreateOutputs_svrng(float *In, int Length)
{
	

	for (int ii = 0; ii < Length; ii += 16)
	{
		svrng_float16_t z = svrng_generate16_float(*engine, *distr1);
		std::memcpy(In+ii,&z,16);
	}
}