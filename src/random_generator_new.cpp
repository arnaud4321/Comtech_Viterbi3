#include "random_generator_new.h"


random_generator_new::random_generator_new()
{
}


random_generator_new::~random_generator_new()
{
}

void random_generator_new::randn(float *output, unsigned int Length)
{
	for (int ii = 0; ii < Length; ii++)
	{
		output[ii] = d(gen);
	}
}
void random_generator_new::randn(double* output, unsigned int Length)
{
	for (int ii = 0; ii < Length; ii++)
	{
		output[ii] = d1(gen);
	}
}

void random_generator_new::random(float *output, unsigned int Length)
{
	float k = 2.3283064365386962890625e-10;
	__m256 mk = _mm256_set1_ps(k);
	__m256 mhalf = _mm256_set1_ps(0.5);
	unsigned int Length8 = (Length >> 3) << 3;
	if (Length > Length8)
		Length8 += 8;
		
	ALIGNED_(32) unsigned int Temp[8];
	
	for (int ii = 0; ii < Length8; ii+=8)
	{
		for (int jj = 0; jj < 8; jj++)
		{
			Temp[jj] = gen();
		}

		__m256i mRand = _mm256_load_si256((__m256i *)Temp);
		__m256 mRandf = _mm256_cvtepi32_ps(mRand);
		mRandf = _mm256_fmadd_ps(mRandf, mk,mhalf);
		_mm256_storeu_ps(output + ii, mRandf);
	}

	//for (int ii = Length8; ii < Length; ii++)
	//{
	//	output[ii] = float(gen())*2.3283064365386962890625e-10;
	//}

}
void random_generator_new::rand_byte(unsigned char *output, unsigned int length)
{
	unsigned int tmp;
	for (int ii = 0; ii < length; ii++)
	{
		tmp = gen();
		tmp = (tmp >> 24);
		output[ii] = tmp;
	}
}

void random_generator_new::randb(unsigned char *output, unsigned int length)
{
	unsigned int tmp;
	for (int ii = 0; ii < length; ii++)
	{
		tmp = gen();
		tmp = (tmp >> 31);
		output[ii] = (char)tmp;
	}


};