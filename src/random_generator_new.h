#pragma once
#include <random>
#include <immintrin.h>
#if defined(_MSC_VER)
#define ALIGNED_(x) __declspec(align(x))
typedef unsigned __int64 uint64;
#else
#if defined(__GNUC__)
#define ALIGNED_(x) __attribute__ ((aligned(x)))
typedef unsigned long uint64;

#endif
#endif


class random_generator_new
{
	std::mt19937 gen{};
	std::normal_distribution<float> d{ 0,1.0 };
	std::normal_distribution<double> d1{ 0,1.0 };


public:
	random_generator_new();
	~random_generator_new();
	void set_seed(unsigned int Seed)
	{
		gen.seed(Seed);
	}
	void randn(float * output, unsigned int length);
	void randn(double * output, unsigned int length);

	void random(float * output, unsigned int length);

	void randb(unsigned char  *output, unsigned int Len);
	void rand_byte(unsigned char  *output, unsigned int Len);
};
