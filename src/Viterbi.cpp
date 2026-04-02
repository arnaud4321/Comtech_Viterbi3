/**
 * @file Viterbi.cpp
 * @brief Rate-1/2 K=7 trellis: branch metrics, SIMD ACS, survivor memory, traceback; @c Decode entry.
 *
 * @details Constructor allocates path metrics, survivor bits, decoding tables, polynomials @c 0171/@c 0133,
 * and AVX shuffle metadata per @c DecoderID (phase alignment). **Decode** — @c CalcMetrics2 soft metrics,
 * @c AcsAll / @c AcsAllNew per symbol, periodic @c traceback_mid; parameters @c ExchangeIQ, @c SignI, @c SignQ select
 * the four metric-alignment trials used before @ref Receiver locks @c ViterbiParams; @c MetricsGrowth is reported
 * for that comparison. Puncturing tables are not initialized here.
 */

#include "Viterbi.h"
#include <emmintrin.h>
#include <immintrin.h>
#include <xmmintrin.h>


Viterbi::Viterbi(int Idx):DecoderID(Idx)
{

	//Derived Constants

	intNParam = 2;
	intCLength = 7;
	CLengthM1 = intCLength - 1;
	FirstTracebackPeriod = 64;
	NumTransitions = (1 << intNParam);

	BranchMetrics = (float *)  _mm_malloc(64*sizeof(float),32); //Each time we calculate 8I 8 Q and spread the metrics across an 256 

	intNoStates = 1 << (intCLength -1);
	intMask = intNoStates - 1;
	//Memory Allocation

	vecPathMetMem = (float *) _mm_malloc (sizeof(float) * 2 * intNoStates,32 );
	//Allign to 16


	bmatSurvMem = new  uint64_t[MAX_BURST_LENGTH]; //One int64 fits a row
	MaskSurv = MAX_BURST_LENGTH -1;
	for(int i = MaskSurv; i >= MAX_BURST_LENGTH-256; i--)
		bmatSurvMem[i] = 0;
	imatDecodingTable = new int [intNoStates*2];
	ivecGenPolys = new unsigned int [intNParam];
	//Branch Metrics Calculation 

	AllBinVecs = new int [ NumTransitions *intNParam];

	//Depuncturing
	//DepuncOutputBuffer = (float *)_mm_malloc(BranchMetricsPeriod*sizeof(float), 32);

	__m256 mZero = _mm256_setzero_ps();

	//for (unsigned int ii = 0; ii < BranchMetricsPeriod; ii += 8)
	//	_mm256_store_ps(DepuncOutputBuffer + ii, mZero);


	/*for(unsigned int ii = 0; ii < intNParam; ii++)
		ivecGenPolys[ii] = GenPolysInp[ii];*/

	ivecGenPolys[0] = 0171;
	ivecGenPolys[1] = 0133;
	calc_decoding_table();
	generate_binvecs(); //Generate all possible combinations of binary vectors of length intNParam, store into AllBinVecs



	mArranged =  (__m256i *) _mm_malloc(8*32,32);
	/*mArranged[0] = _mm256_set_epi32(0,1,0,1,1,0,1,0);//1st 2 are used for phase 0
	mArranged[1] = _mm256_set_epi32(1,0,1,0,0,1,0,1);
	mArranged[2] =  _mm256_set_epi32(6,7,6,7,3,2,3,2);//Used for phase 1
	mArranged[3] =  _mm256_set_epi32(7,6,7,6,2,3,2,3);//Used for phase 1
*/
	mArranged[0] = _mm256_setr_epi32(0,1,0,1,3,2,3,2);//1st 2 are used for phase 0
	mArranged[1] = _mm256_setr_epi32(1,0,1,0,2,3,2,3);//1st 2 are used for phase 0


	mPermCalc =  (__m128i *) _mm_malloc(32,16);
	mPermCalc[0] = _mm_set_epi8(13,13,12,12,9,9,8,8,5,5,4,4,1,1,0,0);//Odd
	mPermCalc[1] = _mm_set_epi8(15,14,15,14,11,10,11,10,7,6,7,6,3,2,3,2);//Even


	//BranchMet2323, BranchMet1010, BranchMet1010, BranchMet2323
	//BranchMet1010, BranchMet2323, BranchMet2323, BranchMet1010


	mOffset = (__m256i *) _mm_malloc(2 * sizeof(__m256i), 32);
	mOffset[0] = _mm256_set_epi32(14, 12, 10, 8, 6, 4, 2, 0);
	mOffset[1] = _mm256_set_epi32(15, 13, 11, 9, 7, 5, 3, 1);

}

Viterbi::~Viterbi(void)
{
	delete [] bmatSurvMem;
	delete [] imatDecodingTable;

	delete [] AllBinVecs;
	//delete [] bvecPuncPatt;
	//_mm_free(DepuncOutputBuffer) ;
	_mm_free(BranchMetrics);
	//_mm_free(mDecsMask);
	_mm_free(vecPathMetMem);
	_mm_free(mArranged);
	_mm_free(mPermCalc);
	if (DiffPunc)
		delete[] DiffPunc;
	delete [] ivecGenPolys;
	_mm_free(mOffset);
}

void Viterbi::calc_decoding_table()
{
	for(unsigned  int CurrState0 = 0; CurrState0 < intNoStates/2; CurrState0++)
	{
		int PrevState0 = (CurrState0 << 1) & intMask; 
		int PrevState1 = PrevState0 + 1;

		int CurrState1 = CurrState0 + intNoStates/2;

		int Transition = 0;

		// What is the transition from PrevState0 to CurrState0 (the input must be zero)

		for(unsigned int jj = 0; jj < intNParam; jj++)
		{
			int TransitionBit = is_odd_ones(PrevState0 & ivecGenPolys[jj]);
			Transition = (Transition << 1) + TransitionBit;
		}
		imatDecodingTable[CurrState0 << 1] = Transition;

		//What is the transition from PrevState1 to CurrState0 (the input must be zero)

		Transition = 0;
		for(unsigned int jj = 0; jj < intNParam; jj++)
		{
			int TransitionBit = is_odd_ones(PrevState1 & ivecGenPolys[jj]);
			Transition = (Transition << 1) + TransitionBit;
		}
		imatDecodingTable[(CurrState0 << 1) +1] = Transition;

		//What is the transition from PrevState0 to CurrState1 (the input must be one)
		Transition = 0;
		for(unsigned int jj = 0; jj < intNParam; jj++)
		{
			int TransitionBit = is_odd_ones((PrevState0 + intNoStates) & ivecGenPolys[jj]);
			Transition = (Transition << 1) + TransitionBit;
		}
		imatDecodingTable[CurrState1<<1] = Transition;


		//What is the transition from PrevState1 to CurrState1 (the input must be one)
		Transition = 0;
		for( unsigned int jj = 0; jj < intNParam; jj++)
		{
			int TransitionBit = is_odd_ones((PrevState1 + intNoStates) & ivecGenPolys[jj]);
			Transition = (Transition << 1) + TransitionBit;
		}
		imatDecodingTable[(CurrState1 << 1) + 1] = Transition;

	}

}

int Viterbi::is_odd_ones(int input)
{

	int ones = 0;
	for(unsigned int ii = 0; ii < intCLength ;ii++)
	{
		ones = ones ^ (input & 0x1);
		input = input >> 1;
	}

	return ones;
}




void Viterbi::traceback_mid( unsigned char *Out, int BestMetInd, int Length)
{
	



		//Start the Traceback 
		intSurvMemRdRowAddr = (intSurvMemWrRowAddr - 1)& MaskSurv; // Init the pointer to compensate for the unecessary advance at end 
		

		int intSurvMemRdColAddr = BestMetInd; //Starting the Traceback from best State

		//Traceback Read

		unsigned int IndexBit;
		for( IndexBit = 0; IndexBit < FirstTracebackPeriod; IndexBit++ )
		{
			uint64_t AllDecisions = bmatSurvMem[intSurvMemRdRowAddr];
			AllDecisions = AllDecisions >> (intSurvMemRdColAddr);
			unsigned int NewDecision = (unsigned int) (AllDecisions &1);
			intSurvMemRdColAddr = ((intSurvMemRdColAddr << 1) & intMask) | NewDecision;
			intSurvMemRdRowAddr--;
 			intSurvMemRdRowAddr &= MaskSurv;
		}

		//Traceback Decode

	

		for(int i = Length-1; i >= 0; i--)
		{
			uint64_t AllDecisions = bmatSurvMem[intSurvMemRdRowAddr];
			AllDecisions = AllDecisions >> (intSurvMemRdColAddr);
			unsigned int NewDecision = (unsigned int) (AllDecisions &1);
			intSurvMemRdColAddr = ((intSurvMemRdColAddr << 1) & intMask) | NewDecision;
			Out[i] = (unsigned char) NewDecision;
			intSurvMemRdRowAddr--;
			intSurvMemRdRowAddr &= MaskSurv;

		}

		//Update the counter of Decoded bits

}


void Viterbi::generate_binvecs(void)
{
	for(int ii = 0; ii < NumTransitions; ii++)
	{
		dec2bin(ii,&AllBinVecs[ii*intNParam],intNParam);
	}
	for(unsigned int ii = 0; ii < (NumTransitions * intNParam); ii++ )
		if(AllBinVecs[ii] == 0)
			AllBinVecs[ii] = -1;
}

void Viterbi::dec2bin(int index, int *v, int length) 
{
	int i, bintemp = index;

	for (i= length-1; i>=0; i--) 
	{
		v[i] = (bintemp & 1);
		bintemp = (bintemp >> 1);
	}
}




void Viterbi::add_comp_select(int intNewState, float * BranchMetrics,  int *intDecisions)
{
	//We handle in parallel NewState,NewState + 1, NewState + NumStates/2, NewState + NumStates/2+1

	unsigned int intOldState0 = (intNewState << 1) & intMask;


	// BranchMetrics[0] = 0; BranchMetrics[1]=2; BranchMetrics[2] = 1; BranchMetrics[3] = 3;




	 
	__m128 BranchMetricsA = _mm_load_ps(BranchMetrics);//Starts from [3]..[0]
	
	__m128 PathMetrics = _mm_load_ps(&vecPathMetMem[intOldState0+intOldPathMetMemPtr]);

	__m128 BranchMetricsB = BranchMetricsA;

	__m128 PathMetricsA = _mm_shuffle_ps(PathMetrics,PathMetrics,_MM_SHUFFLE(2,0,2,0));

	  BranchMetricsB = _mm_shuffle_ps(BranchMetricsB,BranchMetricsA,_MM_SHUFFLE(1,0,3,2));//
		 		
	 __m128 PathMetricsB = _mm_shuffle_ps(PathMetrics,PathMetrics,_MM_SHUFFLE(3,1,3,1));

	 __m128 SumA = _mm_add_ps(PathMetricsA,BranchMetricsA);
	 __m128 SumB = _mm_add_ps(PathMetricsB,BranchMetricsB);

	__m128 NewPathMetrics = _mm_min_ps(SumA,SumB);
	__m128 Decisions4 = _mm_cmpeq_ps(SumB,NewPathMetrics);

	//Store the Path Metrics

	_mm_storel_pi((__m64*) &vecPathMetMem[intNewState+intNewPathMetMemPtr], NewPathMetrics );
	_mm_storeh_pi((__m64*) &vecPathMetMem[intNewState+intNewPathMetMemPtr+(intNoStates>>1)], NewPathMetrics );


	//Store the Decisions

	_mm_storel_pi((__m64*) (intDecisions+intNewState), Decisions4 );
	_mm_storeh_pi((__m64*) (intDecisions+intNewState+(intNoStates>>1)), Decisions4 );

}


/**
 * @brief Soft-input decode of @a InputLength QPSK symbols; writes bits to @a Output; updates @a MetricsGrowth.
 */
void Viterbi::Decode(float *InputI, float *InputQ, unsigned int InputLength, unsigned char *Output, bool ExchangeIQ, float SignI, float SignQ, float &MetricsGrowth)
{
	
	int Ctr = 0;
	
	bvecDecodedBits = Output;
	float BestMetric1;
	float *pInputI = InputI;
	float *pInputQ = InputQ;
	if(ExchangeIQ)
	{
		pInputI = InputQ;
		pInputQ = InputI;
	}

	__m256 mSignI = _mm256_set1_ps(SignI);
	__m256 mSignQ = _mm256_set1_ps(SignQ);
	for (unsigned int ii = 0; ii < InputLength; ii += 8)
	{		
		CalcMetrics2(pInputI+ii,pInputQ+ii,mSignI,mSignQ);
		float *pBranchMet = BranchMetrics;

		for( int jj = 0; jj < 8; jj++)
		{
			AcsAllNew(pBranchMet, bmatSurvMem + intSurvMemWrRowAddr);
			intNewPathMetMemPtr = intOldPathMetMemPtr;
			intOldPathMetMemPtr = intNoStates - intOldPathMetMemPtr;
			intSurvMemWrRowAddr = (intSurvMemWrRowAddr+1)&(MaskSurv);
			pBranchMet += 8;
		}
	}

	int BestIndex;
	float BestMetric2 = CalcBestMetric(vecPathMetMem + intOldPathMetMemPtr, BestIndex);
	//Normalize
	if(BestMetric2 > 0)
	{
		__m256 mBest = _mm256_set1_ps(BestMetric2);
		for(int i = 0; i < 64; i+=8)
		{
			__m256 mIn = _mm256_load_ps(vecPathMetMem + intOldPathMetMemPtr+i);
			mIn = _mm256_sub_ps(mIn,mBest);
			_mm256_store_ps(vecPathMetMem + intOldPathMetMemPtr+i,mIn);
		}
	}

	traceback_mid(Output, BestIndex, InputLength);

	MetricsGrowth = BestMetric2;
	NumBits += InputLength;
	//MetricsGrowth = MetricsGrowth / float(Ctr - intTrcbckLen);

}

void  Viterbi::CalcMetrics2(float *InputI, float *InputQ, __m256 SignI, __m256 SignQ)
{
	
	__m256 A = _mm256_loadu_ps(InputI);
	__m256 B = _mm256_loadu_ps(InputQ);
	A = _mm256_mul_ps(A, SignI);
	B = _mm256_mul_ps(B, SignQ);
	
	__m256 mZero = _mm256_setzero_ps();


	__m256 mMetANeg = _mm256_max_ps(mZero, A);// Every negative is 0
	__m256 mMetBNeg = _mm256_max_ps(mZero, B);// Every negative is 0
	__m256 mMetBPos = _mm256_sub_ps(mMetBNeg, B);//Every positive is 0;
	__m256 mMetAPos = _mm256_sub_ps(mMetANeg, A);//Every positive is 0;
	

	__m256 M00 = _mm256_add_ps(mMetANeg,mMetBNeg);
	__m256 M01 = _mm256_add_ps(mMetANeg, mMetBPos);
	__m256 M11 = _mm256_add_ps(mMetAPos, mMetBPos);
	__m256 M10 = _mm256_add_ps(mMetAPos, mMetBNeg);

	//M00 = _mm256_setr_ps(0,1,2,3,4,5,6,7);
	//M01 = _mm256_setr_ps(10,11,12,13,14,15,16,17);
	//M10 = _mm256_setr_ps(20,21,22,23,24,25,26,27);
	//M11 = _mm256_setr_ps(30,31,32,33,34,35,36,37);

	__m256 ABLow = _mm256_unpacklo_ps(M00,M01);
	__m256 CDLow = _mm256_unpacklo_ps(M10,M11);

	__m256 X01 = _mm256_permute2f128_ps(ABLow,CDLow,0x20);//A0B0A1B1C0D0C1D1
	__m256 X45 = _mm256_permute2f128_ps(ABLow,CDLow,0x31);//A4B4A5B5C4D4C5D5
	__m256 Y0 = _mm256_castsi256_ps(_mm256_permute4x64_epi64(_mm256_castps_si256(X01),0x88));
	__m256 Y1 = _mm256_castsi256_ps(_mm256_permute4x64_epi64(_mm256_castps_si256(X01),0xDD));
	__m256 Y4 = _mm256_castsi256_ps(_mm256_permute4x64_epi64(_mm256_castps_si256(X45),0x88));
	__m256 Y5 = _mm256_castsi256_ps(_mm256_permute4x64_epi64(_mm256_castps_si256(X45),0xDD));
	
	_mm256_store_ps(BranchMetrics, Y0);
	_mm256_store_ps(BranchMetrics+8, Y1);
	_mm256_store_ps(BranchMetrics+32, Y4);
	_mm256_store_ps(BranchMetrics+40, Y5);
	
	
	
	__m256 ABHigh = _mm256_unpackhi_ps(M00,M01);//A2B2A3B3....
	__m256 CDHigh = _mm256_unpackhi_ps(M10,M11);

	__m256 X23 = _mm256_permute2f128_ps(ABHigh,CDHigh,0x20);//A2B2....
	__m256 X67 = _mm256_permute2f128_ps(ABHigh,CDHigh,0x31);
	__m256 Y2 = _mm256_castsi256_ps(_mm256_permute4x64_epi64(_mm256_castps_si256(X23),0x88));
	__m256 Y3 = _mm256_castsi256_ps(_mm256_permute4x64_epi64(_mm256_castps_si256(X23),0xDD));
	__m256 Y6 = _mm256_castsi256_ps(_mm256_permute4x64_epi64(_mm256_castps_si256(X67),0x88));
	__m256 Y7 = _mm256_castsi256_ps(_mm256_permute4x64_epi64(_mm256_castps_si256(X67),0xDD));
	_mm256_store_ps(BranchMetrics+16, Y2);
	_mm256_store_ps(BranchMetrics+24, Y3);
	_mm256_store_ps(BranchMetrics+48, Y6);
	_mm256_store_ps(BranchMetrics+56, Y7);
	
	

}

void Viterbi::CalcMetrics2(float *Input, unsigned int InputLength)
{
	

	
	__m256 *pInput = (__m256 *) Input;

	__m256 *pBranchMetrics = (__m256 *) BranchMetrics;

	__m256 mZero = _mm256_setzero_ps();

	for(unsigned ii =0; ii < InputLength; ii+=16)
	{
	

	

		__m256 mInputA = _mm256_load_ps((float*)pInput);
		pInput++; //Advance the input
		__m256 mInputB = _mm256_load_ps((float*) pInput);
		pInput++; //Advance the input

		__m256 X0_8 = _mm256_permute2f128_ps(mInputA,mInputB,0x20);//0123,8,9,10,11
		__m256 X4_12 = _mm256_permute2f128_ps(mInputA,mInputB,0x31);//4,5,6,7,12,13,14,15
		__m256 A = _mm256_shuffle_ps(X0_8,X4_12,0x88);//0 2 4 6,8,10,12,14
		__m256 B = _mm256_shuffle_ps(X0_8,X4_12,0xDD);//1 3 5 7, 9,11,13,15




		__m256 mMetANeg = _mm256_max_ps(mZero, A);// Every negative is 0
		__m256 mMetBNeg = _mm256_max_ps(mZero, B);// Every negative is 0
		__m256 mMetBPos = _mm256_sub_ps(mMetBNeg, B);//Every positive is 0;
		__m256 mMetAPos = _mm256_sub_ps(mMetANeg, A);//Every positive is 0;
		

		__m256 M00 = _mm256_add_ps(mMetANeg,mMetBNeg);
		__m256 M01 = _mm256_add_ps(mMetANeg, mMetBPos);
		__m256 M11 = _mm256_add_ps(mMetAPos, mMetBPos);
		__m256 M10 = _mm256_add_ps(mMetAPos, mMetBNeg);








		__m256 ABLow = _mm256_unpacklo_ps(M00,M01);
		__m256 CDLow = _mm256_unpacklo_ps(M10,M11);

		__m256 X01 = _mm256_permute2f128_ps(ABLow,CDLow,0x20);//A0B0A1B1C0D0C1D1
		__m256 X45 = _mm256_permute2f128_ps(ABLow,CDLow,0x31);//A4B4A5B5C4D4C5D5

		_mm256_store_ps((float *) pBranchMetrics, X01);
		pBranchMetrics++;
	
		__m256 ABHigh = _mm256_unpackhi_ps(M00,M01);//A2B2A3B3....
		__m256 CDHigh = _mm256_unpackhi_ps(M10,M11);

		__m256 X23 = _mm256_permute2f128_ps(ABHigh,CDHigh,0x20);//A2B2....
		__m256 X67 = _mm256_permute2f128_ps(ABHigh,CDHigh,0x31);
		_mm256_store_ps((float *) pBranchMetrics, X23);
		pBranchMetrics++;


		_mm256_store_ps((float *) pBranchMetrics, X45);
		pBranchMetrics++;
		_mm256_store_ps((float *) pBranchMetrics, X67);
		pBranchMetrics++;

		//A0B0A1B1C0D0C1D1,A2B2A3B3C2D2C3D3,


	}

}


void Viterbi::AcsAllNew(float *BranchMet, uint64_t *Decisions)
{

	//for(int i = 0; i < 64; i++)
	//vecPathMetMem[intOldPathMetMemPtr + i] = i;
	
	alignas (32) unsigned int TempDecisions[4];
	__m256 Input = _mm256_load_ps(BranchMet);
	//Input = _mm256_setr_ps(0,1,2,3,0,1,2,3);

	__m256 BranchMet4[4];

	BranchMet4[0] = _mm256_permutevar_ps(Input, mArranged[0]);
	BranchMet4[1] = _mm256_permute2f128_ps(BranchMet4[0],Input,0x01);
	BranchMet4[2] = _mm256_permutevar_ps(Input, mArranged[1]);
	BranchMet4[3] = _mm256_permute2f128_ps(BranchMet4[2],Input,0x01);
	int TableA[8] = {0,1,3,2,1,0,2,3};
	int TableB[8] = {1,0,2,3,0,1,3,2};
	
	alignas(32) int LocalDecisionsLow[4], LocalDecisionsHigh[4];

	int PtrTbl = 0;
	int State = 0;
	int StateOut = 0;
	int PtrOut = 0;
	for(; State < 64; )
	{
		__m256 mInputA = _mm256_load_ps(&vecPathMetMem[intOldPathMetMemPtr + State]);
		State += 8;
		__m256 mInputB = _mm256_load_ps(&vecPathMetMem[intOldPathMetMemPtr + State]);
		State += 8;
		__m256 X0_8 = _mm256_insertf128_ps(mInputA, _mm256_castps256_ps128(mInputB), 1);
		__m256 X4_12 = _mm256_permute2f128_ps(mInputA, mInputB, 0x31);//4,5,6,7,12,13,14,15
		__m256 mEven = _mm256_shuffle_ps(X0_8, X4_12, 0x88);//0 2 4 6,8,10,12,14
		__m256 mOdd = _mm256_shuffle_ps(X0_8, X4_12, 0xDD);//1 3 5 7, 9,11,13,15
		__m256 SumA = _mm256_add_ps(mEven,BranchMet4[TableA[PtrTbl]]);
		__m256 SumB = _mm256_add_ps(mOdd,BranchMet4[TableB[PtrTbl]]);
		__m256 SumC = _mm256_add_ps(mEven,BranchMet4[TableB[PtrTbl]]);
		__m256 SumD = _mm256_add_ps(mOdd,BranchMet4[TableA[PtrTbl]]);
		PtrTbl++;
		__m256 CmpAB = _mm256_cmp_ps(SumB, SumA, _CMP_LT_OQ);
		__m256 CmpCD = _mm256_cmp_ps(SumD, SumC, _CMP_LT_OQ);
		__m256 NewMetLow = _mm256_blendv_ps(SumA, SumB, CmpAB);
		__m256 NewMetHigh = _mm256_blendv_ps(SumC, SumD, CmpCD);
		LocalDecisionsLow[PtrOut] = _mm256_movemask_ps(CmpAB);
		LocalDecisionsHigh[PtrOut] = _mm256_movemask_ps(CmpCD);
		PtrOut++;
		_mm256_store_ps(vecPathMetMem + intNewPathMetMemPtr + StateOut, NewMetLow);
		_mm256_store_ps(vecPathMetMem + intNewPathMetMemPtr + StateOut+32, NewMetHigh);
		StateOut += 8;
	}
	
	

	*Decisions = (uint64_t)(((uint64_t) LocalDecisionsLow[0] ) | (((uint64_t) LocalDecisionsLow[1] ) << 8) | (((uint64_t) LocalDecisionsLow[2] ) << 16) | (((uint64_t) LocalDecisionsLow[3] ) << 24) | (((uint64_t) LocalDecisionsHigh[0] ) << 32) | (((uint64_t) LocalDecisionsHigh[1] ) << 40) | (((uint64_t) LocalDecisionsHigh[2] ) << 48) | (((uint64_t) LocalDecisionsHigh[3] ) << 56));   
}


void Viterbi::AcsAll(float *BranchMet, uint64_t* Decisions, unsigned int Phase)
{
	alignas (32) unsigned int TempDecisions[4];

	//Create the required structures

	__m256 BranchMet4[4];

	//A0B0A1B1C0D0C1D1
	__m256 Input = _mm256_load_ps(BranchMet);
	if (Phase == 0)
	{
		BranchMet4[1] = _mm256_permutevar_ps(Input, mArranged[1]);//10102323
		BranchMet4[0] = _mm256_permutevar_ps(Input, mArranged[0]);//01013232      0100 0100
		BranchMet4[2] = _mm256_permute2f128_ps(BranchMet4[1], BranchMet4[1], 0x01);
		BranchMet4[3] = _mm256_permute2f128_ps(BranchMet4[0], BranchMet4[0], 0x01);

	}
	else
	{
		BranchMet4[1] = _mm256_permutevar_ps(Input, mArranged[3]);//10102323
		BranchMet4[0] = _mm256_permutevar_ps(Input, mArranged[2]);//01013232      0100 0100
		BranchMet4[2] = _mm256_permute2f128_ps(BranchMet4[1], BranchMet4[1], 0x01);
		BranchMet4[3] = _mm256_permute2f128_ps(BranchMet4[0], BranchMet4[0], 0x01);


	}


	alignas(32) int LocalDecisionsA[4], LocalDecisionsB[4];

	//	IACA_START
	//First 8

	__m256 mInputA = _mm256_load_ps(&vecPathMetMem[intOldPathMetMemPtr + 0]);

	__m256 mInputB = _mm256_load_ps(&vecPathMetMem[intOldPathMetMemPtr + 8]);
	//IACA_END
	//__m256 X0_8 = _mm256_permute2f128_ps(mInputA,mInputB,0x20);//0123,8,9,10,11
	__m256 X0_8 = _mm256_insertf128_ps(mInputA, _mm_load_ps(&vecPathMetMem[intOldPathMetMemPtr + 8]), 1);
	__m256 X4_12 = _mm256_permute2f128_ps(mInputA, mInputB, 0x31);//4,5,6,7,12,13,14,15

	__m256 A = _mm256_shuffle_ps(X0_8, X4_12, 0x88);//0 2 4 6,8,10,12,14
	__m256 B = _mm256_shuffle_ps(X0_8, X4_12, 0xDD);//1 3 5 7, 9,11,13,15

	__m256 SumA = _mm256_add_ps(A, BranchMet4[0]);
	__m256 SumB = _mm256_add_ps(B, BranchMet4[3]);

	mInputA = _mm256_load_ps(&vecPathMetMem[intOldPathMetMemPtr + 16]);
	mInputB = _mm256_load_ps(&vecPathMetMem[intOldPathMetMemPtr + 24]);

	__m256 CmpAB = _mm256_cmp_ps(SumB, SumA, 2);
	__m256 NewMet = _mm256_blendv_ps(SumA, SumB, CmpAB);
	LocalDecisionsA[0] = _mm256_movemask_ps(CmpAB);
	_mm256_store_ps(vecPathMetMem + intNewPathMetMemPtr, NewMet);



	X0_8 = _mm256_insertf128_ps(mInputA, _mm_load_ps(&vecPathMetMem[intOldPathMetMemPtr + 24]), 1);
	X4_12 = _mm256_permute2f128_ps(mInputA, mInputB, 0x31);//4,5,6,7,12,13,14,15


	SumA = _mm256_add_ps(A, BranchMet4[3]);
	SumB = _mm256_add_ps(B, BranchMet4[0]);

	A = _mm256_shuffle_ps(X0_8, X4_12, 0x88);//0 2 4 6,8,10,12,14
	B = _mm256_shuffle_ps(X0_8, X4_12, 0xDD);//1 3 5 7, 9,11,13,15


	__m256 CmpCD = _mm256_cmp_ps(SumB, SumA, 2);
	NewMet = _mm256_blendv_ps(SumA, SumB, CmpCD);
	LocalDecisionsB[0] = _mm256_movemask_ps(CmpCD);
	_mm256_store_ps(vecPathMetMem + intNewPathMetMemPtr + 32, NewMet);


	//Next 8 

	//X0_8 = _mm256_permute2f128_ps(mInputA,mInputB,0x20);//0123,8,9,10,11





	SumA = _mm256_add_ps(A, BranchMet4[3]);
	SumB = _mm256_add_ps(B, BranchMet4[0]);

	mInputA = _mm256_load_ps(&vecPathMetMem[intOldPathMetMemPtr + 32]);
	mInputB = _mm256_load_ps(&vecPathMetMem[intOldPathMetMemPtr + 40]);


	CmpAB = _mm256_cmp_ps(SumB, SumA, 2);
	NewMet = _mm256_blendv_ps(SumA, SumB, CmpAB);
	LocalDecisionsA[1] = _mm256_movemask_ps(CmpAB);
	_mm256_store_ps(vecPathMetMem + intNewPathMetMemPtr + 8, NewMet);

	X0_8 = _mm256_insertf128_ps(mInputA, _mm_load_ps(&vecPathMetMem[intOldPathMetMemPtr + 40]), 1);
	X4_12 = _mm256_permute2f128_ps(mInputA, mInputB, 0x31);//4,5,6,7,12,13,14,15


	SumA = _mm256_add_ps(A, BranchMet4[0]);
	SumB = _mm256_add_ps(B, BranchMet4[3]);

	A = _mm256_shuffle_ps(X0_8, X4_12, 0x88);//0 2 4 6,8,10,12,14
	B = _mm256_shuffle_ps(X0_8, X4_12, 0xDD);//1 3 5 7, 9,11,13,15

	CmpCD = _mm256_cmp_ps(SumB, SumA, 2);
	NewMet = _mm256_blendv_ps(SumA, SumB, CmpCD);
	LocalDecisionsB[1] = _mm256_movemask_ps(CmpCD);


	_mm256_store_ps(vecPathMetMem + intNewPathMetMemPtr + 40, NewMet);



	//Next 8 

	//X0_8 = _mm256_permute2f128_ps(mInputA,mInputB,0x20);//0123,8,9,10,11




	SumA = _mm256_add_ps(A, BranchMet4[2]);
	SumB = _mm256_add_ps(B, BranchMet4[1]);

	mInputA = _mm256_load_ps(&vecPathMetMem[intOldPathMetMemPtr + 48]);
	mInputB = _mm256_load_ps(&vecPathMetMem[intOldPathMetMemPtr + 56]);


	CmpAB = _mm256_cmp_ps(SumB, SumA, 2);
	NewMet = _mm256_blendv_ps(SumA, SumB, CmpAB);

	X0_8 = _mm256_insertf128_ps(mInputA, _mm_load_ps(&vecPathMetMem[intOldPathMetMemPtr + 56]), 1);
	X4_12 = _mm256_permute2f128_ps(mInputA, mInputB, 0x31);//4,5,6,7,12,13,14,15


	LocalDecisionsA[2] = _mm256_movemask_ps(CmpAB);
	_mm256_store_ps(vecPathMetMem + intNewPathMetMemPtr + 16, NewMet);

	SumA = _mm256_add_ps(A, BranchMet4[1]);
	SumB = _mm256_add_ps(B, BranchMet4[2]);

	A = _mm256_shuffle_ps(X0_8, X4_12, 0x88);//0 2 4 6,8,10,12,14
	B = _mm256_shuffle_ps(X0_8, X4_12, 0xDD);//1 3 5 7, 9,11,13,15

	CmpCD = _mm256_cmp_ps(SumB, SumA, 2);
	NewMet = _mm256_blendv_ps(SumA, SumB, CmpCD);
	LocalDecisionsB[2] = _mm256_movemask_ps(CmpCD);


	_mm256_store_ps(vecPathMetMem + intNewPathMetMemPtr + 48, NewMet);



	//Next 8 

	//X0_8 = _mm256_permute2f128_ps(mInputA,mInputB,0x20);//0123,8,9,10,11




	SumA = _mm256_add_ps(A, BranchMet4[1]);
	SumB = _mm256_add_ps(B, BranchMet4[2]);
	CmpAB = _mm256_cmp_ps(SumB,SumA,2);
	NewMet = _mm256_blendv_ps(SumA,SumB,CmpAB);
	LocalDecisionsA[3] = _mm256_movemask_ps(CmpAB);
	_mm256_store_ps(vecPathMetMem+intNewPathMetMemPtr+24,NewMet);

	SumA = _mm256_add_ps(A,BranchMet4[2]);
	SumB = _mm256_add_ps(B,BranchMet4[1]);
	CmpCD = _mm256_cmp_ps(SumB,SumA,2);
	NewMet = _mm256_blendv_ps(SumA,SumB,CmpCD);
	LocalDecisionsB[3] = _mm256_movemask_ps(CmpCD);
	_mm256_store_ps(vecPathMetMem+intNewPathMetMemPtr+56,NewMet);

	__m128i mDecisionsA = _mm_load_si128((__m128i*) LocalDecisionsA);
	__m128i mDecisionsB = _mm_load_si128((__m128i*) LocalDecisionsB);


	mDecisionsA = _mm_packus_epi32(mDecisionsA, mDecisionsB); //cvt to 16 bit
	mDecisionsA = _mm_packus_epi16(mDecisionsA, mDecisionsA); //cvt to 8 bit
	*Decisions = _mm_extract_epi64(mDecisionsA, 0);


	//	IACA_END
}


float Viterbi::CalcBestMetric(float *Input, int& Position)
{
	__m256i currentindex = _mm256_setr_epi32(0,1,2,3,4,5,6,7);
	__m256i bestindex = currentindex;
	
	__m256 BestMetrics = _mm256_load_ps(Input);
	__m256i Delta = _mm256_set1_epi32(8);

	for(int i= 8; i < 64; i+= 8)
	{
		__m256 New = _mm256_load_ps(Input+i);
		currentindex = _mm256_add_epi32(currentindex,Delta);
		__m256 mask = _mm256_cmp_ps(New, BestMetrics, _CMP_LT_OQ);
		BestMetrics = _mm256_blendv_ps(BestMetrics, New, mask);
		bestindex = _mm256_blendv_epi8(bestindex, currentindex, _mm256_castps_si256(mask));
	}

	__m128 BestMetrics0 = _mm256_castps256_ps128(BestMetrics);
	__m128i bestindex0 = _mm256_castsi256_si128(bestindex);
	__m128 BestMetrics1 = _mm256_extractf128_ps(BestMetrics, 1);
	__m128i bestindex1 = _mm256_extracti128_si256(bestindex, 1);
	__m128 mask = _mm_cmp_ps(BestMetrics1, BestMetrics0, _CMP_LT_OQ);
	BestMetrics0 = _mm_blendv_ps(BestMetrics0, BestMetrics1, mask);
	bestindex0 = _mm_blendv_epi8(bestindex0, bestindex1, _mm_castps_si128(mask));
	BestMetrics1 = _mm_permute_ps (BestMetrics0, 0x1B);// 00011011
	bestindex1 = _mm_shuffle_epi32(bestindex0, 0x1B);
	mask = _mm_cmp_ps(BestMetrics1, BestMetrics0, _CMP_LT_OQ);
	BestMetrics0 = _mm_blendv_ps(BestMetrics0, BestMetrics1, mask);
	bestindex0 = _mm_blendv_epi8(bestindex0, bestindex1, _mm_castps_si128(mask));
	
	alignas(32) float TwoMetrics[2];
	alignas(32) int TwoIndexes[2];
	_mm_store_sd((double*)TwoMetrics, _mm_castps_pd(BestMetrics0));
	_mm_storel_epi64((__m128i*)TwoIndexes, bestindex0);
	Position = TwoIndexes[0];
	float BestAll = TwoMetrics[0];
	if(TwoMetrics[1] < BestAll)
	{
		BestAll = TwoMetrics[1];
		Position = TwoIndexes[1];
	}
	return BestAll;
}


float Viterbi::CalcBestMetric(float *Input)
{

	__m256 mBestMet8 = _mm256_load_ps(Input);
	for (unsigned int ii = 8; ii < intNoStates; ii += 8)
	{
		__m256 mNew = _mm256_load_ps(Input + ii);
		mBestMet8 = _mm256_min_ps(mBestMet8, mNew);
	}

	__m256 mBestMet8_2 = _mm256_insertf128_ps(mBestMet8, _mm256_extractf128_ps(mBestMet8, 1), 0);
	 mBestMet8 = _mm256_min_ps(mBestMet8, mBestMet8_2);
	__m128 mBestMet4 = _mm256_extractf128_ps(mBestMet8, 0);
	__m128 mBestMet4_2 = _mm_permute_ps(mBestMet4, 0x4E);
	__m128 mBestMet2 = _mm_min_ps(mBestMet4, mBestMet4_2);
	__m128 mBestMet2_2 = _mm_permute_ps(mBestMet2, 0x11);
	__m128 mBestMet = _mm_min_ps(mBestMet2, mBestMet2_2);
	
	float Result;

	_mm_store_ss(&Result, mBestMet);

	return Result;
}




void Viterbi::reset_decoder(void)
{
	//reset the number of decoded bits
	NumBits = 0;

	//reset the Survivor Row Address Pointers
	intSurvMemWrRowAddr = 0;
	intSurvMemRdRowAddr = 0;
	//reset the Old and New Path Metrics pointers
	intOldPathMetMemPtr = 0;
	intNewPathMetMemPtr = intNoStates; 

	// reset the Path Metrics all 0  0 the preffered one
	for(unsigned int ii = 0; ii < 2*intNoStates; ii++)
	{
		vecPathMetMem[ii] = 0;
	}

	for(int i = MaskSurv; i >= MAX_BURST_LENGTH-256; i--)
		bmatSurvMem[i] = 0;
	
}