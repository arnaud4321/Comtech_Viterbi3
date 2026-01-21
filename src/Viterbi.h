#pragma once
#include <immintrin.h>
#include <cstdint>

#define MAX_BURST_LENGTH 65536
#define AVX2_CODE
class Viterbi
{
	


protected:
	float *BranchMetrics;
	__m128i *mDecsMask;
	__m256i *mArranged;
	__m128i *mPermCalc;
	__m256i *mOffset;
	//Branch Metrics Calculation

	int NumTransitions;
	int *AllBinVecs; // Auxiliary matrix containing all n-tuples vector
	int CLengthM1, FirstTracebackPeriod;
	//Depuncturing

	//unsigned int * bvecPuncPatt; //Puncturing Pattern
	//unsigned int PuncPattCycle; // The length of the cycle of the puncturing pattern
	unsigned int MaxOutLen; //Maximal Expected Output Length
	unsigned int OutputCtr; //Counter of the length of the output
	//float *DepuncOutputBuffer;

	unsigned int FirstPositionPunc, *DiffPunc, DiffPuncLength;


	//Private Methods

	void generate_binvecs(void);
	void dec2bin(int index, int *v, int length);

	unsigned int *ivecGenPolys; 
	unsigned int intCLength; //Constraint Length
	unsigned int intNParam; //the n parameter of 1/n
	unsigned int intMaxDataLen; // The Maximal Output Decoded Data Length (does not include the tail bits)
	unsigned int intTrcbckLen;// The maximal possible traceback length
	unsigned int intNoStates ; //The Number of States of the Decoder
	unsigned int intMask; 
	unsigned int intOldPathMetMemPtr; // Pointers to the old and new banks of Path Metrics
	unsigned int intNewPathMetMemPtr;
	unsigned int intSurvMemRdRowAddr; //Pointer to Rd and Wr of the Row Address to the Survivor Memory
	unsigned int intSurvMemWrRowAddr;
	unsigned int MaskSurv;
	unsigned int intNumDecodedBits; //Number of Bits that had been decoded 
	float *vecPathMetMem, *TempVec; //Path Metrics Memory Vector
	__m128i PathMet[8];
	uint64_t *bmatSurvMem;  //Survivor Memory Matrix
	int *imatDecodingTable;
	unsigned char *bvecDecodedBits;

	//Private Functions
	void calc_decoding_table(void);
	int is_odd_ones(int input);
	void add_comp_select(int intNewState, float *CurrBranchMetrics,  int *intDecisions); // add compare select implementation
	void AcsAll(float *BranchMet, uint64_t *Decisions, unsigned int Phase);
	void AcsAllNew(float *BranchMet, uint64_t *Decisions);

	//void AcsAll(unsigned int BranchMet, __int64 *Decisions);
	void traceback_mid( unsigned char *Out, int BestMetInd, int Length);
	void CalcMetrics2(float *Input, unsigned int InputLength);


	float CalcBestMetric(float *Input);
	float CalcBestMetric(float *Input, int &Position);
	void CalcMetrics2(float *InputI, float *InputQ, __m256 SignI, __m256 SignQ);

	bool ReversedPoly;
	int DecoderID;
	uint64_t NumBits = 0;
	public:
	Viterbi(int Idx);
	//GenPolysInp - vector containing the Generator Polynomials
	//MaxDataLen - The maximal possible Data Length not including the tail bits
	//TrcbckLenInp - the Traceback Length. If 0 - Traceback is performed once at end of packet
	 ~Viterbi(void);
		//matBranchMetrics length data x NParam
	//Traceback Length - can overide the existing Traceback length, however does not allow to change
	//Zero to Non-Zero Length and otherwise
	

	void Decode(float *InputI, float *InputQ, unsigned int InputLength, unsigned char *Output, bool ExchangeIQ, float SignI, float SignQ, float &MetricsGrowth);
	void reset_decoder(void);


};



