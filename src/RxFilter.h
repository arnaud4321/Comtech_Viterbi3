/**
 * @file RxFilter.h
 * @brief Receive matched filter: same sqrt-RC @c ippsFIRSR structure as @ref TxFilter, gain @c K = 1/8192.
 */
#pragma once
#include "TxFilter.h"
class RxFilter: public TxFilter
{
    const double K = 1.0/8192.0;
private:
    /* data */
public:
    RxFilter(void);
    ~RxFilter();
    void CreateObjects(Ipp32f RollofIn);
    void CreateOutputs(Ipp32f *InputI, Ipp32f *InputQ, Ipp32f *OutI, Ipp32f *OutQ,  unsigned int InputLength);

};


