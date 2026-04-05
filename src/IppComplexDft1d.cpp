/**
 * @file IppComplexDft1d.cpp
 */

#include "IppComplexDft1d.h"

#include <cassert>
#include <cstdio>

IppComplexDft1d::~IppComplexDft1d()
{
    release();
}

IppComplexDft1d::IppComplexDft1d(IppComplexDft1d&& o) noexcept
    : length_(o.length_), pSpec_(o.pSpec_), pWork_(o.pWork_)
{
    o.length_ = 0;
    o.pSpec_ = nullptr;
    o.pWork_ = nullptr;
}

IppComplexDft1d& IppComplexDft1d::operator=(IppComplexDft1d&& o) noexcept
{
    if (this != &o)
    {
        release();
        length_ = o.length_;
        pSpec_ = o.pSpec_;
        pWork_ = o.pWork_;
        o.length_ = 0;
        o.pSpec_ = nullptr;
        o.pWork_ = nullptr;
    }
    return *this;
}

void IppComplexDft1d::release()
{
    if (pWork_)
    {
        ippsFree(pWork_);
        pWork_ = nullptr;
    }
    if (pSpec_)
    {
        ippsFree(pSpec_);
        pSpec_ = nullptr;
    }
    length_ = 0;
}

void IppComplexDft1d::ensureLength(int length)
{
    if (length <= 0)
    {
        release();
        return;
    }
    if (length == length_ && pSpec_ && pWork_)
        return;

    release();
    length_ = length;

    int sizeSpec = 0;
    int sizeInit = 0;
    int sizeBuf = 0;
    IppStatus st =
        ippsDFTGetSize_C_32fc(length_, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast, &sizeSpec, &sizeInit, &sizeBuf);
    assert(st == ippStsNoErr && sizeSpec > 0);
    (void)st;

    Ipp8u* pInit = ippsMalloc_8u(sizeInit);
    pSpec_ = ippsMalloc_8u(sizeSpec);
    pWork_ = ippsMalloc_8u(sizeBuf);
    assert(pSpec_ && pWork_ && pInit);

    st = ippsDFTInit_C_32fc(length_, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                            reinterpret_cast<IppsDFTSpec_C_32fc*>(pSpec_), pInit);
    ippsFree(pInit);
    assert(st == ippStsNoErr);
    (void)st;
}

void IppComplexDft1d::forward(const Ipp32fc* src, Ipp32fc* dst) const
{
    assert(valid() && src && dst);
    const IppStatus st = ippsDFTFwd_CToC_32fc(src, dst, reinterpret_cast<const IppsDFTSpec_C_32fc*>(pSpec_), pWork_);
    assert(st == ippStsNoErr);
    (void)st;
}
