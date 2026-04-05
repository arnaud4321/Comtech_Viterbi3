/**
 * @file IppComplexDft1d.h
 * @brief RAII forward complex DFT (float) via Intel IPP @c ippsDFTFwd_CToC_32fc (any length @f$\ge 1@f$).
 */
#pragma once

#include <ipps.h>

#include <cstddef>

/**
 * @brief One plan: length @f$N@f$, forward @f$\mathbb{C}^N \to \mathbb{C}^N@f$ with @c IPP_FFT_NODIV_BY_ANY
 *        (same scaling convention as typical unnormalized “forward FFT” usage).
 */
class IppComplexDft1d
{
public:
    IppComplexDft1d() = default;
    ~IppComplexDft1d();

    IppComplexDft1d(const IppComplexDft1d&) = delete;
    IppComplexDft1d& operator=(const IppComplexDft1d&) = delete;
    IppComplexDft1d(IppComplexDft1d&& o) noexcept;
    IppComplexDft1d& operator=(IppComplexDft1d&& o) noexcept;

    /// (Re)allocate spec/work buffers if @a length changed; no-op if already @a length.
    void ensureLength(int length);

    int length() const { return length_; }

    /// Forward DFT; @a src and @a dst must hold @c length() complex samples.
    void forward(const Ipp32fc* src, Ipp32fc* dst) const;

    bool valid() const { return length_ > 0 && pSpec_ != nullptr && pWork_ != nullptr; }

private:
    void release();

    int length_ = 0;
    Ipp8u* pSpec_ = nullptr;
    Ipp8u* pWork_ = nullptr;
};
