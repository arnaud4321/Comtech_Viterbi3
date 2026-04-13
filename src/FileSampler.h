/**
 * @file FileSampler.h
 * @brief Replay interleaved int16 IQ from disk into the same @c BufferShort ring as @ref Sampler (no USRP).
 */
#pragma once

#include <string>
#include <thread>
#include "Sampler.h"

/**
 * @brief Feeds the RX matched-filter thread from a binary IQ file (short I,Q interleaved, native endian).
 *
 * @details Reads @c 2*SPB shorts per batch (same layout as @ref UHDSampler). Optional real-time pacing matches
 * @ref Sampler. On EOF: if @a loop is true, rewinds the file; otherwise stops the producer thread after padding
 * the last partial batch with zeros.
 */
class FileSampler : public Sampler
{
public:
    /**
     * @param iqFilePath Path to binary file (e.g. @c InputRecording.bin from USRP capture).
     * @param loopPlayback If true, rewind at EOF; if false, exit thread after EOF.
     * @param debugIn Same meaning as @ref Sampler(bool): multiplies pacing interval by 10 when true.
     */
    FileSampler(const std::string& iqFilePath, bool loopPlayback, bool debugIn);

    void StartThread() override;

private:
    std::string path_;
    bool loop_{true};

    void OperateFileSampler();
};
