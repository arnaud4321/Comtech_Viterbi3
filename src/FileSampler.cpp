/**
 * @file FileSampler.cpp
 * @brief Disk IQ replay into the RX short ring (see @ref FileSampler).
 */

#include "FileSampler.h"
#include "ConsoleAlert.h"
#include "definitions.h"
#include <chrono>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <thread>

extern std::atomic<bool> Finish;

FileSampler::FileSampler(const std::string& iqFilePath, bool loopPlayback, bool debugIn)
    : Sampler(debugIn), path_(iqFilePath), loop_(loopPlayback)
{
}

void FileSampler::StartThread()
{
    oBuffer.Reset();
    StopAll = false;
    SamplingThread = std::thread(&FileSampler::OperateFileSampler, this);
}

void FileSampler::OperateFileSampler()
{
    FILE* fid = std::fopen(path_.c_str(), "rb");
    if (!fid) {
        std::cerr << "[FileSampler] ERROR: cannot open IQ file \"" << path_ << "\"" << std::endl;
        Finish.store(true, std::memory_order_relaxed);
        return;
    }

    NumBatches = 0;
    auto Start = std::chrono::high_resolution_clock::now();

    auto ThroughputWindowStart = std::chrono::steady_clock::now();
    auto SamplerConsoleWindowStart = std::chrono::steady_clock::now();
    uint64_t SamplesInWindow = 0;
    auto lastSamplerRingSatLog = std::chrono::steady_clock::now();

    std::cout << "[FileSampler] Reading int16 IQ from \"" << path_ << "\""
              << " batchShorts=" << (2 * SPB) << " loop=" << (loop_ ? 1 : 0) << std::endl;

    while (!StopAll) {
        alignas(32) short chTmp[SPB * 2];
        const size_t toRead = static_cast<size_t>(SPB * 2);
        size_t n = std::fread(chTmp, sizeof(short), toRead, fid);

        if (n < toRead) {
            if (loop_) {
                std::clearerr(fid);
                if (std::fseek(fid, 0, SEEK_SET) != 0) {
                    std::cerr << "[FileSampler] ERROR: fseek failed, stopping." << std::endl;
                    break;
                }
                const size_t n2 = std::fread(chTmp + n, sizeof(short), toRead - n, fid);
                n += n2;
            }
            if (n < toRead) {
                if (n == 0 && !loop_) {
                    std::cout << "[FileSampler] EOF (no loop), stopping producer." << std::endl;
                    StopAll.store(true, std::memory_order_relaxed);
                    Finish.store(true, std::memory_order_relaxed);
                    CvSamplerUser.notify_all();
                    break;
                }
                std::memset(chTmp + n, 0, (toRead - n) * sizeof(short));
            }
        }

        {
            std::unique_lock<std::mutex> lk(mtxSamplerBuffer_);
            while (!StopAll && oBuffer.AlmostFull()) {
                const auto nowSat = std::chrono::steady_clock::now();
                if (nowSat - lastSamplerRingSatLog >= std::chrono::seconds(1)) {
                    lastSamplerRingSatLog = nowSat;
                    CONSOLE_ALERT_STMT(std::cout << ConsoleAlert::kRedOpen
                                                 << "[FileSampler] RX short ring almost full; fill="
                                                 << oBuffer.GetSizeInBuffer() << "/"
                                                 << oBuffer.GetBufferSize() - 1 << ConsoleAlert::kReset
                                                 << std::endl;);
                }
                CvSamplerUser.wait(lk);
            }

            short* BufOut = oBuffer.GetWriteBuffer(SPB * 2);
            std::memcpy(BufOut, chTmp, sizeof(short) * toRead);
            oBuffer.AdvancePtrWr(SPB * 2);
        }

        SamplesInWindow += SPB;

        auto nowTp = std::chrono::steady_clock::now();
        const double dtMeas = std::chrono::duration<double>(nowTp - ThroughputWindowStart).count();
        if (dtMeas >= throughputMeasurePeriodSec_) {
            const double msps = static_cast<double>(SamplesInWindow) / dtMeas / 1e6;
            lastThroughputMsps_.store(msps, std::memory_order_relaxed);
            ThroughputWindowStart = std::chrono::steady_clock::now();
            SamplesInWindow = 0;
        }
        const double dtConsole = std::chrono::duration<double>(nowTp - SamplerConsoleWindowStart).count();
        if (displayPeriodSec_ > 0.0 && dtConsole >= displayPeriodSec_) {
            std::cout << "[FileSampler] throughput=" << lastThroughputMsps_.load(std::memory_order_relaxed)
                      << " Msps" << std::endl;
            SamplerConsoleWindowStart = std::chrono::steady_clock::now();
        }

        auto Now = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = Now - Start;
        const double ExpectedTime = static_cast<double>(NumBatches) * TimeBatch;
        if (elapsed.count() < ExpectedTime) {
            std::chrono::duration<double> remaining(ExpectedTime - elapsed.count());
            std::this_thread::sleep_for(remaining);
        }
        CvSamplerUser.notify_one();
        ++NumBatches;
    }

    std::fclose(fid);
}
