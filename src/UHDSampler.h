#pragma once
#include <uhd/exception.hpp>
#include <uhd/types/tune_request.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/thread.hpp>
#include <immintrin.h>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <filesystem>
#include <string>
#include <thread>
#include "BufferShort.h"
#include "definitions.h"
#include "Sampler.h"
//#define WRITE_INPUT1
using namespace std;


typedef std::function<uhd::sensor_value_t(const std::string &)> get_sensor_fn_t;

struct SamplerParams
{
    BufferShort *pBuff;
    double Requiredgain;
    double Actualgain;
    uhd::usrp::multi_usrp::sptr usrp;
};

class UHDSampler: public Sampler
{
private:
    /* data */

    double FreqIn;
    int recmode = 0;

SamplerParams cSamplerParams;
std::thread *EttusThread;
void recv_to_buffer(uhd::usrp::multi_usrp::sptr usrp, BufferShort *pBuff);


int CalculatePower(short *Buffer, int len)
{
    long long sum = 0;
    for (int i = 0; i < len; i++)
    {
        long long val = Buffer[i];
        sum += (val * val);
    }
    long long mean = sum / len;
    return (int)(mean >> 15);
};

double OperateAGC(SamplerParams *p, double &gain);
bool check_locked_sensor(std::vector<std::string> sensor_names,
                         const char *sensor_name,
                         get_sensor_fn_t get_sensor_fn,
                         double setup_time)
{
    if (std::find(sensor_names.begin(), sensor_names.end(), sensor_name) == sensor_names.end())
        return false;

    auto setup_timeout = std::chrono::steady_clock::now() + std::chrono::milliseconds(int64_t(setup_time * 1000));
    bool lock_detected = false;

    std::cout << boost::format("Waiting for \"%s\": ") % sensor_name;
    std::cout.flush();

    while (true) {
        if (lock_detected and (std::chrono::steady_clock::now() > setup_timeout)) {
            std::cout << " locked." << std::endl;
            break;
        }
        if (get_sensor_fn(sensor_name).to_bool()) {
            std::cout << "+";
            std::cout.flush();
            lock_detected = true;
        } else {
            if (std::chrono::steady_clock::now() > setup_timeout) {
                std::cout << std::endl;
                throw std::runtime_error(
                    str(boost::format("timed out waiting for consecutive locks on sensor \"%s\"") % sensor_name));
            }
            std::cout << "_";
            std::cout.flush();
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    std::cout << std::endl;
    return true;
};

void OperateSampler(SamplerParams *p)
{
        recv_to_buffer(p->usrp, p->pBuff);
};

bool EndRecording = false;

public:
    UHDSampler(uhd::usrp::multi_usrp::sptr usrp, int recmode = 0);
    ~UHDSampler();
    void StartThread() override;
    void StopThread() override;
    void SetFrequency(double Frequency)
    {
        FreqIn = Frequency;
    }

};


