/**
 * @file USRPInit.h
 * @brief USRP hardware initialization via the UHD API.
 */
#pragma once

#include <uhd/usrp/multi_usrp.hpp>
#include <string>
#include "definitions.h"
using namespace std;
typedef std::function<uhd::sensor_value_t(const std::string &)> get_sensor_fn_t;

/**
 * @brief Handles the initial configuration and setup of a USRP device.
 * 
 * @details This class is responsible for opening the device via `uhd::usrp::multi_usrp::make()`,
 * configuring sampling rates, center frequencies, gains, and checking PLL lock status
 * before the TX and RX DSP pipelines are started.
 */
class USRPInit
{
private:
    /* data */
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

    while (true)
    {
        if (lock_detected and (std::chrono::steady_clock::now() > setup_timeout))
        {
            std::cout << " locked." << std::endl;
            break;
        }
        if (get_sensor_fn(sensor_name).to_bool())
        {
            std::cout << "+";
            std::cout.flush();
            lock_detected = true;
        }
        else
        {
            if (std::chrono::steady_clock::now() > setup_timeout)
            {
                std::cout << std::endl;
                throw std::runtime_error(
                    str(boost::format(
                            "timed out waiting for consecutive locks on sensor \"%s\"") %
                        sensor_name));
            }
            std::cout << "_";
            std::cout.flush();
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    std::cout << std::endl;
    return true;
};
public:
    uhd::usrp::multi_usrp::sptr usrp;

    USRPInit(double RxFreq, double TxFreq, double TxBW, double TxGain, double lo_offset, string ref, double rxRate, double txRate);
    ~USRPInit();
};

