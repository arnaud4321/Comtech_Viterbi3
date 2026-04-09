
#include <boost/algorithm/string.hpp>

#include "USRPInit.h"
#include "definitions.h"
USRPInit::USRPInit(double RxFreq, double TxFreq, double TxBW, double TxGain, double lo_offset, string ref, double rxRate, double txRate)
{
    std::string args = "num_recv_frames=256";
    std::cout << std::endl;
    std::cout << boost::format("Creating the usrp device with: %s...") % args
              << std::endl;
    usrp = uhd::usrp::multi_usrp::make(args);
    double rate = rxRate;
    string ant = "RX2";

    string wirefmt = "sc12";
    double setup_time = 1.0;
    std::cout << std::endl;
    std::cout << boost::format("Creating the usrp device with: %s...") % args
              << std::endl;
    usrp = uhd::usrp::multi_usrp::make(args);
    int n = usrp->get_rx_num_channels();
    std::cout << boost::format("Lock mboard clocks: %f") % ref << std::endl;
    usrp->set_clock_source(ref);
    // always select the subdevice first, the channel mapping affects the other settings

    std::vector<size_t> channel_list;
    std::vector<std::string> channel_strings;
    std::string channels = "0";
    boost::split(channel_strings, channels, boost::is_any_of("\"',"));
    for (size_t ch = 0; ch < channel_strings.size(); ch++) {
        try {
            int chan = std::stoi(channel_strings[ch]);
            if (chan >= static_cast<int>(usrp->get_rx_num_channels()) || chan < 0) {
                throw std::runtime_error("Invalid channel(s) specified.");
            } else {
                channel_list.push_back(static_cast<size_t>(chan));
            }
        } catch (std::invalid_argument const& c) {
            throw std::runtime_error("Invalid channel(s) specified.");
        } catch (std::out_of_range const& c) {
            throw std::runtime_error("Invalid channel(s) specified.");
        }
    }
  // set the sample rate
    std::cout << boost::format("Setting RX Rate: %f Msps...") % (rate / 1e6) << std::endl;
    usrp->set_rx_rate(rate);
   double ActualRate = usrp->get_rx_rate() / 1e6;
    std::cout << boost::format("Actual RX Rate: %f Msps...") % (ActualRate)
              << std::endl
              << std::endl;

    // with default of 0.0 this will always be true
    std::cout << boost::format("Setting RX Freq: %f MHz...") % (RxFreq / 1e6)
              << std::endl;
    std::cout << boost::format("Setting RX LO Offset: %f MHz...") % (lo_offset / 1e6)
              << std::endl;
    uhd::tune_request_t tune_request(RxFreq, lo_offset);
    //       if (vm.count("int-n"))
    //           tune_request.args = uhd::device_addr_t("mode_n=integer");
    usrp->set_rx_freq(tune_request);
    std::cout << boost::format("Actual RX Freq: %f MHz...") % (usrp->get_rx_freq(0) / 1e6)
              << std::endl
              << std::endl;
    double gain = RX_INITIAL_GAIN;
    std::cout << boost::format("Setting RX Gain: %f dB...") % gain << std::endl;
    usrp->set_rx_gain(gain, 0);
    gain = usrp->get_rx_gain(0);
    std::cout << boost::format("Actual RX Gain: %f dB...") % gain
              << std::endl
              << std::endl;
    //std::cout << boost::format("Setting RX Bandwidth: %f MHz...") % (bw / 1e6)
     //         << std::endl;
    //usrp->set_rx_bandwidth(bw, 0);
    //std::cout << boost::format("Actual RX Bandwidth: %f MHz...") % (usrp->get_rx_bandwidth(0) / 1e6)
     //         << std::endl
      //        << std::endl;
    std::cout << boost::format("Setting RX Antenna: %s") % ant << std::endl;
    usrp->set_rx_antenna(ant, 0);


    std::string subdev("A:A");

    usrp->set_tx_subdev_spec(subdev);
    rate = txRate;
 // set sample rate
    std::cout << boost::format("Setting TX Rate: %f Msps...") % (rate / 1e6) << std::endl;
    usrp->set_tx_rate(rate);
    std::cout << boost::format("Actual TX Rate: %f Msps...") % (usrp->get_tx_rate() / 1e6) << std::endl << std::endl;

    // set freq
    std::cout << boost::format("Setting TX Freq: %f MHz...") % (TxFreq / 1e6) << std::endl;

    tune_request = uhd::tune_request_t(TxFreq, lo_offset);
    usrp->set_tx_freq(tune_request);
    std::cout << boost::format("Actual TX Freq: %f MHz...") % (usrp->get_tx_freq() / 1e6) << std::endl << std::endl;

    // set the rf gain
    std::cout << boost::format("Setting TX Gain: %f dB...") % TxGain << std::endl;
    usrp->set_tx_gain(TxGain);
    std::cout << boost::format("Actual TX Gain: %f dB...") % usrp->get_tx_gain() << std::endl << std::endl;

    // set the IF filter bandwidth
    if(TxBW > 1e6)
    {
        std::cout << boost::format("Setting TX Bandwidth: %f MHz...") % (TxBW / 1e6) << std::endl;
        usrp->set_tx_bandwidth(TxBW);
        std::cout << boost::format("Actual TX Bandwidth: %f MHz...") % (usrp->get_tx_bandwidth() / 1e6) << std::endl << std::endl;
    }
    ant = "TX/RX";

  	std::cout << boost::format("Setting TX Antenna: %s") % ant << std::endl;
    usrp->set_tx_antenna(ant);

    std::this_thread::sleep_for(std::chrono::milliseconds(int64_t(1000 * setup_time)));

    if ( ref == "external" ) {
            check_locked_sensor(
                usrp->get_mboard_sensor_names(0),
                "ref_locked",
                [&](const std::string& sensor_name) {
                    return usrp->get_mboard_sensor(sensor_name);
                },
                setup_time);
        }

}

USRPInit::~USRPInit()
{
}
