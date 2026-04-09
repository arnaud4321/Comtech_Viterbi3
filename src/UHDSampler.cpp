#include "UHDSampler.h"
#include "AWGNChannel.h"

extern std::chrono::system_clock::time_point StartOperationRx;
extern uint64_t InputTime;
extern mutex mtxfilethr;
extern bool WriteDebugParam;
UHDSampler::UHDSampler(uhd::usrp::multi_usrp::sptr usrp, int recmodex):Sampler(false)
{
    recmode = recmodex;
    std::string args = "num_recv_frames=256";
    std::string type = "short";
    size_t spb = SPB;

    cSamplerParams.pBuff = &oBuffer;
    cSamplerParams.usrp = usrp;
    cSamplerParams.Requiredgain = RX_INITIAL_GAIN;
    cSamplerParams.Actualgain = RX_INITIAL_GAIN;
}

UHDSampler::~UHDSampler()
{
}

void UHDSampler::recv_to_buffer(uhd::usrp::multi_usrp::sptr usrp, BufferShort *pBuff)
{
    mtxfilethr.lock();
    FILE *fidthr = fopen("LogThreadsInfo.txt","at");
    fprintf(fidthr,"UHD Sampler %d\n",gettid());
    fclose(fidthr);
    mtxfilethr.unlock();

    const std::string &cpu_format = "sc16";
    const std::string &wire_format = "sc12";
    const size_t &channel = 0;
    size_t samps_per_buff = SPB;
    bool continue_on_bad_packet = false;

    unsigned long long num_total_samps = 0;
    // create a receive streamer
    uhd::stream_args_t stream_args(cpu_format, wire_format);
    std::vector<size_t> channel_nums;
    channel_nums.push_back(channel);
    stream_args.channels = channel_nums;
    uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(stream_args);

    uhd::rx_metadata_t md;
    bool overflow_message = true;
    int num_requested_samples = 0;
    // setup streaming
    uhd::stream_cmd_t stream_cmd((num_requested_samples == 0)
                                     ? uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS
                                     : uhd::stream_cmd_t::STREAM_MODE_NUM_SAMPS_AND_DONE);
    stream_cmd.num_samps = size_t(num_requested_samples);
    stream_cmd.stream_now = true;
    stream_cmd.time_spec = uhd::time_spec_t();
    rx_stream->issue_stream_cmd(stream_cmd);

    // Run this loop until either time expired (if a duration was given), until
    // the requested number of samples were collected (if such a number was
    // given), or until Ctrl-C was pressed.
    bool BufferOverflow = false;
    #ifdef WRITE_INPUT1
    int CtrLook=0;
    #endif
    int Ptr1 = 0;

    std::vector<short> samples_for_recording;
    samples_for_recording.reserve(8192 * 1024);


    


    std::vector<short> tmp_buffer(2 * SPB);
    size_t samples_accumulated = 0;

    while (not StopAll) {
        size_t num_rx_samps =
            rx_stream->recv(tmp_buffer.data() + 2 * samples_accumulated, SPB - samples_accumulated, md, 0.1, false);
        
        if ( md.error_code == uhd::rx_metadata_t::ERROR_CODE_TIMEOUT ) {
            if (StopAll) break;
            // Transient timeout, just continue and try again
            continue;
        }
        if ( md.error_code == uhd::rx_metadata_t::ERROR_CODE_OVERFLOW ) {
            if ( overflow_message ) {
                overflow_message = false;
                std::cerr
                    << boost::format("Got an overflow indication.\n  This message will not appear again.\n");
            }
            continue;
        }
        if ( md.error_code != uhd::rx_metadata_t::ERROR_CODE_NONE ) {
            std::string error = str(boost::format("Receiver error: %s") % md.strerror());
            if ( continue_on_bad_packet ) {
                std::cerr << error << std::endl;
                continue;
            }
            else
                throw std::runtime_error(error);
        }

        samples_accumulated += num_rx_samps;

        if (samples_accumulated >= SPB) {
            short *Buff = pBuff->GetWriteBuffer(BufferOverflow);
            memcpy(Buff, tmp_buffer.data(), 2 * SPB * sizeof(short));

            if ( EndRecording ) {
                Ptr1 += SPB;

                int PtrWr = pBuff->GetWrPtr();
                if ( (PtrWr == 0) && (Ptr1 >= 16777216) ) {
                    FILE *fid = fopen("InputRecording.bin","wb");
                    fwrite(Buff,sizeof(short),pBuff->GetBufferSize(),fid);
                    fclose(fid);
                    cout<<"Finished writing InputRecording.bin"<<endl;
                    exit(-1);
                }
            }

            if (std::filesystem::exists("dumpit.txt"))
            {
                std::filesystem::remove("dumpit.txt");
                auto end = Buff + SPB * 2;
                samples_for_recording.insert(samples_for_recording.end(), Buff, end);
            }
            else
            {
                if (not samples_for_recording.empty())
                {
                    auto end = Buff + SPB * 2;
                    samples_for_recording.insert(samples_for_recording.end(), Buff, end);
                }

                if (samples_for_recording.size() >= 8192 * 1024)
                {
                    FILE *fid = fopen("InputRecording.bin", "wb");
                    fwrite(samples_for_recording.data(), sizeof(short), samples_for_recording.size(), fid);
                    fclose(fid);
                    cout << "Finished writing InputRecording.bin" << endl;
                    samples_for_recording.clear();
                    std::abort();
                }
            }

            {
                std::lock_guard<std::mutex> lock(mtxSamplerBuffer_); // <-- LOCK ACQUIRED HERE
                pBuff->AdvancePtrWr(2*SPB);
                CvSamplerUser.notify_one();
            }
            
            samples_accumulated = 0;
        }
    }

    stream_cmd.stream_mode = uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
    rx_stream->issue_stream_cmd(stream_cmd);
}

void UHDSampler::StartThread()
{
    cout<<"Started Sampler Thread"<<endl;
   
    EttusThread = new std::thread(&UHDSampler::OperateSampler, this, &cSamplerParams);
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    //StartOperationRx = std::chrono::system_clock::now();

    double gain = cSamplerParams.Actualgain;
    double Delta = OperateAGC(&cSamplerParams, gain);
    if ( Delta == -999999 )
        return;

    if ( Delta >= 4.0 ) {
        cout << "Warning: Signal too strong - Delta: " << Delta << "\n" << endl;
    } else if ( Delta <= -6.0 ) {
        cout << "Warning: Signal too weak - Delta: " << Delta << " (Max Gain reached)\n" << endl;
    } else
        cout << "Signal is OK - Delta: " << Delta << "\n" << endl;

    if ( recmode == 1 ) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        EndRecording = true;
    }

    cout<<"Done AGC"<<endl;
    return;
    // return 0;
}

void UHDSampler::StopThread()
{
    StopAll = true;
    CvSamplerUser.notify_all();
    if(EttusThread && EttusThread->joinable())
        EttusThread->detach();
}

double UHDSampler::OperateAGC(SamplerParams *p, double &gain)
{

    double TargetBackoff = AGCBACKOFF;                                  // db
    double TargetPower = 32768.0 * pow(10.0, -TargetBackoff * 0.1); 
    TargetPower = 10.0 * log10(TargetPower);
    cout<<"Target Backoff "<<TargetBackoff<<" Target Power "<<TargetPower<<endl;
    double Delta;
    short *Buffer;
    short *LastBuffer = 0;
    const int MaxNumSteps = 10;
    int NumSteps = 0;
    const double MaxGain = 76;
    const int framesToFlush = 10; // Adjust based on your buffer size
    do
    {

        for(int i = 0; i < framesToFlush; i++)
        {
            while (!StopAll && (Buffer = p->pBuff->GetReadBuffer(2*SPB)) == 0)
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
            if (StopAll) return 0.0;
        }
        p->pBuff->Reset();
        while (!StopAll && (Buffer = p->pBuff->GetReadBuffer(2*SPB)) == 0)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
        if (StopAll) return 0.0;

        // Make sure to advance to the last Buffer
        LastBuffer = Buffer;
        

        int NewPower = CalculatePower(LastBuffer, 2 * SPB);
        double lPower = -100;
        if(NewPower <= 2)
            Delta = -6;
        else
        {
            lPower = 10 * log10(double(NewPower));
            Delta = lPower - TargetPower;
        }
        if(Delta > 6)
            Delta = 6;
        else
            if(Delta < -6)
                Delta = -6;
        cout<<"Power "<<lPower<<" Delta Gain "<<Delta<<endl;
        
        if(gain == MaxGain)
            break;
        bool ToBreak = false;
        if(abs(Delta) < 0.5)
        {
            //Verify
             while ((Buffer = p->pBuff->GetReadBuffer(2*SPB)) == 0)
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
            }
            LastBuffer = Buffer;
            NewPower = CalculatePower(LastBuffer, 2 * SPB);
            lPower = -100;
            if(NewPower <= 2)
                Delta = -6;
            else
            {
                lPower = 10 * log10(double(NewPower));
                Delta = lPower - TargetPower;
            }
            if(abs(Delta)<= 0.5)
            {   
                ToBreak = true;
                double gain1 = p->usrp->get_rx_gain(0);
                std::cout<<"Verified with Delta"<<Delta<<" Gain "<<gain1;
            }
        }
        if(ToBreak)
            break;

        // Make sure to advance to the last Buffer
        
        gain -= Delta;
        gain = round(gain);
        std::cout << boost::format("Setting RX Gain: %f dB...") % gain << std::endl;
        p->usrp->set_rx_gain(gain, 0);
        gain = p->usrp->get_rx_gain(0);
        std::cout << boost::format("Actual RX Gain: %f dB...") % gain << std::endl
                  << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        NumSteps++;
        if (NumSteps == MaxNumSteps)
        {
            break;
        }

    } while (abs(Delta) >= 0.5);
    /*
    for(int i = 0; i < framesToFlush; i++)
    {
        while ((Buffer = p->pBuff->GetReadBuffer()) == 0)
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    }
    //
    */
        p->pBuff->Reset();
//    WriteDebugParam = true;
    return Delta;
}
