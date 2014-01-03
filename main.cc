/*
 *  SoftFM - Software decoder for FM broadcast radio with RTL-SDR
 *
 *  Copyright (C) 2013, Joris van Rantwijk.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, see http://www.gnu.org/licenses/gpl-2.0.html
 */

#include <cstdlib>
#include <cstdio>
#include <climits>
#include <cmath>
#include <csignal>
#include <cstring>
#include <atomic>
#include <condition_variable>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <unistd.h>
#include <getopt.h>

#include "SoftFM.h"
#include "RtlSdrSource.h"
#include "FmDecode.h"
#include "AudioOutput.h"

using namespace std;


/** Flag is set on SIGINT / SIGTERM. */
static atomic_bool stop_flag(false);


/** Buffer to move sample data between threads. */
template <class Element>
class DataBuffer
{
public:
    /** Constructor. */
    DataBuffer()
        : m_qlen(0)
        , m_end_marked(false)
    { }

    /** Add samples to the queue. */
    void push(vector<Element>&& samples)
    {
        if (!samples.empty()) {
            unique_lock<mutex> lock(m_mutex);
            m_qlen += samples.size();
            m_queue.push(move(samples));
            lock.unlock();
            m_cond.notify_all();
        }
    }

    /** Mark the end of the data stream. */
    void push_end()
    {
        unique_lock<mutex> lock(m_mutex);
        m_end_marked = true;
        lock.unlock();
        m_cond.notify_all();
    }

    /** Return number of samples in queue. */
    size_t queued_samples()
    {
        unique_lock<mutex> lock(m_mutex);
        return m_qlen;
    }

    /**
     * If the queue is non-empty, remove a block from the queue and
     * return the samples. If the end marker has been reached, return
     * an empty vector. If the queue is empty, wait until more data is pushed
     * or until the end marker is pushed.
     */
    vector<Element> pull()
    {
        vector<Element> ret;
        unique_lock<mutex> lock(m_mutex);
        while (m_queue.empty() && !m_end_marked)
            m_cond.wait(lock);
        if (!m_queue.empty()) {
            m_qlen -= m_queue.front().size();
            swap(ret, m_queue.front());
            m_queue.pop();
        }
        return ret;
    }

    /** Return true if the end has been reached at the Pull side. */
    bool pull_end_reached()
    {
        unique_lock<mutex> lock(m_mutex);
        return m_qlen == 0 && m_end_marked;
    }

    /** Wait until the buffer contains minfill samples or an end marker. */
    void wait_buffer_fill(size_t minfill)
    {
        unique_lock<mutex> lock(m_mutex);
        while (m_qlen < minfill && !m_end_marked)
            m_cond.wait(lock);
    }

private:
    size_t              m_qlen;
    bool                m_end_marked;
    queue<vector<Element>> m_queue;
    mutex               m_mutex;
    condition_variable  m_cond;
};


/** Simple linear gain adjustment. */
void adjust_gain(SampleVector& samples, double gain)
{
    for (unsigned int i = 0, n = samples.size(); i < n; i++) {
        samples[i] *= gain;
    }
}


/**
 * Read data from source device and put it in a buffer.
 *
 * This code runs in a separate thread.
 * The RTL-SDR library is not capable of buffering large amounts of data.
 * Running this in a background thread ensures that the time between calls
 * to RtlSdrSource::get_samples() is very short. 
 */
void read_source_data(RtlSdrSource *rtlsdr, DataBuffer<IQSample> *buf)
{
    IQSampleVector iqsamples;

    while (!stop_flag.load()) {

        if (!rtlsdr->get_samples(iqsamples)) {
            fprintf(stderr, "ERROR: RtlSdr: %s\n", rtlsdr->error().c_str());
            exit(1);
        }

        buf->push(move(iqsamples));
    }

    buf->push_end();
}


/**
 * Get data from output buffer and write to output stream.
 *
 * This code runs in a separate thread.
 */
void write_output_data(AudioOutput *output, DataBuffer<Sample> *buf,
                       unsigned int buf_minfill)
{
    while (!stop_flag.load()) {

        if (buf->queued_samples() == 0) {
            // The buffer is empty. Perhaps the output stream is consuming
            // samples faster than we can produce them. Wait until the buffer
            // is back at its nominal level to make sure this does not happen
            // too often.
            buf->wait_buffer_fill(buf_minfill);
        }

        if (buf->pull_end_reached()) {
            // Reached end of stream.
            break;
        }

        // Get samples from buffer and write to output.
        SampleVector samples = buf->pull();
        output->write(samples);
        if (!(*output)) {
            fprintf(stderr, "ERROR: AudioOutput: %s\n", output->error().c_str());
        }
    }
}


/** Handle Ctrl-C and SIGTERM. */
static void handle_sigterm(int sig)
{
    stop_flag.store(true);

    string msg = "Got signal ";
    msg += strsignal(sig);
    msg += ", stopping ...\n";

    const char *s = msg.c_str();
    write(STDERR_FILENO, s, strlen(s));
}


void usage()
{
    fprintf(stderr,
            "Usage: softfm -f freq [options]\n"
            "  -f freq       Frequency of radio station in Hz\n"
            "  -d devidx     RTL-SDR device index (default 0)\n"
            "  -s ifrate     IF sample rate in Hz (default 1000000)\n"
            "  -r pcmrate    Audio sample rate in Hz (default 48000 Hz)\n"
            "  -M            Disable stereo decoding\n"
            "  -R filename   Write audio data as raw S16_LE samples\n"
            "                use filename '-' to write to stdout\n"
            "  -W filename   Write audio data to .WAV file\n"
            "  -P [device]   Play audio via ALSA device (default 'default')\n"
            "  -b seconds    Set audio buffer size in seconds\n"
            "\n");
}


void badarg(const char *label)
{
    usage();
    fprintf(stderr, "ERROR: Invalid argument for %s\n", label);
    exit(1);
}


bool parse_opt(const char *s, int& v)
{
    char *endp;
    long t = strtol(s, &endp, 10);
    if (endp == s || *endp != '\0' || t < INT_MIN || t > INT_MAX)
        return false;
    v = t;
    return true;
}


bool parse_opt(const char *s, double& v)
{
    char *endp;
    v = strtod(s, &endp);
    return (endp != s && *endp == '\0');
}


int main(int argc, char **argv)
{
    double  freq    = -1;
    int     devidx  = 0;
    double  ifrate  = 1.0e6;
    int     pcmrate = 48000;
    bool    stereo  = true;
    enum OutputMode { MODE_RAW, MODE_WAV, MODE_ALSA };
    OutputMode outmode = MODE_ALSA;
    string  filename;
    string  devname("default");
    double  bufsecs = -1;

    fprintf(stderr,
            "SoftFM - Software decoder for FM broadcast radio with RTL-SDR\n");

    const struct option longopts[] = {
        { "freq",       1, NULL, 'f' },
        { "dev",        1, NULL, 'd' },
        { "ifrate",     1, NULL, 's' },
        { "pcmrate",    1, NULL, 'r' },
        { "mono",       0, NULL, 'M' },
        { "raw",        1, NULL, 'R' },
        { "wav",        1, NULL, 'W' },
        { "play",       2, NULL, 'P' },
        { "buffer",     1, NULL, 'b' },
        { NULL,         0, NULL, 0 } };

    int c, longindex;
    while ((c = getopt_long(argc, argv,
                            "f:d:s:r:MR:W:P::b:",
                            longopts, &longindex)) >= 0) {
        switch (c) {
            case 'f':
                if (!parse_opt(optarg, freq) || freq <= 0) {
                    badarg("-f");
                }
                break;
            case 'd':
                if (!parse_opt(optarg, devidx) || devidx < 0) {
                    badarg("-d");
                }
                break;
            case 's':
                if (!parse_opt(optarg, ifrate) || ifrate <= 0) {
                    badarg("-s");
                }
                break;
            case 'r':
                if (!parse_opt(optarg, pcmrate) || pcmrate < 1) {
                    badarg("-r");
                }
                break;
            case 'M':
                stereo = false;
                break;
            case 'R':
                outmode = MODE_RAW;
                filename = optarg;
                break;
            case 'W':
                outmode = MODE_WAV;
                filename = optarg;
                break;
            case 'P':
                outmode = MODE_ALSA;
                if (optarg != NULL)
                    devname = optarg;
                break;
            case 'b':
                if (!parse_opt(optarg, bufsecs) || bufsecs < 0) {
                    badarg("-b");
                }
                break;
            default:
                usage();
                fprintf(stderr, "ERROR: Unknown option\n");
                exit(1);
        }
    }

    if (freq <= 0) {
        usage();
        fprintf(stderr, "ERROR: Specify a tuning frequency\n");
        exit(1);
    }

    if (3 * FmDecoder::default_bandwidth_if > ifrate) {
        fprintf(stderr, "ERROR: IF sample rate must be at least %.0f Hz\n",
                        3 * FmDecoder::default_bandwidth_if);
        exit(1);
    }

    // Catch Ctrl-C and SIGTERM
    struct sigaction sigact;
    sigact.sa_handler = handle_sigterm;
    sigemptyset(&sigact.sa_mask);
    sigact.sa_flags = SA_RESETHAND;
    if (sigaction(SIGINT, &sigact, NULL) < 0) {
        fprintf(stderr, "WARNING: can not install SIGINT handler (%s)\n",
                strerror(errno));
    }
    if (sigaction(SIGTERM, &sigact, NULL) < 0) {
        fprintf(stderr, "WARNING: can not install SIGTERM handler (%s)\n",
                strerror(errno));
    }

    // Intentionally tune at a higher frequency to avoid DC offset.
    double tuner_freq = freq;
    if (ifrate >= 5 * FmDecoder::default_bandwidth_if) {
        tuner_freq += 0.25 * ifrate;
    }

    // Open RTL-SDR device.
    RtlSdrSource rtlsdr(devidx);
    if (!rtlsdr) {
        fprintf(stderr, "ERROR: RtlSdr: %s\n", rtlsdr.error().c_str());
        exit(1);
    }

    // Configure RTL-SDR device and start streaming.
    rtlsdr.configure(ifrate, tuner_freq, -1);
    if (!rtlsdr) {
        fprintf(stderr, "ERROR: RtlSdr: %s\n", rtlsdr.error().c_str());
        exit(1);
    }

    tuner_freq = rtlsdr.get_frequency();
    fprintf(stderr, "device tuned for %.6f MHz\n", tuner_freq * 1.0e-6);

    ifrate = rtlsdr.get_sample_rate();
    fprintf(stderr, "IF sample rate %.0f Hz\n", ifrate);

    // Create source data queue.
    DataBuffer<IQSample> source_buffer;

    // Start reading from device in separate thread.
    thread source_thread(read_source_data, &rtlsdr, &source_buffer);

    // The baseband signal is empty above 100 kHz, so we can
    // downsample to ~ 200 kS/s without loss of information.
    // This will speed up later processing stages.
    unsigned int downsample = max(1, int(ifrate / 215.0e3));
    fprintf(stderr, "baseband downsampling factor %u\n", downsample);

    // Prevent aliasing at very low output sample rates.
    double bandwidth_pcm = min(FmDecoder::default_bandwidth_pcm,
                               0.45 * pcmrate);
    fprintf(stderr, "audio sample rate %u Hz\n", pcmrate);
    fprintf(stderr, "audio bandwidth %.3f kHz\n", bandwidth_pcm * 1.0e-3);

    // Prepare decoder.
    FmDecoder fm(ifrate,                            // sample_rate_if
                 freq - tuner_freq,                 // tuning_offset
                 pcmrate,                           // sample_rate_pcm
                 stereo,                            // stereo
                 FmDecoder::default_deemphasis,     // deemphasis,
                 FmDecoder::default_bandwidth_if,   // bandwidth_if
                 FmDecoder::default_freq_dev,       // freq_dev
                 bandwidth_pcm,                     // bandwidth_pcm
                 downsample);                       // downsample

    // Calculate number of samples in audio buffer.
    unsigned int outputbuf_samples = 0;
    if (bufsecs < 0 &&
        (outmode == MODE_ALSA || (outmode == MODE_RAW && filename == "-"))) {
        // Set default buffer to 1 second for interactive output streams.
        outputbuf_samples = pcmrate;
    } else if (bufsecs > 0) {
        // Calculate nr of samples for configured buffer length.
        outputbuf_samples = (unsigned int)(bufsecs * pcmrate);
    }
    if (outputbuf_samples > 0) {
        fprintf(stderr, "output buffer %.1f seconds\n",
                outputbuf_samples / double(pcmrate));
    }

    // Prepare output writer.
    unique_ptr<AudioOutput> audio_output;
    switch (outmode) {
        case MODE_RAW:
            audio_output.reset(new RawAudioOutput(filename));
            break;
        case MODE_WAV:
            audio_output.reset(new WavAudioOutput(filename, pcmrate, stereo));
            break;
        case MODE_ALSA:
            audio_output.reset(new AlsaAudioOutput(devname, pcmrate, stereo));
            break;
    }

    if (!(*audio_output)) {
        fprintf(stderr, "ERROR: AudioOutput: %s\n",
                        audio_output->error().c_str());
        exit(1);
    }

    // If buffering enabled, start background output thread.
    DataBuffer<Sample> output_buffer;
    thread output_thread;
    if (outputbuf_samples > 0) {
        unsigned int nchannel = stereo ? 2 : 1;
        output_thread = thread(write_output_data,
                               audio_output.get(),
                               &output_buffer,
                               outputbuf_samples * nchannel);
    }

    SampleVector audiosamples;
    bool inbuf_length_warning = false;
    double audio_level = 0;
    bool got_stereo = false;

    // Main loop.
    for (unsigned int block = 0; !stop_flag.load(); block++) {

        // Check for overflow of source buffer.
        if (!inbuf_length_warning &&
            source_buffer.queued_samples() > 10 * ifrate) {
            fprintf(stderr,
                    "\nWARNING: Input buffer is growing (system too slow)\n");
            inbuf_length_warning = true;
        }

        // Pull next block from source buffer.
        IQSampleVector iqsamples = source_buffer.pull();
        if (iqsamples.empty())
            break;

        // Decode FM signal.
        fm.process(iqsamples, audiosamples);

        // Measure audio level.
        double audio_mean, audio_rms;
        samples_mean_rms(audiosamples, audio_mean, audio_rms);
        audio_level = 0.95 * audio_level + 0.05 * audio_rms;

        adjust_gain(audiosamples, 0.5);

// TODO : investigate I/Q imbalance to fix Radio4 noise

        fprintf(stderr,
                "\rblk=%6d  freq=%8.4fMHz  IF=%+5.1fdB  BB=%+5.1fdB  audio=%+5.1fdB ",
                block,
                (tuner_freq + fm.get_tuning_offset()) * 1.0e-6,
                20*log10(fm.get_if_level()),
                20*log10(fm.get_baseband_level()) + 3.01,
                20*log10(audio_level) + 3.01);
        if (outputbuf_samples > 0) {
            unsigned int nchannel = stereo ? 2 : 1;
            size_t buflen = output_buffer.queued_samples();
            fprintf(stderr,
                    " buf=%.1fs ",
                    buflen / nchannel / double(pcmrate));
        }
        fflush(stderr);

        if (fm.stereo_detected() != got_stereo) {
            got_stereo = fm.stereo_detected();
            if (got_stereo)
                fprintf(stderr, "\ngot stereo signal (pilot level = %f)\n",
                        fm.get_pilot_level());
            else
                fprintf(stderr, "\nlost stereo signal\n");
        }

        // Throw away first block. It is noisy because IF filters
        // are still starting up.
        if (block > 0) {

            // Write samples to output.
            if (outputbuf_samples > 0) {
                // Buffered write.
                output_buffer.push(move(audiosamples));
            } else {
                // Direct write.
                audio_output->write(audiosamples);
            }
        }

    }

    fprintf(stderr, "\n");

    // Join background threads.
    source_thread.join();
    if (outputbuf_samples > 0) {
        output_buffer.push_end();
        output_thread.join();
    }

    // No cleanup needed; everything handled by destructors.

    return 0;
}

/* end */
