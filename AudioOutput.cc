
#define _FILE_OFFSET_BITS 64

#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <fcntl.h>
#include <algorithm>

#include <alsa/asoundlib.h>

#include "SoftFM.h"
#include "AudioOutput.h"

using namespace std;


/* ****************  class AudioOutput  **************** */

// Encode a list of samples as signed 16-bit little-endian integers.
void AudioOutput::samplesToInt16(const SampleVector& samples,
                                 vector<uint8_t>& bytes)
{
    bytes.resize(2 * samples.size());

    SampleVector::const_iterator i = samples.begin();
    SampleVector::const_iterator n = samples.end();
    vector<uint8_t>::iterator k = bytes.begin();

    while (i != n) {
        Sample s = *(i++);
        s = max(Sample(-1.0), min(Sample(1.0), s));
        long v = lrint(s * 32767);
        unsigned long u = v;
        *(k++) = u & 0xff;
        *(k++) = (u >> 8) & 0xff;
    }
}


/* ****************  class RawAudioOutput  **************** */

// Construct raw audio writer.
RawAudioOutput::RawAudioOutput(const string& filename)
{
    if (filename == "-") {

        m_fd = STDOUT_FILENO;

    } else {

        m_fd = open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
        if (m_fd < 0) {
            m_error  = "can not open '" + filename + "' (" +
                       strerror(errno) + ")";
            m_zombie = true;
            return;
        }

    }
}


// Destructor.
RawAudioOutput::~RawAudioOutput()
{
    // Close file descriptor.
    if (m_fd >= 0 && m_fd != STDOUT_FILENO) {
        close(m_fd);
    }
}


// Write audio data.
bool RawAudioOutput::write(const SampleVector& samples)
{
    if (m_fd < 0)
        return -1;

    // Convert samples to bytes.
    samplesToInt16(samples, m_bytebuf);

    // Write data.
    size_t p = 0;
    size_t n = m_bytebuf.size();
    while (p < n) {

        ssize_t k = ::write(m_fd, m_bytebuf.data() + p, n - p);
        if (k <= 0) {
            if (k == 0 || errno != EINTR) {
                m_error = "write failed (";
                m_error += strerror(errno);
                m_error += ")";
                return false;
            }
        } else {
            p += k;
        }
    }

    return true;
}


#if 0

/** Write audio data as .WAV file. */
class WavAudioOutput
{
public:

    /**
     * Construct .WAV writer.
     *
     * filename     :: file name (including path) or "-" to write to stdout
     * samplerate   :: audio sample rate in Hz
     * stereo       :: true if the output stream contains stereo data
     */
    WavAudioOutput(const std::string& filename,
                   unsigned int samplerate,
                   bool stereo);

    ~WavAudioOutput();
    bool write(const SampleVector& samples);
    std::string error();

private:
// TODO
};

#endif

/* ****************  class AlsaAudioOutput  **************** */

// Construct ALSA output stream.
AlsaAudioOutput::AlsaAudioOutput(const std::string& devname,
                                 unsigned int samplerate,
                                 bool stereo)
{
    m_pcm = NULL;
    m_nchannels = stereo ? 2 : 1;

    int r = snd_pcm_open(&m_pcm, devname.c_str(),
                         SND_PCM_STREAM_PLAYBACK, SND_PCM_NONBLOCK);

    if (r < 0) {
        m_error = "can not open PCM device '" + devname + "' (" +
                  strerror(-r) + ")";
        m_zombie = true;
        return;
    }

    snd_pcm_nonblock(m_pcm, 0);

    r = snd_pcm_set_params(m_pcm,
                           SND_PCM_FORMAT_S16_LE,
                           SND_PCM_ACCESS_RW_INTERLEAVED,
                           m_nchannels,
                           samplerate,
                           1,               // allow soft resampling
                           500000);         // latency in us

    if (r < 0) {
        m_error = "can not set PCM parameters (";
        m_error += strerror(-r);
        m_error += ")";
        m_zombie = true;
    }
}


// Destructor.
AlsaAudioOutput::~AlsaAudioOutput()
{
    // Close device.
    if (m_pcm != NULL) {
        snd_pcm_close(m_pcm);
    }
}


// Write audio data.
bool AlsaAudioOutput::write(const SampleVector& samples)
{
    if (m_zombie)
        return false;

    // Convert samples to bytes.
    samplesToInt16(samples, m_bytebuf);

    // Write data.
    unsigned int p = 0;
    unsigned int n = samples.size() / m_nchannels;
    unsigned int framesize = 2 * m_nchannels;
    while (p < n) {

        int k = snd_pcm_writei(m_pcm,
                               m_bytebuf.data() + p * framesize, n - p);
        if (k < 0) {
            m_error = "write failed (";
            m_error += strerror(errno);
            m_error += ")";
            // After an underrun, ALSA keeps returning error codes until we
            // explicitly fix the stream.
            snd_pcm_recover(m_pcm, k, 0);
            return false;
        } else {
            p += k;
        }
    }

    return true;
}

/* end */
