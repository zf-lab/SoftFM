
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <fcntl.h>
#include <algorithm>

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


/** Write audio data to ALSA device. */
class AlsaAudioOutput
{
public:

    /**
     * Construct ALSA output stream.
     *
     * dename       :: ALSA PCM device
     * samplerate   :: audio sample rate in Hz
     * stereo       :: true if the output stream contains stereo data
     */
    AlsaAudioOutput(const std::string& devname,
                    unsigned int samplerate,
                    bool stereo);

    ~AlsaAudioOutput();
    bool write(const SampleVector& samples);
    std::string error();

private:
    // TODO
};
#endif
