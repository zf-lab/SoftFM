#ifndef SOFTFM_FMDECODE_H
#define SOFTFM_FMDECODE_H

#include "SoftFM.h"
#include "Filter.h"


/* Detect frequency by phase discrimination between successive samples. */
class PhaseDiscriminator
{
public:

    /**
     * Construct phase discriminator.
     *
     * max_freq_dev :: Full scale frequency deviation relative to the
     *                 full sample frequency.
     */
    PhaseDiscriminator(double max_freq_dev);

    /**
     * Process samples.
     * Output is a sequence of frequency estimates, scaled such that
     * output value +/- 1.0 represents the maximum frequency deviation.
     */
    void process(const IQSampleVector& samples_in, SampleVector& samples_out);

private:
    const Sample m_freq_scale_factor;
    IQSample     m_last_sample;
};


/** Phase-locked loop for stereo pilot. */
class PilotPhaseLock
{
public:

    /**
     * Construct phase-locked loop.
     *
     * freq       :: center frequency of capture range relative to sample freq
     *               (0.5 is Nyquist)
     * bandwidth  :: approximate bandwidth relative to sample frequency
     * minsignal  :: minimum pilot amplitude
     */
    PilotPhaseLock(double freq, double bandwidth, double minsignal);

    /** Process samples and extract pilot tone at unit amplitude. */
    void process(const SampleVector& samples_in, SampleVector& samples_out);

    /** Return true if the phase-locked loop is locked. */
    bool locked() const
    {
        return m_lock_cnt >= m_lock_delay;
    }

private:
    Sample  m_minfreq, m_maxfreq;
    Sample  m_phasor_b0, m_phasor_a1, m_phasor_a2;
    Sample  m_phasor_i1, m_phasor_i2, m_phasor_q1, m_phasor_q2;
    Sample  m_loopfilter_b0, m_loopfilter_b1;
    Sample  m_loopfilter_x1;
    Sample  m_freq, m_phase;
    Sample  m_minsignal;
    int     m_lock_delay;
    int     m_lock_cnt;
};


/** Complete decoder for FM broadcast signal. */
class FmDecoder
{
public:
    static constexpr double default_deemphasis    =     50;
    static constexpr double default_bandwidth_if  = 100000;
    static constexpr double default_freq_dev      =  75000;
    static constexpr double default_bandwidth_pcm =  15000;

    /**
     * Construct FM decoder.
     *
     * sample_rate_if   :: IQ sample rate in Hz.
     * tuning_offset    :: Frequency offset in Hz of radio station with respect
     *                     to receiver LO frequency (positive value means
     *                     station is at higher frequency than LO).
     * sample_rate_pcm  :: Audio sample rate.
     * stereo           :: True to enable stereo decoding.
     * deemphasis       :: Time constant of de-emphasis filter in microseconds
     *                     (50 us for broadcast FM, 0 to disable de-emphasis).
     * bandwidth_if     :: Half bandwidth of IF signal in Hz
     *                     (~ 100 kHz for broadcast FM)
     * freq_dev         :: Full scale carrier frequency deviation
     *                     (75 kHz for broadcast FM)
     * bandwidth_pcm    :: Half bandwidth of audio signal in Hz
     *                     (15 kHz for broadcast FM)
     * downsample       :: Downsampling factor to apply after FM demodulation.
     *                     Set to 1 to disable.
     */
    FmDecoder(double sample_rate_if,
              double tuning_offset,
              double sample_rate_pcm,
              bool   stereo=true,
              double deemphasis=50,
              double bandwidth_if=default_bandwidth_if,
              double freq_dev=default_freq_dev,
              double bandwidth_pcm=default_bandwidth_pcm,
              unsigned int downsample=1);

    /**
     * Process IQ samples and return audio samples.
     * 
     * If the decoder is set in stereo mode, samples for left and right
     * channels are interleaved in the output vector (even if no stereo
     * signal is detected). If the decoder is set in mono mode, the output
     * vector only contains samples for one channel.
     */
    void process(const IQSampleVector& samples_in,
                 SampleVector& audio);

    /** Return true if a stereo signal is detected. */
    bool stereo_detected() const
    {
        return m_stereo_detected;
    }

    /** Return actual frequency offset in Hz with respect to receiver LO. */
    double get_tuning_offset() const
    {
        double tuned = - m_tuning_shift * m_sample_rate_if /
                       double(m_tuning_table_size);
        return tuned + m_baseband_mean * m_freq_dev;
    }

    /** Return RMS IF level (where full scale IQ signal is 1.0). */
    double get_if_level() const
    {
        return m_if_level;
    }

    /** Return RMS baseband signal level (where nominal level is 0.707). */
    double get_baseband_level() const
    {
        return m_baseband_level;
    }

// TODO : stuff for stereo pilot locking

private:
    const double    m_sample_rate_if;
    const int       m_tuning_table_size;
    const int       m_tuning_shift;
    const double    m_freq_dev;
    const unsigned int m_downsample;
    const bool      m_stereo_enabled;
    bool            m_stereo_detected;
    double          m_if_level;
    double          m_baseband_mean;
    double          m_baseband_level;

    IQSampleVector  m_buf_iftuned;
    IQSampleVector  m_buf_iffiltered;
    SampleVector    m_buf_baseband;
    SampleVector    m_buf_mono;

    FineTuner           m_finetuner;
    LowPassFilterFirIQ  m_iffilter;
    PhaseDiscriminator  m_phasedisc;
    DownsampleFilter    m_resample_baseband;
    DownsampleFilter    m_resample_mono;
    HighPassFilterIir   m_dcblock_mono;
    LowPassFilterRC     m_deemph_mono;
};

#endif
