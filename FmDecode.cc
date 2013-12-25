
#include <cmath>

#include "FmDecode.h"

using namespace std;


/** Fast approximation of atan function. */
static inline Sample fast_atan(Sample x)
{
    // http://stackoverflow.com/questions/7378187/approximating-inverse-trigonometric-funcions

    Sample y = 1;
    Sample p = 0;

    if (x < 0) {
        x = -x;
        y = -1;
    }

    if (x > 1) {
        p = y;
        y = -y;
        x = 1 / x;
    }

    const Sample b = 0.596227;
    y *= (b*x + x*x) / (1 + 2*b*x + x*x);

    return (y + p) * Sample(M_PI_2);
}


/** Compute RMS level over a small prefix of the specified sample vector. */
static Sample rms_level_approx(const IQSampleVector& samples)
{
    unsigned int n = samples.size();
    n = (n + 63) / 64;

    Sample level = 0;
    for (unsigned int i = 0; i < n; i++) {
        const IQSample& s = samples[i];
        IQSample::value_type re = s.real(), im = s.imag();
        level += re * re + im * im;
    }

    return sqrt(level / n);
}


/* ****************  class PhaseDiscriminator  **************** */

// Construct phase discriminator.
PhaseDiscriminator::PhaseDiscriminator(double max_freq_dev)
    : m_freq_scale_factor(1.0 / (max_freq_dev * 2.0 * M_PI))
{ }


// Process samples.
void PhaseDiscriminator::process(const IQSampleVector& samples_in,
                                 SampleVector& samples_out)
{
    unsigned int n = samples_in.size();
    IQSample s0 = m_last_sample;

    samples_out.resize(n);

    for (unsigned int i = 0; i < n; i++) {
        IQSample s1(samples_in[i]);
        IQSample d(conj(s0) * s1);
// TODO : implement fast approximation of atan2
        Sample w = atan2(d.imag(), d.real());
        samples_out[i] = w * m_freq_scale_factor;
        s0 = s1;
    }

    m_last_sample = s0;
}


/* ****************  class FmDecoder  **************** */

FmDecoder::FmDecoder(double sample_rate_if,
                     double tuning_offset,
                     double sample_rate_pcm,
                     bool   stereo,
                     double deemphasis,
                     double bandwidth_if,
                     double freq_dev,
                     double bandwidth_pcm,
                     unsigned int downsample)
    : m_sample_rate_if(sample_rate_if)
    , m_tuning_table_size(64)
    , m_tuning_shift(lrint(-64.0 * tuning_offset / sample_rate_if))
    , m_freq_dev(freq_dev)
    , m_downsample(downsample)
    , m_stereo_enabled(stereo)
    , m_stereo_detected(false)
    , m_if_level(0)
    , m_baseband_mean(0)
    , m_baseband_level(0)
    , m_finetuner(m_tuning_table_size, m_tuning_shift)
    , m_iffilter(10, bandwidth_if / sample_rate_if)
    , m_phasedisc(freq_dev / sample_rate_if)
    , m_resample_baseband(6 * downsample,
                          0.5 / downsample,
                          downsample, true)
    , m_resample_mono(int(15 * sample_rate_if / downsample / bandwidth_pcm),
                      bandwidth_pcm * downsample / sample_rate_if,
                      sample_rate_if / downsample / sample_rate_pcm, false)
    , m_dcblock_mono(30.0 / sample_rate_pcm)
    , m_deemph_mono((deemphasis == 0) ? 1.0 : (deemphasis * sample_rate_pcm * 1.0e-6))
{
}


void FmDecoder::process(const IQSampleVector& samples_in,
                        SampleVector& audio)
{
    // Fine tuning.
    m_finetuner.process(samples_in, m_buf_iftuned);

    // Low pass filter to isolate station.
    m_iffilter.process(m_buf_iftuned, m_buf_iffiltered);

    // Measure IF level.
    Sample if_rms = rms_level_approx(m_buf_iffiltered);
    m_if_level = 0.95 * m_if_level + 0.05 * if_rms;

    // Extract carrier frequency.
    m_phasedisc.process(m_buf_iffiltered, m_buf_baseband);

    // Downsample baseband signal to reduce processing.
    if (m_downsample > 1) {
        SampleVector tmp(move(m_buf_baseband));
        m_resample_baseband.process(tmp, m_buf_baseband);
    }

    // Measure baseband level.
    Sample baseband_mean, baseband_rms;
    samples_mean_rms(m_buf_baseband, baseband_mean, baseband_rms);
    m_baseband_mean  = 0.95 * m_baseband_mean + 0.05 * baseband_mean;
    m_baseband_level = 0.95 * m_baseband_level + 0.05 * baseband_rms;

// TODO : stereo decoding

    // Extract mono audio signal.
    m_resample_mono.process(m_buf_baseband, m_buf_mono);

    // DC blocking and de-emphasis.
    m_dcblock_mono.processInPlace(m_buf_mono);
    m_deemph_mono.processInPlace(m_buf_mono);

// TODO : stereo mixing
    audio = move(m_buf_mono);
}

/* end */
