
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
static IQSample::value_type rms_level_approx(const IQSampleVector& samples)
{
    unsigned int n = samples.size();
    n = (n + 63) / 64;

    IQSample::value_type level = 0;
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


/* ****************  class PilotPhaseLock  **************** */

// Construct phase-locked loop.
PilotPhaseLock::PilotPhaseLock(double freq, double bandwidth, double minsignal)
{
    /*
     * This is a type-2, 4th order phase-locked loop.
     *
     * Open-loop transfer function:
     *   G(z) = K * (z - Qz) / ((z - Qp) * (z - Qp) * (z - 1) * (z - 1))
     *   K  = 3.125 * (bandwidth * 2 * PI)**3
     *   Qz = exp(-0.2 * bandwidth * 2*PI)
     *   Qp = exp(-2.5 * bandwidth * 2*PI)
     *
     * I don't understand what I'm doing; hopefully it just works.
     */

    // Set min/max locking frequencies.
    m_minfreq = (freq - bandwidth) * 2.0 * M_PI;
    m_maxfreq = (freq + bandwidth) * 2.0 * M_PI;

    // Set valid signal threshold.
    m_minsignal  = minsignal;
    m_lock_delay = int(10.0 / bandwidth);
    m_lock_cnt   = 0;

    // Create 2nd order filter for I/Q representation of phase error.
    // Filter has both poles at z = exp(-2.5 * bandwidth * 2*PI).
    double t = exp(-2.5 * bandwidth * 2.0 * M_PI);
    m_phasor_a1 = -2.0 * t;
    m_phasor_a2 = t * t;
    m_phasor_b0 = 1 + m_phasor_a1 + m_phasor_a2;

    // Create loop filter to stabilize the loop.
    // Filter has one Zero at z = exp(-0.2 * bandwidth * 2*PI).
    m_loopfilter_b0 = 0.5 * bandwidth * 2.0 * M_PI;
    m_loopfilter_b1 = - m_loopfilter_b0 * exp(-0.2 * bandwidth * 2.0 * M_PI);

    // After the loop filter, the phase error is integrated to produce
    // the frequency. Then the frequency is integrated to produce the phase.
    // These two integrators form the two remaining poles, both at z = 1.

    // Reset frequency and phase.
    m_freq  = freq * 2.0 * M_PI;
    m_phase = 0;

    m_phasor_i1 = 0;
    m_phasor_i2 = 0;
    m_phasor_q1 = 0;
    m_phasor_q2 = 0;
}


// Process samples.
void PilotPhaseLock::process(const SampleVector& samples_in,
                             SampleVector& samples_out)
{
    unsigned int n = samples_in.size();

    samples_out.resize(n);

    for (unsigned int i = 0; i < n; i++) {

        // Generate locked pilot tone.
        Sample psin = sin(m_phase);
        Sample pcos = cos(m_phase);
        samples_out[i] = pcos;

        // Multiply locked tone with input.
        Sample x = samples_in[i];
        Sample phasor_i = pcos * x;
        Sample phasor_q = psin * x;

        // Run IQ phase error through low-pass filter.
        phasor_i = m_phasor_b0 * phasor_i
                   - m_phasor_a1 * m_phasor_i1
                   - m_phasor_a2 * m_phasor_i2;
        phasor_q = m_phasor_b0 * phasor_q
                   - m_phasor_a1 * m_phasor_q1
                   - m_phasor_a2 * m_phasor_q2;
        m_phasor_i2 = m_phasor_i1;
        m_phasor_i1 = phasor_i;
        m_phasor_q2 = m_phasor_q1;
        m_phasor_q1 = phasor_q;

        // Convert I/Q ratio to estimate of phase error.
        Sample phase_err;
        if (phasor_i > abs(phasor_q)) {
            // We are within +/- 45 degrees from lock.
            // Use simple linear approximation of arctan.
            phase_err = phasor_q / phasor_i;
        } else if (phasor_q > 0) {
            // We are more than 45 degrees ahead of the input.
            phase_err = 1;
        } else {
            // We are lagging more than 45 degrees behind the input.
            phase_err = -1;
        }

        // Detect signal threshold.
        if (phasor_i > m_minsignal) {
            m_lock_cnt++;
        } else {
            m_lock_cnt = 0;
        }

        // Run phase error through loop filter and update frequency estimate.
        m_freq -= m_loopfilter_b0 * phase_err
                  + m_loopfilter_b1 * m_loopfilter_x1;
        m_loopfilter_x1 = phase_err;

        // Limit frequency to allowable range.
        m_freq = max(m_minfreq, min(m_maxfreq, m_freq));

        // Update locked phase.
        m_phase += m_freq;
        if (m_phase < -2.0 * M_PI)
            m_phase += 2.0 * M_PI;
        else if (m_phase > 2.0 * M_PI)
            m_phase -= 2.0 * M_PI;
    }

    m_lock_cnt = min(m_lock_delay, m_lock_cnt);
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

    // Initialize member fields
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

    // Construct FineTuner
    , m_finetuner(m_tuning_table_size, m_tuning_shift)

    // Construct LowPassFilterFirIQ
    , m_iffilter(10, bandwidth_if / sample_rate_if)

    // Construct PhaseDiscriminator
    , m_phasedisc(freq_dev / sample_rate_if)

    // Construct DownsampleFilter for baseband
    , m_resample_baseband(8 * downsample, 0.4 / downsample, downsample, true)

    // Construct DownsampleFilter for mono channel
    , m_resample_mono(int(sample_rate_if / downsample / 1000.0),
                      bandwidth_pcm * downsample / sample_rate_if,
                      sample_rate_if / downsample / sample_rate_pcm, false)

    // Construct HighPassFilterIir
    , m_dcblock_mono(30.0 / sample_rate_pcm)

    // Construct LowwPassFilterRC
    , m_deemph_mono((deemphasis == 0) ? 1.0 : (deemphasis * sample_rate_pcm * 1.0e-6))

{
    // nothing more to do
}


void FmDecoder::process(const IQSampleVector& samples_in,
                        SampleVector& audio)
{
    // Fine tuning.
    m_finetuner.process(samples_in, m_buf_iftuned);

    // Low pass filter to isolate station.
    m_iffilter.process(m_buf_iftuned, m_buf_iffiltered);

    // Measure IF level.
    double if_rms = rms_level_approx(m_buf_iffiltered);
    m_if_level = 0.95 * m_if_level + 0.05 * if_rms;

    // Extract carrier frequency.
    m_phasedisc.process(m_buf_iffiltered, m_buf_baseband);

    // Downsample baseband signal to reduce processing.
    if (m_downsample > 1) {
        SampleVector tmp(move(m_buf_baseband));
        m_resample_baseband.process(tmp, m_buf_baseband);
    }

    // Measure baseband level.
    double baseband_mean, baseband_rms;
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
