#ifndef SOFTFM_H
#define SOFTFM_H

#include <complex>
#include <vector>

typedef std::complex<float> IQSample;
typedef std::vector<IQSample> IQSampleVector;

typedef float Sample;
typedef std::vector<Sample> SampleVector;

typedef float Coeff;
typedef std::vector<Coeff> CoeffVector;


/** Compute mean and RMS over a sample vector. */
inline void samples_mean_rms(const SampleVector& samples, Sample& mean, Sample& rms)
{
    Sample vsum = 0;
    Sample vsumsq = 0;

    unsigned int n = samples.size();
    for (unsigned int i = 0; i < n; i++) {
        Sample v = samples[i];
        vsum   += v;
        vsumsq += v * v;
    }

    mean = vsum / n;
    rms  = sqrt(vsumsq / n);
}

#endif
