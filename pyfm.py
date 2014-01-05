"""
Test lab for FM decoding algorithms.
"""

import sys
import types
import numpy
import numpy.fft


def readRawSamples(fname):
    """Read raw sample file from rtl_sdr."""

    d = numpy.fromfile(fname, dtype=numpy.uint8)
    d = d.astype(numpy.float64)
    d = (d - 128) / 128.0

    return d[::2] - 1j * d[1::2]


def lazyRawSamples(fname, blocklen):
    """Return generator over blocks of raw samples."""

    f = file(fname, 'rb')
    while 1:
        d = f.read(2 * blocklen)
        if len(d) < 2 * blocklen:
            break
        d = numpy.fromstring(d, dtype=numpy.uint8)
        d = d.astype(numpy.float64)
        d = (d - 128) / 128.0
        yield d[::2] - 1j * d[1::2]


def freqShiftIQ(d, freqshift):
    """Shift frequency by multiplication with complex phasor."""

    def g(d, freqshift):
        p = 0
        for b in d:
            n = len(b)
            w = numpy.exp((numpy.arange(n) + p) * (2j * numpy.pi * freqshift))
            p += n
            yield b * w

    if isinstance(d, types.GeneratorType):
        return g(d, freqshift)
    else:
        n = len(d)
        w = numpy.exp(numpy.arange(n) * (2j * numpy.pi * freqshift))
        return d * w


def spectrum(d, fs=1, nfft=None, sortfreq=False):
    """Calculate Welch-style power spectral density.

    fs       :: sample rate, default to 1
    nfft     :: FFT length, default to block length
    sortfreq :: True to put negative freqs in front of positive freqs

    Use Hann window with 50% overlap.

    Return (freq, Pxx)."""

    if not isinstance(d, types.GeneratorType):
        d = [ d ]

    prev = None

    if nfft is not None:
        assert nfft > 0
        w = numpy.hanning(nfft)
        q = numpy.zeros(nfft)

    pos = 0
    i = 0
    for b in d:

        if nfft is None:
            nfft = len(b)
            assert nfft > 0
            w = numpy.hanning(nfft)
            q = numpy.zeros(nfft)

        while pos + nfft <= len(b):

            if pos < 0:
                t = numpy.concatenate((prev[pos:], b[:pos+nfft]))
            else:
                t = b[pos:pos+nfft]

            t *= w
            tq = numpy.fft.fft(t)
            tq *= numpy.conj(tq)
            q += numpy.real(tq)

            del t
            del tq

            pos += (nfft+(i%2)) // 2
            i += 1

        pos -= len(b)
        if pos + len(b) > 0:
            prev = b
        else:
            prev = numpy.concatenate((prev[pos+len(b):], b))

    if i > 0:
        q /= (i * numpy.sum(numpy.square(w)) * fs)

    f = numpy.arange(nfft) * (fs / float(nfft))
    f[nfft//2:] -= fs

    if sortfreq:
        f = numpy.concatenate((f[nfft//2:], f[:nfft//2]))
        q = numpy.concatenate((q[nfft//2:], q[:nfft//2]))

    return f, q


def pll(d, centerfreq, bandwidth):

    minfreq = (centerfreq - bandwidth) * 2 * numpy.pi
    maxfreq = (centerfreq + bandwidth) * 2 * numpy.pi

    w = bandwidth * 2 * numpy.pi
    phasor_a = numpy.poly([ numpy.exp(-1.146*w), numpy.exp(-5.331*w) ])
    phasor_b = numpy.array([ sum(phasor_a) ])

    loopfilter_b = numpy.poly([ numpy.exp(-0.1153*w) ])
    loopfilter_b *= 0.62 * w

    n = len(d)
    y = numpy.zeros(n)
    phasei = numpy.zeros(n)
    phaseq = numpy.zeros(n)
    phaseerr = numpy.zeros(n)
    freq = numpy.zeros(n)
    phase = numpy.zeros(n)
    freq[0] = centerfreq * 2 * numpy.pi

    phasor_i1 = phasor_i2 = 0
    phasor_q1 = phasor_q2 = 0
    loopfilter_x1 = 0

    for i in xrange(n):

        psin = numpy.sin(phase[i])
        pcos = numpy.cos(phase[i])
        y[i] = pcos

        pi = pcos * d[i]
        pq = psin * d[i]

        pi = phasor_b[0] * pi - phasor_a[1] * phasor_i1 - phasor_a[2] * phasor_i2
        pq = phasor_b[0] * pq - phasor_a[1] * phasor_q1 - phasor_a[2] * phasor_q2
        phasor_i2 = phasor_i1
        phasor_i1 = pi
        phasor_q2 = phasor_q1
        phasor_q1 = pq

        phasei[i] = pi
        phaseq[i] = pq

        if pi > abs(pq):
            perr = pq / pi
        elif pq > 0:
            perr = 1
        else:
            perr = -1
        phaseerr[i] = perr

        dfreq = loopfilter_b[0] * perr + loopfilter_b[1] * loopfilter_x1
        loopfilter_x1 = perr

        if i + 1 < n:
            freq[i+1] = min(maxfreq, max(minfreq, freq[i] - dfreq))
            p = phase[i] + freq[i+1]
            if p > 2 * numpy.pi:  p -= 2 * numpy.pi
            if p < -2 * numpy.pi: p += 2 * numpy.pi
            phase[i+1] = p

    return y, phasei, phaseq, phaseerr, freq, phase


