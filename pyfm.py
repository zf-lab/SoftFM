"""
Test lab for FM decoding algorithms.
"""

import sys
import types
import numpy
import numpy.fft
import scipy.signal


def readRawSamples(fname):
    """Read raw sample file from rtl_sdr."""

    d = numpy.fromfile(fname, dtype=numpy.uint8)
    d = d.astype(numpy.float64)
    d = (d - 128) / 128.0

    return d[::2] + 1j * d[1::2]


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
        yield d[::2] + 1j * d[1::2]


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


def firFilter(d, coeff):
    """Apply FIR filter to sample stream."""

    # lazy version
    def g(d, coeff):
        prev = None
        for b in d:
            if prev is None:
                yield scipy.signal.lfilter(coeff, 1, b)
                prev = b
            else:
                k = min(len(prev), len(coeff))
                x = numpy.concatenate((prev[-k:], b))
                y = scipy.signal.lfilter(coeff, 1, x)
                yield y[k:]
                if len(coeff) > len(b):
                    prev = x
                else:
                    prev = b

    if isinstance(d, types.GeneratorType):
        return g(d, coeff)
    else:
        return scipy.signal.lfilter(coeff, 1, d)


def quadratureDetector(d):
    """FM frequency detector based on quadrature demodulation."""

    # lazy version
    def g(d):
        prev = None
        for b in d:
            if prev is None:
                yield numpy.angle(b[1:] * b[:-1].conj())
            else:
                x = numpy.concatenate((prev[-1:], b[:-1]))
                yield numpy.angle(b * x.conj())
            prev = b

    if isinstance(d, types.GeneratorType):
        return g(d)
    else:
        return numpy.angle(d[1:] * d[:-1].conj())


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
    """Simulate the stereo pilot PLL."""

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


def pilotLevel(d, fs, freqshift, bw=150.0e3):
    """Calculate level of the 19 kHz pilot vs noise floor in the guard band.

    d         :: block of raw I/Q samples or lazy I/Q sample stream
    fs        :: sample frequency in Hz
    freqshift :: frequency offset in Hz
    bw        :: half-bandwidth of IF signal in Hz

    Return (pilot_power, guard_floor, noise)
    where pilot_power is the power of the pilot tone in dB
          guard_floor is the noise floor in the guard band in dB/Hz
          noise       is guard_floor - pilot_power in dB/Hz
    """

    # Shift frequency
    if freqshift != 0:
        d = freqShiftIQ(d, freqshift / float(fs))

    # Filter
    b = scipy.signal.firwin(31, 2.0 * bw / fs, window='nuttall')
    d = firFilter(d, b)

    # Demodulate FM.
    d = quadratureDetector(d)

    # Power spectral density.
    f, q = spectrum(d, fs=fs, sortfreq=False)

    # Locate 19 kHz bin.
    k19 = int(19.0e3 * len(q) / fs)
    kw  = 5 + int(100.0 * len(q) / fs)
    k19 = k19 - kw + numpy.argmax(q[k19-kw:k19+kw])

    # Calculate pilot power.
    p19 = numpy.sum(q[k19-1:k19+2]) * fs * 0.75 / len(q)

    # Calculate noise floor in guard band.
    k17 = int(17.0e3 * len(q) / fs)
    k18 = int(18.0e3 * len(q) / fs)
    guard = numpy.mean(q[k17:k18])

    p19db   = 10 * numpy.log10(p19)
    guarddb = 10 * numpy.log10(guard)

    return (p19db, guarddb, guarddb - p19db)

