"""
Test PLL algorithm.

Usage: testpll.py baseband.dat centerfreq bandwidth > output.dat

  baseband.dat   Raw 16-bit signed little-endian sample stream
  centerfreq     Center frequency relative to sample frequency (0.5 = Nyquist)
  bandwidth      Approximate bandwidth of PLL relative to sample frequency
  output.dat     ASCII file with space-separated data
"""

import sys
import numpy


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


def main():

    if len(sys.argv) != 4:
        print >>sys.stderr, __doc__
        sys.exit(1)

    infile      = sys.argv[1]
    centerfreq  = float(sys.argv[2])
    bandwidth   = float(sys.argv[3])

    d = numpy.fromfile(infile, '<i2')
    d = d.astype(numpy.float64) / 32767.0

    (y, phasei, phaseq, phaseerr, freq, phase) = pll(d, centerfreq, bandwidth)

    print '#output phasei, phaseq, phaseerr freq phase'
    for i in xrange(len(y)):
        print y[i], phasei[i], phaseq[i], phaseerr[i], freq[i], phase[i]


if __name__ == '__main__':
    main()

