import sys
import math
import random
import time
import splat
import splat.gen
import splat.data
import splat.filters
from splat import lin2dB, dB2lin

FMIN = lin2dB(20.0)
FMAX = lin2dB(20000.0)
FREQ_0 = lin2dB(440.0)
IMAX = 300000
END = 10.0

def curve(s, t, freq=FREQ_0):
    return s - 0.5 * ((2 * s) - (2 * freq)) * math.exp(
        -((s - freq) * (s - freq)) / dB2lin(120.0 * t / END))

def main(argv):
    gen = splat.gen.SineGenerator(splat.data.Fragment(2, 48000, END),
                                  [(splat.filters.linear_fade, (0.005,))])
    L = 0.2
    T = END - L
    for i in range(IMAX):
        t = random.uniform(0, T)
        flog = random.uniform(FMIN, FMAX)
        f = dB2lin(curve(flog, t))
        if (i % 1000) == 0:
            print((i, t, dB2lin(flog), f))
        gen.run(f, t, (t + L))

    rev = splat.data.Fragment(2, 48000)
    rev.mix(gen.frag)
    splat.filters.reverse(rev)
    gen.frag.mix(rev, gen.frag.duration)
    splat.filters.linear_fade(gen.frag, 0.1)
    padded = splat.data.Fragment(2, 48000, gen.frag.duration + 1.0)
    padded.mix(gen.frag, 0.5)
    padded.normalize()
    padded.save('splat_001.wav')

    return True

if __name__ == '__main__':
    ret = main(sys.argv)
    sys.exit(0 if ret is True else 1)
