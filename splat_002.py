import sys
import time
import math
import random
import splat.data
import splat.gen
import splat.interpol
import splat.scales
from splat import lin2dB, dB2lin

class GaussianGenerator(splat.gen.Generator):

    def __init__(self, subgen, *args, **kw):
        super(GaussianGenerator, self).__init__(subgen.frag, *args, **kw)
        self.subgen = subgen
        self.notes = []
        self.start = None
        self.end = None

    def run(self, freq, start, end, k=80.0, levels=(0.0, 0.0)):
        print("run {0} [{1}] {2} -> {3}".format(freq, k, start, end))
        self.notes.append((freq, start, end, k, levels))

        if self.start is None:
            self.start = start
            self.end = end
        else:
            self.start = min(self.start, start)
            self.end = max(self.end, end)

    def add_particle(self, freq, start, end, levels):
        for n_freq, n_start, n_end, k, n_levels in self.notes:
            if (start >= n_start) and (start <= n_end):
                f = dB2lin(self.focus(lin2dB(n_freq), freq, k))
                l = (0.0, 0.0)
                self.subgen.run(f, start, end, l)

    @classmethod
    def focus(cls, freq, s, k):
        return s - 0.5 * ((2 * s) - (2 * freq)) * math.exp(
            -((s - freq) * (s - freq)) / dB2lin(k))


class SplineGaussianGenerator(GaussianGenerator):

    def __init__(self, *args, **kw):
        super(SplineGaussianGenerator, self).__init__(*args, **kw)
        self._freq_pts = []
        self._freq_spline = []
        self._k_pts = []
        self._k_spline = []

    def add_freq(self, time, freq, slope=None):
        if slope is None:
            pt = (time, freq)
        else:
            pt = (time, freq, slope)
        self._freq_pts.append(pt)

    def add_k(self, time, k, slope=None):
        if slope is None:
            pt = (time, k)
        else:
            pt = (time, k, slope)
        self._k_pts.append(pt)

    def set_f(self, freqs):
        self._freq_pts = freqs

    def set_k(self, ks):
        self._k_pts = ks

    def make_splines(self):
        self._freq_spline = splat.interpol.Spline(self._freq_pts)
        self._k_spline = splat.interpol.Spline(self._k_pts)
        self.start = self._freq_pts[0][0]
        self.end = self._freq_pts[-1][0]

        print("Freq spline:")
        for x0, x1, pol in self._freq_spline.polynomials:
            print('[{0} {1}] {2}'.format(
                    x0, x1, ' '.join([("%04.3f" % x) for x in pol.pol])))

        print("K spline:")
        for x0, x1, pol in self._k_spline.polynomials:
            print('[{0} {1}] {2}'.format(
                    x0, x1, ' '.join([("%04.3f" % x) for x in pol.pol])))

    def add_particle(self, freq, start, end, levels):
        f0 = self._freq_spline.value(start)
        k = self._k_spline.value(start)
        f = dB2lin(self.focus(lin2dB(f0), freq, k))
        l = (0.0, 0.0)
        self.subgen.run(f, start, end, l)


def main(argv):
    T = 0.4
    m = 10000

    subgen = splat.gen.OvertonesGenerator(
        splat.data.Fragment(),
        [(splat.filters.linear_fade, (0.0017,)),
         (splat.filters.dec_envelope, (100.0, 1.2))])
    subgen.overtones = { 1.0: 0.0, 2.0: -9.0, 3.0: -18.0 }
    delays = splat.filters.reverb_delays(length=0.5, dec=0.3)
    subgen.filters.append(splat.filters.reverb, (delays, 0.2, 4.0))

    # Main generator
    gen = SplineGaussianGenerator(subgen)

    # Notes
    scale = splat.scales.LogScale(fund=440.0)
    gen.set_f([(0.0, scale['C-1']), (1.5, scale['C-1']),
               (2.0, scale['E-1'], 0.0), (2.9, scale['E-1']),
               (3.0, scale['G-1'], 0.0), (3.9, scale['G-1']),
               (4.0, scale['C'], 0.0), (4.9, scale['C']),
               (5.0, scale['E'], 0.0), (5.9, scale['E']),
               (6.0, scale['G'], 0.0), (6.9, scale['G']),
               (7.0, scale['C'], 0.0), (7.5, scale['C-1'], 0.0),
               (8.0, scale['C-1'])])
    gen.set_k([(0.0, 50.0), (3.0, 60.0), (6.0, 75.0), (8.0, 60.0)])
    gen.make_splines()

    # Equalisation
    # ToDo: use spline
    eq = [(20.0, -36.0), (150.0, -12.0), (500.0, -3.0), (200.0, 0.0),
          (1000.0, -12.0), (3000.0, -38.0), (8000.0, -60.0)]
    f_min_log = splat.lin2dB(eq[0][0])
    f_max_log = splat.lin2dB(eq[-1][0])

    # Envelope
    time_pts = [(0.0, 0.0), (0.01, 0.0, 0.0), (0.5, 0.0, 0.0), (2.0, 0.3),
                (3.0, 1.0, 0.0), (4.5, 0.45), (6.0, 0.6),
                (7.8, 0.0, 0.0), (8.0, 0.0)]
    time_spline = splat.interpol.Spline(time_pts)

    # Generating now...

    n = 100
    density = m / n

    for y0 in range(1, n + 1):
        if (y0 % 10) == 0:
            print(y0)
        y0 = float(y0) / n
        seg = time_spline.slices(y0, gen.start, gen.end)

        for x0, x1 in seg:
            n_events = int((x1 - x0) * density)

            for i in range(n_events):
                start = random.uniform(x0, x1)
                end = random.uniform((start + 0.05), (start + 0.10))
                freq = random.uniform(f_min_log, f_max_log)

                for i, (f, g) in enumerate(eq):
                    if freq < f:
                        break

                eq0 = eq[i - 1]
                eq1 = eq[i]
                g = ((freq - eq0[0]) * ((eq1[1] - eq0[1]) / (eq1[0] - eq0[0]))
                     + eq0[1])
                g += random.uniform(-18.0, 18.0)
                levels = tuple((g + random.uniform(-3.0, 3.0)
                                for c in range(gen.frag.channels)))

                gen.add_particle(freq, start, end, levels)

    print("saving file...")
    gen.frag.normalize(-0.1, False)
    splat.filters.linear_fade(gen.frag, 0.1)
    gen.frag.save('splat_002.wav')

    return True

if __name__ == '__main__':
    ret = main(sys.argv)
    sys.exit(0 if ret is True else 1)
