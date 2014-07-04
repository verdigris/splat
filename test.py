# Splat - test.py
#
# Copyright (C) 2012, 2013, 2014 Guillaume Tucker <guillaume@mangoz.org>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import md5
import math
import unittest
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
import splat
import splat.data
import splat.gen
import splat.sources
import splat.interpol

splat.check_version((1, 3))

# -----------------------------------------------------------------------------
# Splat test base class

class SplatTest(unittest.TestCase):

    def setUp(self):
        self._places = 12 if splat.NATIVE_SAMPLE_WIDTH == 64 else 4

    def assert_md5(self, frags, hexdigest):
        if isinstance(frags, splat.data.Fragment):
            frags = [frags]
        for frag in frags:
            md5sum = frag.md5(sample_type=splat.SAMPLE_INT, sample_width=16)
            self.assertEqual(md5sum, hexdigest,
                             "MD5 mismatch: {0} {1}".format(md5sum, hexdigest))

    def assert_samples(self, frag, samples, places=None):
        if places is None:
            places = self._places
        for n, s in samples.iteritems():
            for c, d in zip(frag[n], s):
                self.assertAlmostEqual(
                    c, d, places,
                    "Samples mismatch {0}: {1} {2}".format(n, frag[n], s))

# -----------------------------------------------------------------------------
# test cases

class FragmentTest(SplatTest):

    def test_frag(self):
        """Fragment"""
        frag = splat.data.Fragment(duration=1.0)
        self.assert_samples(frag, {int(len(frag) / 2): (0.0, 0.0)})
        self.assert_md5(frag, 'fe384f668da282694c29a84ebd33481d')

    def test_frag_md5(self):
        """Fragment.md5"""
        frag = splat.data.Fragment()
        splat.gen.SineGenerator(frag=frag).run(0.0, 0.345, 123.0)
        for st, sw in splat.audio_formats:
            md5sum = md5.new(frag.export_bytes(st, sw)).hexdigest()
            self.assertEqual(md5sum, frag.md5(st, sw))

    def test_frag_offset(self):
        """Fragment.offset"""
        offset = 1234.5
        frag = splat.data.Fragment(duration=1.0)
        n = int(0.5678 * frag.duration * frag.rate)
        frag.offset(offset)
        frag_twice = splat.data.Fragment(duration=1.0)
        gen = splat.gen.SineGenerator(frag=frag_twice)
        gen.run(0.0, 1.0, 1000.0)
        frag_twice.offset(offset)
        frag_twice_ref = frag_twice[n][0]
        offset2 = -57.25
        offset_check = frag_twice_ref + offset2
        frag_twice.offset(offset2)
        frag_sig = splat.data.Fragment(duration=1.0)
        frag_sig.offset(lambda x: offset)
        self.assert_samples(frag, {n: (offset, offset)})
        self.assert_samples(frag_sig, {n: (offset, offset)})
        self.assert_samples(frag_twice, {n: (offset_check, offset_check)})
        self.assert_md5([frag, frag_sig], '248070c79f99014cf800d05ea81e0679')

    def test_frag_import_bytes(self):
        """Fragment.import_bytes"""
        frag = splat.data.Fragment()
        splat.gen.SineGenerator(frag).run(0.1, 2.7, 3456.7)
        for st, sw in splat.audio_formats:
            ref_bytes = frag.export_bytes(st, sw)
            ref_md5 = md5.new(ref_bytes).hexdigest()
            imp = splat.data.Fragment()
            imp.import_bytes(ref_bytes, frag.rate, frag.channels, st, sw)
            exp_bytes = imp.export_bytes(st, sw)
            exp_md5 = md5.new(exp_bytes).hexdigest()
            self.assertEqual(ref_md5, exp_md5,
                             "Import/export MD5 mismatch (type={}, width={})"
                             .format(st, sw))

    def test_frag_export_bytes(self):
        """Fragment.export_bytes"""
        duration = 0.1
        sample_type = splat.SAMPLE_INT
        sample_width = 16
        frag = splat.data.Fragment(duration=duration, channels=1)
        length = len(frag)
        for i in range(length):
            frag[i] = (float(i) / length,)
        x = list(int(length * r) for r in (0.0, 0.12, 0.34, 0.62, 0.987, 0.99))
        y = list((float(i) / length) for i in x)
        frag_bytes = frag.export_bytes(sample_type, sample_width)
        max16 = (2 ** 15) - 1
        for i, j in zip(x, y):
            val = frag_bytes[i * 2] + (frag_bytes[(i * 2) + 1] * 256)
            ref = int(j * max16)
            self.assertEqual(
                val, ref,
                "Incorrect 16-bit value: {} instead of {}".format(val, ref))
        i0, i1 = (int(length * r) for r in (0.23, 0.83))
        frag_bytes = frag.export_bytes(sample_type, sample_width, i0, i1)
        self.assertEqual(len(frag_bytes), ((i1 - i0) * sample_width / 8),
                         "Incorrect length of bytes array")
        for i, j in zip(x, y):
            if (i < i0) or (i > i1):
                continue
            k = i - i0
            val = frag_bytes[k * 2] + (frag_bytes[(k * 2) + 1] * 256)
            ref = int(j * max16)
            self.assertEqual(
                val, ref,
                "Incorrect 16 bits value: {} instead of {}".format(val, ref))
        start = frag.s2n(-2.3)
        end = frag.s2n(duration * 1.2)
        frag_bytes = frag.export_bytes(sample_type, sample_width,
                                       frag.s2n(-2.3),
                                       frag.s2n(duration * 1.2))
        b_len = len(frag_bytes)
        ref_len = len(frag) * sample_width / 8
        self.assertEqual(
            b_len, ref_len,
            "Incorrect data length: {} instead of {}".format(b_len, ref_len))

    def test_frag_dup(self):
        """Fragment.dup"""
        frag1 = splat.data.Fragment(channels=1)
        splat.gen.SineGenerator(frag=frag1).run(0.0, 0.93, 1234.0)
        frag2 = frag1.dup()
        self.assertEqual(frag1.md5(), frag2.md5(),
                         "Duplicated fragment MD5 mismatch")

    def test_frag_save(self):
        """Fragment.save"""
        duration = 0.1
        frag = splat.data.Fragment(duration=duration, channels=1)
        length = len(frag)
        for i in range(length):
            frag[i] = (float(i) * 0.9 / length,)
        for fmt in ['wav', 'saf']:
            f = StringIO()
            frag.save(f, fmt)
            f.reset()
            frag2 = splat.data.Fragment.open(f, fmt)
            if fmt == 'saf':
                self.assertEqual(frag.md5(), frag2.md5())
            for i in range(length):
                a = frag[i][0]
                b = frag2[i][0]
                self.assertAlmostEqual(
                    a, b, 3,
                    "Fragment data mismatch [{}] {} {}".format(i, a, b))


class InterpolTest(SplatTest):

    def test_polynomial(self):
        """interpol.Polynomial"""
        k0, k1, k2, k3 = coefs = (2.345, 3.6, 6.5, 100)
        p = splat.interpol.Polynomial(coefs)
        self.assertEqual(coefs, p.coefs,"Polynomial coefs mismatch")
        d = p.derivative()
        dcoefs = (k1, (k2 * 2), (k3 * 3))
        self.assertEqual(
            d.coefs, dcoefs,
            "Derivative error: {} instead of {}".format(d.coefs, dcoefs))
        i = d.integral(k0)
        self.assertEqual(
            i.coefs, coefs,
            "Integral error: {} instead of {}".format(i.coefs, coefs))

    def test_spline(self):
        """interpol.Spline"""
        pts = [(1.45, 6.78), (2.56, -3.56), (12.65, 9.6), (18.9, 5.0)]
        n = len(pts) - 1
        s0 = splat.interpol.Spline(pts)
        s1 = splat.interpol.Spline(pts, n=n)
        for p in pts:
            x, y = p[0], p[1]
            for s in (s0, s1):
                y0 = s.value(x)
                self.assertAlmostEqual(
                    y, y0, 6,
                    msg="Spline error, order: {}, s[{}] = {} != {}".format(
                        s._n, x, y0, y))


class SignalTest(SplatTest):

    def test_signal(self):
        """Signal"""
        duration = 0.0123
        x1 = int(duration * 0.234)
        x2 = int(duration * 0.789)
        frag = splat.data.Fragment(duration=duration)
        float_value = 1.234
        float_tuple = (float_value,)
        sig_float = splat.Signal(frag, float_value)
        self.assertEqual(
            len(sig_float), len(frag),
            "Signal and Fragment lengths mismatch: {} {}".format(
                len(sig_float), len(frag)))
        for i, s in enumerate(splat.Signal(frag, float_value)):
            self.assertIsInstance(s, tuple)
            self.assertEqual(len(s), 1)
            self.assertAlmostEqual(
                s[0], float_value, self._places,
                "Incorrect float signal value[{}]: {} {}".format(
                    i, s, float_tuple))
            for x in [x1, x2]:
                sig = sig_float[x]
                self.assertIsInstance(sig, tuple)
                self.assertEqual(len(sig), 1)
                self.assertAlmostEqual(
                    sig[0], float_value, self._places,
                    "Incorrect float signal indexed values: [{}] {} {}".format(
                        x, sig, float_tuple))
        func = lambda x: x * 0.1469
        (y1, y2) = (func(x) for x in (x1, x2))
        for i, (y,) in enumerate(splat.Signal(frag, func)):
            x = func(frag.n2s(i))
            self.assertAlmostEqual(
                x, y, self._places,
                "Incorrect function signal value[{}]: {} {}".format(i, x, y))
        frag2 = splat.data.Fragment(duration=duration, channels=1)
        splat.sources.sine(frag2, 0.0, 456.789)
        for i, (x, y) in enumerate(zip(frag2, splat.Signal(frag, frag2))):
            self.assertEqual(
                x, y, "Incorrect fragment signal value[{}]: {} {}".format(
                    i, x, y))
        for i, (y, z) in enumerate(splat.Signal(frag, (func, float_value))):
            x = func(frag.n2s(i))
            self.assertAlmostEqual(
                x, y, self._places,
                "Incorrect mixed signal value (func)[{}]: {} {}".format(
                    i, x, y))
            self.assertAlmostEqual(
                z, float_value, self._places,
                "Incorrect mixed signal value (float)[{}]: {} {}".format(
                    i, z, float_value))


class GeneratorTest(SplatTest):

    def test_gen_frag(self):
        """Generator Fragment"""
        gen = splat.gen.SineGenerator()
        self.assertTrue(isinstance(gen.frag, splat.data.Fragment),
                        "Invalid Fragment object")
        self.assertEqual(gen.frag.duration, 0.0,
                         "Incorrect initial Fragment length")

    def test_sine(self):
        """sources.sine"""
        freq = 1237.9
        frag_float = splat.data.Fragment(duration=1.0)
        splat.sources.sine(frag_float, -0.5, freq)
        frag_signal = splat.data.Fragment(duration=1.0)
        splat.sources.sine(frag_signal, -0.5, freq, lambda x: 0.0)
        frag_freq = splat.data.Fragment(duration=1.0, channels=1)
        frag_freq.offset(freq)
        frag_frag = splat.data.Fragment(duration=1.0)
        splat.sources.sine(frag_frag, -0.5, frag_freq)
        self.assert_md5([frag_float, frag_signal, frag_frag],
                        '46a8962a759033371f45c4ade9f2bfbd')

    def test_sine_gen(self):
        """gen.SineGenerator"""
        gen = splat.gen.SineGenerator()
        f = 1000.0
        gen.run(0.0, 1.0, f)
        n = int(0.1234 * gen.frag.duration * gen.frag.rate)
        s = math.sin(2 * math.pi * f * float(n) / gen.frag.rate)
        self.assert_samples(gen.frag, {n: (s, s)})
        self.assert_md5(gen.frag, 'ec18389e198ee868d61c9439343a3337')

    def test_square(self):
        """sources.square"""
        freq = 1237.9
        frag_float = splat.data.Fragment(duration=1.0)
        splat.sources.square(frag_float, -0.5, freq)
        frag_signal = splat.data.Fragment(duration=1.0)
        splat.sources.square(frag_signal, -0.5, freq, lambda x: 0.0)
        frag_freq = splat.data.Fragment(duration=1.0, channels=1)
        frag_freq.offset(freq)
        frag_frag = splat.data.Fragment(duration=1.0)
        splat.sources.square(frag_frag, -0.5, frag_freq)
        self.assert_md5([frag_float, frag_signal, frag_frag],
                        '6a6ab2e991baf48a6fe2c1d18700e40e')

    def test_square_gen(self):
        """gen.SquareGenerator"""
        gen = splat.gen.SquareGenerator()
        f = 1000.0
        gen.run(0.0, 1.0, f)
        nf = gen.frag.rate / f
        samples = {int(nf * 0.1): (1.0, 1.0), int(nf * 0.9): (-1.0, -1.0)}
        self.assert_samples(gen.frag, samples)
        self.assert_md5(gen.frag, '0ca047e998f512280800012b05107c63')

    def test_triangle(self):
        """sources.triangle"""
        freq = 1237.5
        frag_float = splat.data.Fragment(duration=1.0)
        splat.sources.triangle(frag_float, -0.5, freq)
        frag_signal = splat.data.Fragment(duration=1.0)
        splat.sources.triangle(frag_signal, -0.5, freq, lambda x: 0.0)
        frag_freq = splat.data.Fragment(duration=1.0, channels=1)
        frag_freq.offset(freq)
        frag_frag = splat.data.Fragment(duration=1.0)
        splat.sources.triangle(frag_frag, -0.5, frag_freq)
        self.assert_md5([frag_float, frag_signal, frag_frag],
                        '4bce3885732ba2f5450e79e42155adaa')

    def test_triangle_gen(self):
        """gen.TriangleGenerator"""
        gen = splat.gen.TriangleGenerator()
        f = 1000.0
        ratio = 0.567
        gen.run(0.0, 1.0, f, 0.0, ratio, levels=0.0)
        nf = gen.frag.rate / f
        x1 = 0.25
        t1 = int(nf * ratio * x1)
        s1 = (t1 * 2.0 / (ratio * nf)) - 1.0
        x2 = 0.75
        ratio2 = 1 - ratio
        t2 = int(nf * (ratio + (ratio2 * x2)))
        a2 = -2.0 / (ratio2 * nf)
        b2 = 1.0 - (a2 * ratio * nf)
        s2 = (t2 * a2) + b2
        samples = {t1: (s1, s1), t2: (s2, s2)}
        self.assert_samples(gen.frag, samples)
        self.assert_md5(gen.frag, 'b6d9eb000b328134cd500173b24f1c88')

    def test_overtones(self):
        """sources.overtones"""
        freq = 1237.5
        ot = [(1.3, 0.0, -2.5), (5.7, 10.0, -12.9)]
        frag_float = splat.data.Fragment(duration=1.0)
        splat.sources.overtones(frag_float, -0.5, freq, ot)
        frag_mixed = splat.data.Fragment(duration=1.0)
        splat.sources.overtones(frag_mixed, -0.5, freq, ot, lambda x: 0.0)
        frag_signal = splat.data.Fragment(duration=1.0)
        ot_sig = [(1.3, 0.0, -2.5), (5.7, lambda x: 10.0, -12.9)]
        splat.sources.overtones(frag_signal, -0.5, freq, ot_sig, lambda x: 0.0)
        frag_freq = splat.data.Fragment(duration=1.0, channels=1)
        frag_freq.offset(1237.5)
        frag_frag = splat.data.Fragment(duration=1.0)
        splat.sources.overtones(frag_frag, -0.5, frag_freq, ot)
        self.assert_md5([frag_float, frag_mixed, frag_signal, frag_frag],
                        '8974a1eea0db97af1aa171f531685e9d')

    def test_overtones_gen(self):
        """gen.OvertonesGenerator"""
        gen = splat.gen.OvertonesGenerator()
        gen.ot_decexp(1.0)
        gen.run(0.0, 1.0, 1000.0)
        self.assert_md5(gen.frag, 'ee045e012673ff7ed4ab9bd590b57368')


class ParticleTest(SplatTest):

    def test_particle(self):
        """gen.Particle"""
        start = 1.2
        end = 3.4
        f = 1234.56
        p = splat.gen.Particle(start, end, splat.lin2dB(f))
        self.assertEqual(p.start, start, "Incorrect particle start time")
        self.assertEqual(p.end, end, "Incorrect particle end time")
        self.assertAlmostEqual(
            p.freq, f, 6,
            msg="Frequency conversion error: {} instead of {}".format(
                f, p.freq))

    def test_particle_pool(self):
        """gen.ParticlePool"""
        min_f = 123.45
        min_f_log = splat.lin2dB(min_f)
        max_f = 678.9
        max_f_log = splat.lin2dB(max_f)
        min_len = 0.1
        max_len = 0.3
        envelope = splat.interpol.Spline([(0.5, 0.0), (1.3, 1.0), (2.1, 0.0)])
        n_slices = 20
        density = 200
        count = 196
        half_count = 103
        pool = splat.gen.ParticlePool(min_f_log, max_f_log, min_len, max_len,
                                      envelope, n_slices, density)
        self.assertEqual(pool.start, envelope.start,
                         "Incorrect pool start time")
        self.assertEqual(pool.end, envelope.end, "Incorrect pool end time")
        self.assertEqual(
            pool.count(), count, # plain magic
            "Unexpected number of particles: {}".format(pool.count()))
        for p in pool.iterate(share=0.5):
            self.assertFalse((p.length < min_len) or (p.length > max_len),
                             "Invalid particle length")
            self.assertFalse(
                (p.start < envelope.start) or (p.start > envelope.end),
                "Invalid particle start/end times")
            self.assertFalse((p.freq < min_f) or (p.freq > max_f),
                             "Invalid particle frequency")
        self.assertFalse(
            abs(pool.count() - half_count) > (count / 35.0),
            "Invalid number of particles left: {}".format(pool.count()))

    def test_particle_gen(self):
        """gen.ParticleGenerator"""

        class TestGenerator(splat.gen.Generator):
            def __init__(self, test, start, end, f_min, f_max, *args, **kw):
                super(TestGenerator, self).__init__(*args, **kw)
                self._test = test
                self._start = start
                self._end = end
                self._f_min = f_min
                self._f_max = f_max

            def run(self, start, end, freq, levels):
                self._test.assertFalse(
                    (start < self._start) or (end > self._end),
                    "Invalid subgen start/end run times")
                self._test.assertFalse(
                    (freq < self._f_min) or (freq > self._f_max),
                    "Invalid subgen frequency")

        z_start = 0.5
        z_mid = 1.8
        z_end = 2.5
        eq_min_f = 100.0
        eq_mid_f = 2345.0
        eq_max_f = 8000.0
        subgen = TestGenerator(self, z_start, z_end, eq_min_f, eq_max_f)
        pgen = splat.gen.ParticleGenerator(subgen)
        self.assertIs(pgen.subgen, subgen,
                      "Failed to get ParticleGenerator.subgen")
        pgen.set_z([(z_start, 0.0), (z_mid, 1.0, 0.0), (z_end, 0.0)])
        self.assertEqual(pgen.z.start, z_start,"Incorrect envelope start time")
        self.assertEqual(pgen.z.end, z_end, "Incorrect envelope end time")
        pgen.set_eq([(eq_min_f, -30.0), (eq_mid_f, 10.0), (eq_max_f, -50)])
        self.assertEqual(pgen.eq.start, eq_min_f, "Incorrect eq start freq")
        self.assertEqual(pgen.eq.end, eq_max_f, "Incorrect eq end freq")
        gfuzz = (3.0, 1.0)
        pgen.gain_fuzz = gfuzz
        self.assertEqual(pgen.gain_fuzz, gfuzz,"Failed to get gain fuzz")
        pgen.make_pool()
        self.assertEqual(pgen.pool.count(), 126, # plain magic again
                         "Incorrect number of particles")
        self.assertAlmostEqual(pgen.curve(123.45, 456.78, 56.78), 154.33118, 6,
                               "ParticleGenerator.curve error")
        pgen.do_show_progress = False
        pgen.run(0.0, 5.0, 1234.56)
        self.assertEqual(pgen.pool.count(), 0,
                         "Not all particles were consumed in the run")

# -----------------------------------------------------------------------------
# main function

if __name__ == '__main__':
    print("Sample precision: {}-bit".format(splat.NATIVE_SAMPLE_WIDTH))
    unittest.main()
