# Splat - test.py
#
# Copyright (C) 2012, 2013, 2014, 2015 Guillaume Tucker <guillaume@mangoz.org>
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
import hashlib
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
import splat.scales
import splat.seq
from splat import dB2lin as dB
import splat.tools.compare

splat.check_version((1, 6))

# -----------------------------------------------------------------------------
# Splat test base class

class SplatTest(unittest.TestCase):

    def setUp(self):
        self._places = 10 if splat.SAMPLE_WIDTH == 64 else 4

    def assert_md5(self, frags, hexdigest):
        if isinstance(frags, splat.data.Fragment):
            frags = [frags]
        for frag in frags:
            md5sum = frag.md5(sample_type="int16")
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
        frag_name = "Test Fragment"
        frag = splat.data.Fragment(duration=1.0, name=frag_name)
        self.assertEqual(frag_name, frag.name)
        self.assert_samples(frag, {int(len(frag) / 2): (0.0, 0.0)})
        self.assert_md5(frag, 'fe384f668da282694c29a84ebd33481d')

    def test_frag_md5(self):
        """Fragment.md5"""
        frag = splat.data.Fragment()
        splat.gen.SineGenerator(frag=frag).run(0.0, 0.345, 123.0)
        for sample_type in splat.sample_types.iterkeys():
            md5sum = hashlib.md5(frag.export_bytes(sample_type)).hexdigest()
            self.assertEqual(md5sum, frag.md5(sample_type))

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

    def test_frag_resize(self):
        """Fragment.resize"""
        duration = 1.3
        rate = 24000
        length = int(rate * duration)
        frag = splat.data.Fragment(rate=rate)
        self.assertEqual(len(frag), 0)
        frag.resize(duration=duration)
        self.assertEqual(len(frag), length)
        frag.grow(length=(length / 2))
        self.assertEqual(len(frag), length)
        frag.resize(length=(length / 2))
        self.assertEqual(len(frag), length / 2)
        frag.grow(duration=(duration * 1.5))
        self.assertEqual(len(frag), (length * 1.5))

    def test_frag_normalize(self):
        """Fragment.normalize"""
        levels = dB(-3.0)
        places = int(self._places / 2)
        small_places = 2
        frag = splat.data.Fragment(channels=1)
        splat.gen.SineGenerator(frag=frag).run(0.0, 1.0, 123.4, levels=levels)
        sample_n = (len(frag) / 2)
        x = frag[sample_n][0]
        frag_peak = frag.get_peak()[0]
        peak, avg = (frag_peak[item] for item in ['peak', 'avg'])
        self.assertAlmostEqual(avg, 0.0, small_places)
        self.assertAlmostEqual(levels, peak, places)
        frag.normalize()
        frag_peak = frag.get_peak()[0]
        norm_peak, norm_avg = (frag_peak[item] for item in ['peak', 'avg'])
        ref_peak = splat.dB2lin(-0.05)
        self.assertAlmostEqual(norm_peak, ref_peak, places)
        self.assertAlmostEqual(norm_avg, 0.0, places)
        y = frag[sample_n][0]
        ref = x * ref_peak / levels
        self.assertAlmostEqual(y, ref, small_places)

    def test_frag_mix(self):
        """Fragment.mix"""
        duration = 1.0
        frag1, frag2 = (splat.data.Fragment(channels=1) for i in range(2))
        splat.gen.SineGenerator(frag=frag1).run(0.0, duration, 1234.0, -3.0)
        splat.gen.SineGenerator(frag=frag2).run(0.0, duration, 5678.0, -3.0)
        offset = 0.234
        offset_n = frag1.s2n(offset)
        ref1_n = int(offset_n / 2)
        ref2_n = int(len(frag2) / 2)
        sample = 0.123

        sample_n = frag1.s2n(sample)
        frag1x, frag2x = (frag.dup() for frag in (frag1, frag2))
        x1 = frag1x[ref1_n][0]
        x2 = frag2x[ref2_n][0]
        a = frag1x[offset_n + sample_n][0]
        b = frag2x[sample_n][0]
        frag1x.mix(frag2, offset=offset)
        c = frag1x[offset_n + sample_n][0]
        y1 = frag1x[ref1_n][0]
        y2 = frag2x[ref2_n][0]
        self.assertEqual(x1, y1)
        self.assertEqual(x2, y2)
        self.assertAlmostEqual((a + b), c, self._places)

        frag1x, frag2x = (frag.dup() for frag in (frag1, frag2))
        skip = 0.087
        skip_n = frag1.s2n(skip)
        a = frag1x[offset_n + sample_n][0]
        b = frag2x[skip_n + sample_n][0]
        x1 = frag1x[ref1_n][0]
        x2 = frag2x[ref2_n][0]
        frag1x.mix(frag2x, offset, skip)
        c = frag1x[offset_n + sample_n][0]
        self.assertEqual(x1, y1)
        self.assertEqual(x2, y2)
        self.assertAlmostEqual((a + b), c, self._places)

    def test_frag_import_bytes(self):
        """Fragment.import_bytes"""
        frag = splat.data.Fragment()
        splat.gen.SineGenerator(frag).run(0.1, 2.7, 3456.7)
        for sample_type, sample_width in splat.sample_types.iteritems():
            ref_bytes = frag.export_bytes(sample_type)
            ref_md5 = hashlib.md5(ref_bytes).hexdigest()
            imp = splat.data.Fragment()
            imp.import_bytes(ref_bytes, frag.rate, frag.channels, sample_type)
            exp_bytes = imp.export_bytes(sample_type)
            exp_md5 = hashlib.md5(exp_bytes).hexdigest()
            self.assertEqual(ref_md5, exp_md5,
                             "Import/export MD5 mismatch (type={}, width={})"
                             .format(sample_type, sample_width))

    def test_frag_export_bytes(self):
        """Fragment.export_bytes"""
        duration = 0.1
        sample_type = "int16"
        sample_width = 16
        frag = splat.data.Fragment(duration=duration, channels=1)
        length = len(frag)
        for i in range(length):
            frag[i] = (float(i) / length,)
        x = list(int(length * r) for r in (0.0, 0.12, 0.34, 0.62, 0.987, 0.99))
        y = list((float(i) / length) for i in x)
        frag_bytes = frag.export_bytes(sample_type)
        max16 = (2 ** 15) - 1
        for i, j in zip(x, y):
            val = frag_bytes[i * 2] + (frag_bytes[(i * 2) + 1] * 256)
            ref = int(j * max16)
            self.assertEqual(
                val, ref,
                "Incorrect 16-bit value: {} instead of {}".format(val, ref))
        i0, i1 = (int(length * r) for r in (0.23, 0.83))
        frag_bytes = frag.export_bytes(sample_type, i0, i1)
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
        frag_bytes = frag.export_bytes(sample_type, frag.s2n(-2.3),
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

    def test_frag_resample(self):
        """Fragment.resample"""
        def run_resample_test(source, duration, ratio, freq, rate, thr):
            frag1 = splat.data.Fragment(duration=duration, channels=1)
            source(frag1, 1.0, freq)
            for ratio_value in (ratio, lambda x: ratio):
                frag2 = frag1.dup()
                frag2.resample(ratio=ratio_value, rate=rate)
                new_length = int(len(frag1) * ratio * rate / frag1.rate)
                if isinstance(ratio_value, float):
                    self.assertEqual(len(frag2), new_length,
                                     "Incorrect resampled fragment length")
                    new_duration = float(new_length) / frag2.rate
                    self.assertAlmostEqual(frag2.duration, new_duration,
                                           self._places,
                                       "Incorrect resampled fragment duration")
                else:
                    new_length = min(len(frag1), new_length)
                    delta = abs(len(frag2) - new_length)
                    self.assertLessEqual(delta, 1,
                                  "Incorrect signal resampled fragment length")
                ref = splat.data.Fragment(length=len(frag2), channels=1,
                                          rate=rate)
                source(ref, 1.0, (freq / ratio))
                delta = splat.tools.compare.frag_delta(ref, frag2)
                err = splat.lin2dB(delta.get_peak()[0]['peak'])
                self.assertLess(err, thr,
                                "Interpolation error: {:.1f} dB".format(err))

        for src, thr in [(splat.sources.sine, -60.0),
                         (splat.sources.triangle, -27.0)]:
            for f in [25.0, 1287.6]:
                for d in [0.01, 0.8642]:
                    for ratio in [0.35, (9.0 / 7.0)]:
                        for rate in [48000, 44100, 96000]:
                            run_resample_test(src, d, ratio, f, rate, thr)


class FragmentTest_mmap(FragmentTest):

    def setUp(self):
        super(FragmentTest_mmap, self).setUp()
        splat.use_mmap(True)

    def tearDown(self):
        super(FragmentTest_mmap, self).tearDown()
        splat.use_mmap(False)


class InterpolTest(SplatTest):

    def test_polynomial(self):
        """interpol.Polynomial"""
        k0, k1, k2, k3 = coefs = (2.345, 3.6, 6.5, 100)
        p = splat.interpol.Polynomial(coefs)
        self.assertEqual(coefs, p.coefs, "Polynomial coefs mismatch")
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
        """interpol.spline"""
        pts = [(1.45, 6.78), (2.56, -3.56), (12.65, 9.6), (18.9, 5.0)]
        n0 = 2
        n1 = len(pts) - 1
        s0 = splat.interpol.spline(pts, n=n0)
        s1 = splat.interpol.spline(pts, n=n1)
        for p in pts:
            x, y = p[0], p[1]
            for s, n in [(s0, n0), (s1, n1)]:
                y0 = s.value(x)
                self.assertAlmostEqual(
                    y, y0, 6,
                    msg="spline error, order: {}, s[{}] = {} != {}".format(
                        n, x, y0, y))


class SignalTest(SplatTest):

    def test_signal_callable(self):
        """Callable signals"""
        offset_value = 0.85
        frag = splat.data.Fragment(channels=1, duration=1.0)
        frag.offset(lambda x: offset_value)
        pts = [0.0, 0.1, 0.4, 0.6, 0.9, 0.99]
        for pt in pts:
            self.assertEqual(frag[frag.s2n(pt)], (offset_value,))

    def test_signal_frag(self):
        """Fragment signals"""
        offset_value = 1.5
        sig_frag = splat.data.Fragment(channels=1, duration=1.0)
        sig_frag.offset(offset_value)
        frag = splat.data.Fragment(channels=1, duration=1.0)
        frag.offset(sig_frag)
        pts = [0.0, 0.1, 0.4, 0.6, 0.9, 0.99]
        for pt in pts:
            self.assertEqual(frag[frag.s2n(pt)], (offset_value,))

    def test_signal_spline(self):
        """Spline signals"""
        spline_pts = [(0.0, 0.0), (0.1, 0.5), (0.5, 0.2), (1.0, 1.0)]
        spline = splat.interpol.spline(spline_pts, 1.23)
        frag1 = splat.data.Fragment(channels=1, duration=1.0)
        frag1.offset(spline.value)
        frag2 = splat.data.Fragment(channels=1, duration=1.0)
        frag2.offset(spline.signal)
        frags = [frag1, frag2]
        pts = [0.0, 0.1, 0.4, 0.6, 0.9, 0.99]
        for pt in pts:
            offset_value = spline.value(pt)
            for frag in frags:
                self.assertAlmostEqual(frag[frag.s2n(pt)][0], offset_value,
                                       self._places)

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
        splat.sources.sine(frag2, 1.0, 456.789)
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

    def test_gen_levels(self):
        """Generator and source levels"""
        freq = 123.45
        duration = 1.0
        for levels in [-3.0, -3, int(-3),
                        long(-3), (-3.0, -3.0),
                        (int(-3), int(-3)), (long(-3), long(-3))]:
            frag = splat.data.Fragment(duration=duration)
            splat.sources.sine(frag, levels, freq)
            gen = splat.gen.SineGenerator()
            gen.run(0.0, duration, freq, levels=levels)
            self.assert_md5([frag, gen.frag],
                            '1d6c38f467b1bdb5a9c3e8d573e815f4')

    def test_sine(self):
        """sources.sine"""
        freq = 1237.9
        ph = 0.123
        lvl0 = dB(-0.5)
        frag_float = splat.data.Fragment(duration=1.0)
        splat.sources.sine(frag_float, lvl0, freq, ph)
        frag_signal = splat.data.Fragment(duration=1.0)
        splat.sources.sine(frag_signal, lvl0, freq, lambda x: ph)
        frag_freq = splat.data.Fragment(duration=1.0, channels=1)
        frag_freq.offset(freq)
        frag_frag = splat.data.Fragment(duration=1.0)
        splat.sources.sine(frag_frag, lvl0, frag_freq, ph)
        self.assert_md5([frag_float, frag_signal, frag_frag],
                        'ebd3117927861068cf77af5ed2e7c5d7')

    def test_sine_gen(self):
        """gen.SineGenerator"""
        gen = splat.gen.SineGenerator()
        freq = 1000.0
        ph = 6.789
        gen.run(0.0, 1.0, freq, ph)
        n = int(0.1234 * gen.frag.duration * gen.frag.rate)
        s = math.sin(2 * math.pi * freq * (ph + float(n) / gen.frag.rate))
        self.assert_samples(gen.frag, {n: (s, s)})
        self.assert_md5(gen.frag, 'ec18389e198ee868d61c9439343a3337')

    def test_square(self):
        """sources.square"""
        freq = 1237.9
        lvl0 = dB(-0.5)
        frag_float = splat.data.Fragment(duration=1.0)
        splat.sources.square(frag_float, lvl0, freq)
        frag_signal = splat.data.Fragment(duration=1.0)
        splat.sources.square(frag_signal, lvl0, freq, lambda x: 0.0)
        frag_freq = splat.data.Fragment(duration=1.0, channels=1)
        frag_freq.offset(freq)
        frag_frag = splat.data.Fragment(duration=1.0)
        splat.sources.square(frag_frag, lvl0, frag_freq)
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
        lvl0 = dB(-0.5)
        frag_float = splat.data.Fragment(duration=1.0)
        splat.sources.triangle(frag_float, lvl0, freq)
        frag_signal = splat.data.Fragment(duration=1.0)
        splat.sources.triangle(frag_signal, lvl0, freq, lambda x: 0.0)
        frag_freq = splat.data.Fragment(duration=1.0, channels=1)
        frag_freq.offset(freq)
        frag_frag = splat.data.Fragment(duration=1.0)
        splat.sources.triangle(frag_frag, lvl0, frag_freq)
        self.assert_md5([frag_float, frag_signal, frag_frag],
                        '4bce3885732ba2f5450e79e42155adaa')

    def test_triangle_gen(self):
        """gen.TriangleGenerator"""
        gen = splat.gen.TriangleGenerator()
        f = 1000.0
        ratio = 0.567
        gen.run(0.0, 1.0, f, 0.0, ratio)
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
        lvl0 = dB(-0.5)
        ot = [(1.3, 0.0, dB(-2.5)), (5.7, 10.0, dB(-12.9))]
        ot_sig = [(1.3, 0.0, dB(-2.5)), (5.7, lambda x: 10.0, dB(-12.9))]
        frag_float = splat.data.Fragment(duration=1.0)
        splat.sources.overtones(frag_float, lvl0, freq, ot)
        frag_mixed = splat.data.Fragment(duration=1.0)
        splat.sources.overtones(frag_mixed, lvl0, freq, ot, lambda x: 0.0)
        frag_signal = splat.data.Fragment(duration=1.0)
        splat.sources.overtones(frag_signal, lvl0, freq, ot_sig, lambda x: 0.0)
        frag_freq = splat.data.Fragment(duration=1.0, channels=1)
        frag_freq.offset(1237.5)
        frag_frag = splat.data.Fragment(duration=1.0)
        splat.sources.overtones(frag_frag, lvl0, frag_freq, ot)
        self.assert_md5([frag_float, frag_mixed, frag_signal, frag_frag],
                        '8974a1eea0db97af1aa171f531685e9d')

    def test_overtones_gen(self):
        """gen.OvertonesGenerator"""
        gen = splat.gen.OvertonesGenerator()
        gen.ot_decexp(1.0)
        gen.run(0.0, 1.0, 1000.0, 0.1)
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
        envelope = splat.interpol.spline([(0.5, 0.0), (1.3, 1.0), (2.1, 0.0)])
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


class ScaleTest(SplatTest):

    def _d2f(self, s):
        return lambda d, o: s.f0 * math.pow(2, (d  + (o * 12)) / 12.0)

    def _check_note_freqs(self, s, note_freqs):
        for note, freq in note_freqs:
            self.assertAlmostEqual(s[note], freq, 12,
                                   "Mismatch note '{}', {} != {}".format(
                    note, s[note], freq))

    def test_log_scale_default(self):
        s = splat.scales.LogScale()
        d2f = self._d2f(s)
        note_freqs_default = [
            ('A', d2f(0, 0)), ('A1', d2f(0, 1)), ('A-1', d2f(0, -1)),
            ('A2', d2f(0, 2)), ('D', d2f(5, 0)), ('E', d2f(7, 0)),
            ('E1', d2f(7, 1)), ('B', d2f(2, 0)), ('B3', d2f(2, 3)),
            ('F', d2f(8, 0)), ('G#', d2f(11, 0)), ('G#0', d2f(11, 0)),
            ]
        self._check_note_freqs(s, note_freqs_default)

    def test_log_scale_rel(self):
        s = splat.scales.LogScale(key='E')
        d2f = self._d2f(s)
        note_freqs_rel = [
            ('A', d2f(5, 0)), ('A1', d2f(5, 1)), ('A-1', d2f(5, -1)),
            ('A2', d2f(5, 2)), ('D', d2f(10, 0)), ('E', d2f(0, 0)),
            ('E1', d2f(0, 1)), ('B', d2f(7, 0)), ('B3', d2f(7, 3)),
            ('F', d2f(1, 0)), ('G#', d2f(4, 0)), ('G#0', d2f(4, 0)),
            ]
        self._check_note_freqs(s, note_freqs_rel)

    def test_log_scale_abs(self):
        s = splat.scales.LogScale(key='E')
        s.abs_cycles = True
        d2f = self._d2f(s)
        note_freqs_abs = [
            ('A', d2f(5, -1)), ('A4', d2f(5, -1)), ('A3', d2f(5, -2)),
            ('A5', d2f(5, 0)), ('D', d2f(10, -1)), ('E', d2f(0, 0)),
            ('E5', d2f(0, 1)), ('B', d2f(7, -1)), ('B7', d2f(7, 2)),
            ('F', d2f(1, 0)), ('G#', d2f(4, 0)), ('G#4', d2f(4, 0)),
            ]
        self._check_note_freqs(s, note_freqs_abs)

    def test_harmonic_scale(self):
        s = splat.scales.HarmonicScale()
        f0 = 440.0
        note_freqs = [
            ('A', f0), ('E', f0 * 3 / 2), ('B', f0 * 9 / 8),
            ('F#', f0 * 27 / 16), ('C#', (f0 * 81 / 64)),
            ('D', f0 * 4 / 3), ('G', f0 * 16 / 9),
            ('A1', f0 * 2), ('E1', f0 * 6 / 2), ('B1', f0 * 18 / 8),
            ('D1', f0 * 8 / 3), ('G1', f0 * 32 / 9),
            ('A-1', f0 / 2), ('E-1', f0 * 3 / 4), ('B-1', f0 * 9 / 16),
            ('D-1', f0 * 4 / 6), ('G-1', f0 * 16 / 18),
            ]
        self._check_note_freqs(s, note_freqs)


class SequencerTest(SplatTest):

    class TestPattern(splat.seq.Pattern):

        def __init__(self, gen, *args, **kw):
            super(SequencerTest.TestPattern, self).__init__(*args, **kw)
            self.gen = gen

        def play(self, frag, bar, beat, t, T):
            freq = 220.0 * (beat + 1)
            self.gen.frag = frag
            self.gen.run(t, t + (T * 0.8), freq)


    def test_pattern_sequencer(self):
        p = SequencerTest.TestPattern(splat.gen.SineGenerator())
        frag = splat.data.Fragment(channels=1)
        splat.seq.PatternSequencer(120).run(frag, [(p,), (p,)])
        self.assert_md5(frag, '3b9dbf48691e185339a461a493fc5e7e')

# -----------------------------------------------------------------------------
# main function

if __name__ == '__main__':
    print("Sample precision: {}-bit".format(splat.SAMPLE_WIDTH))
    unittest.main()
