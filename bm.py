# Splat - bm.py
#
# Copyright (C) 2014 Guillaume Tucker <guillaume@mangoz.org>
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
import benchmark
import _splat
import splat
import splat.data
import splat.gen
import splat.sources
import splat.interpol

class RefMixin(object):

    def setUp(self):
        self.frag = splat.data.Fragment(duration=90, channels=1)
        if getattr(self, 'skip', False) is False:
            for test in dir(self):
                if test.startswith('test_'):
                    test_meth = getattr(self, test)
                    setattr(self, '_'.join(['run', test]), test_meth)

    def test_ref(self):
        _splat.gen_ref(self.frag)


class GenMixin(RefMixin):

    def setUp(self):
        super(GenMixin, self).setUp()
        self.freq = 1234.0
        self.mod = splat.data.Fragment(duration=self.frag.duration, channels=1)
        splat.sources.sine(self.mod, -9.0, 1.0)
        self.mod.offset(0.5)

    def test_source(self, *args):
        self.source(self.frag, 0.0, self.freq, *args)

    def test_source_signal(self, *args):
        self.source(self.frag, self.mod, self.freq, *args)

    def test_gen(self, *args, **kw):
        self.gen_cls(frag=self.frag, *args, **kw).run(
            0.0, self.frag.duration, self.freq)

    def test_gen_signal(self, *args, **kw):
        self.gen_cls(frag=self.frag, *args, **kw).run(
            0.0, self.frag.duration, self.freq, levels=self.mod)


class Sine(benchmark.Benchmark, GenMixin):
    source = splat.sources.sine
    gen_cls = splat.gen.SineGenerator


class Triangle(benchmark.Benchmark, GenMixin):
    source = splat.sources.triangle
    gen_cls = splat.gen.TriangleGenerator


class Square(benchmark.Benchmark, GenMixin):
    source = splat.sources.square
    gen_cls = splat.gen.SquareGenerator


class Overtones(benchmark.Benchmark, GenMixin):
    source = splat.sources.overtones
    gen_cls = splat.gen.OvertonesGenerator

    def setUp(self):
        super(Overtones, self).setUp()
        self.ot = [(1.0, 0.0, 0.0), (2.5, 0.0, -12.0)]
        self.sig_ot = [(1.0, 0.0, lambda x: 0.0), (2.5, 0.0, lambda x: -12.0)]

    def test_source_single(self):
        super(Overtones, self).test_source([(1.0, 0.0, 0.0)])

    def test_source(self):
        super(Overtones, self).test_source(self.ot)

    def test_source_mixed(self):
        super(Overtones, self).test_source_signal(self.ot)

    def test_source_signal(self):
        super(Overtones, self).test_source_signal(self.sig_ot)

    def test_gen_single(self):
        super(Overtones, self).test_gen()

    def test_gen(self):
        gen = self.gen_cls(frag=self.frag)
        gen.overtones = self.ot
        gen.run(0.0, self.frag.duration, self.freq)

    def test_gen_mixed(self):
        gen = self.gen_cls(frag=self.frag)
        gen.overtones = self.ot
        gen.run(0.0, self.frag.duration, self.freq, levels=self.mod)

    def test_gen_signal(self):
        gen = self.gen_cls(frag=self.frag)
        gen.overtones = self.sig_ot
        gen.run(0.0, self.frag.duration, self.freq, levels=self.mod)


class Spline(benchmark.Benchmark, RefMixin):

    def setUp(self):
        super(Spline, self).setUp()
        step = self.frag.duration / 10
        t = 0.0
        self.pts = list()
        while t <= self.frag.duration:
            self.pts.append((t, t / self.frag.duration))
            t += step

    def test_ref(self):
        self.frag.offset(1.5)

    def test_offset_frag(self):
        amp = splat.data.Fragment(channels=1, length=len(self.frag))
        amp.offset(1.5)
        self.frag.offset(amp)

    def test_offset_spline(self):
        self.frag.offset(splat.interpol.spline(self.pts).signal)

    def test_offset_spline_value(self):
        self.frag.offset(splat.interpol.spline(self.pts).value)


if __name__ == '__main__':
    benchmark.main(format="reST", numberFormat="%.4g", each=15,
                   prefix='run_test_')
    sys.exit(0)
