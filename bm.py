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
import argparse
import benchmark
import copy
import cPickle
import _splat
import splat
import splat.data
import splat.gen
import splat.sources
import splat.interpol
from splat import dB2lin as dB

DEFAULT_ITERATIONS = 15
PREFIX = 'run_test_'

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
        splat.sources.sine(self.mod, dB(-9.0), 1.0)
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
        tmp_gen = splat.gen.OvertonesGenerator()
        tmp_gen.ot_decexp(n=6)
        self.ot = tmp_gen.overtones
        self.sig_ot = copy.copy(self.ot)
        mod = splat.data.Fragment(channels=1)
        splat.gen.TriangleGenerator(frag=mod).run(
            0.0, self.frag.duration, 6.5, 12.0)
        self.sig_ot += [(2.3, 0.0, mod)]

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
        self.frag.offset(splat.interpol.spline(self.pts).value)

    def test_offset_spline_signal(self):
        self.frag.offset(splat.interpol.spline(self.pts).signal)


def run_no_report(iterations=DEFAULT_ITERATIONS):
    """
    Run the benchmark but return the result objects and do not print a report.
    """
    bms = []
    obj = None
    for obj in globals().itervalues():
        if isinstance(obj, type) and issubclass(obj, benchmark.Benchmark):
            bms.append(obj)

    runs = 0
    results = []
    for bm_cls in bms:
        bm = bm_cls(each=iterations, prefix=PREFIX)
        bm.run()
        results.append(bm)
        runs += bm.getTotalRuns()
    return {res.__class__.__name__: res.table for res in results}


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Standard Splat benchmark.")
    parser.add_argument('--iterations', type=int, default=DEFAULT_ITERATIONS,
                        help="number of iterations")
    parser.add_argument('--format', default='reST',
                        choices=['reST', 'csv', 'comma', 'markdown'],
                        help="report output format")
    parser.add_argument('--pickle',
                        help="do not generate report but pickle the results")
    args = parser.parse_args(sys.argv[1:])

    if args.pickle:
        results = run_no_report(iterations=args.iterations)
        cPickle.dump(results, open(args.pickle, 'w'))
    else:
        benchmark.main(format=args.format, numberFormat="%.4g",
                       each=args.iterations, prefix=PREFIX)

    sys.exit(0)
