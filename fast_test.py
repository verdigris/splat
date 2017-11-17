# Splat - fast_test.py
#
# Copyright (C) 2015 Guillaume Tucker <guillaume@mangoz.org>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import math
import cmath
import argparse
import splat.gen
import splat.data
import splat
from splat import dB2lin as dB
import splat.tools.compare

def run_tests(freq=123.45, duration=30.0, phase=23.456, pts=None,
              overtones=None, verbose=False):
    if pts is None:
        pts = [0.0, 0.2, 0.234, 0.456, 0.45602, 0.7124, 0.89, 0.90001]
    if overtones is None:
        overtones = [(1.0, 0.0, 0.45), (2.58, 12.345, 0.45), (200.0, 0.0, 1.0)]

    cases = []

    # -------------------------------------------------------------------------

    frag = splat.data.Fragment(channels=1)
    gen = splat.gen.SineGenerator(frag=frag)
    gen.run(0.0, duration, freq, phase)
    frag.save('sine-{}.wav'.format(splat.SAMPLE_WIDTH), normalize=False)
    cases.append('sine')

    if verbose:
        print('sine')
    for t in pts:
        t = duration * t
        n = frag.s2n(t)
        t = float(n) / frag.rate
        st = 2 * cmath.pi * freq * (t + phase)
        y1 = frag[n][0]
        y2 = math.sin(st)
        delta = abs(y2 - y1)
        if delta != 0.0:
            delta_dB = "{:.3f}".format(splat.lin2dB(delta))
        else:
            delta_dB = 'infinity'
        if verbose:
            print(t, y1, y2, delta_dB)

    # -------------------------------------------------------------------------

    frag = splat.data.Fragment(channels=1)
    gen = splat.gen.OvertonesGenerator(frag=frag)
    gen.overtones = overtones
    gen.run(0.0, duration, freq)
    frag.save('overtones-{}.wav'.format(splat.SAMPLE_WIDTH), normalize=False)
    cases.append('overtones')

    max_ratio = frag.rate / 2.0 / freq
    ot_clipped = []
    for ot in gen.overtones:
        if ot[0] < max_ratio:
            ot_clipped.append(ot)

    if verbose:
        print('overtones')
    for t in pts:
        t = duration * t
        n = frag.s2n(t)
        y1 = frag[n][0]
        y2 = sum(a * math.sin(2 * cmath.pi * freq * r * (t + ph))
                 for r, ph, a in ot_clipped)
        delta = abs(y2 - y1)
        if delta != 0.0:
            delta_dB = "{:.3f}".format(splat.lin2dB(delta))
        else:
            delta_dB = 'infinity'
        if verbose:
            print(t, y1, y2, delta_dB)

    # -------------------------------------------------------------------------

    duration = 1.0
    sig = splat.data.Fragment(channels=1)
    splat.gen.TriangleGenerator(frag=sig).run(0.0, duration, 12.0, levels=0.5)
    sig.offset(0.5)

    mod = splat.data.Fragment(channels=1)
    splat.gen.TriangleGenerator(frag=mod).run(0.0, duration, 4.0, levels=0.002)
    mod.offset(0.5)

    frag = splat.data.Fragment()
    splat.gen.SineGenerator(frag=frag).run(0.0, duration, 456.0, levels=sig,
                                           phase=lambda x: math.sin(x))
    frag.save('sine-signal-{}.wav'.format(splat.SAMPLE_WIDTH), normalize=False)
    cases.append('sine-signal')

    frag = splat.data.Fragment()
    gen = splat.gen.OvertonesGenerator(frag=frag)
    gen.overtones = overtones
    gen.run(0.0, duration, 456.0, levels=sig)
    frag.save('overtones-mixed1-{}.wav'.format(splat.SAMPLE_WIDTH),
              normalize=False)
    cases.append('overtones-mixed1')

    frag = splat.data.Fragment()
    gen = splat.gen.OvertonesGenerator(frag=frag)
    gen.overtones = overtones
    gen.run(0.0, duration, 456.0, levels=sig, phase=mod)
    frag.save('overtones-mixed2-{}.wav'.format(splat.SAMPLE_WIDTH),
              normalize=False)
    cases.append('overtones-mixed2')

    sig2 = splat.data.Fragment(channels=1)
    splat.gen.TriangleGenerator(frag=sig2).run(0.0, duration, 1.8, levels=0.1)
    sig2.offset(-sig2.get_peak()[0]['min'])
    overtones.append((0.54, 0.0, sig2))

    frag = splat.data.Fragment()
    gen = splat.gen.OvertonesGenerator(frag=frag)
    gen.overtones = overtones
    gen.run(0.0, duration, 456.0, levels=sig, phase=mod)
    frag.save('overtones-signal-{}.wav'.format(splat.SAMPLE_WIDTH),
              normalize=False)
    cases.append('overtones-signal')

    frag = splat.data.Fragment()
    gen = splat.gen.SineGenerator(frag=frag)
    gen.run(0.0, 5.678, 1234.56, levels=dB(-3))
    ratio = 1.987
    new_len = int(len(frag) * ratio)
    rem = new_len % 4
    if rem:
        new_len -= rem
    else:
        new_len -= 4
    frag.resample(ratio=ratio)
    frag.resize(length=new_len)
    frag.save('resample-float-{}.wav'.format(splat.SAMPLE_WIDTH))
    cases.append('resample-float')

    frag = splat.data.Fragment()
    gen = splat.gen.SineGenerator(frag=frag)
    gen.run(0.0, 5.678, 1234.56, levels=dB(-3))
    ratio = 1.987
    frag.resample(ratio=lambda x: ratio)
    frag.save('resample-signal-{}.wav'.format(splat.SAMPLE_WIDTH))
    cases.append('resample-signal')

    return cases

def compare_all(cases, thr_dB=-40.0):
    ret = True
    for c in cases:
        f1, f2 = ('-'.join([c, str(w)]) + '.wav' for w in (64, 32))
        delta_dB = splat.tools.compare.file_peak_delta_dB(f1, f2)
        if delta_dB < thr_dB:
            res = 'OK'
        else:
            ret = False
            res = 'ERROR'
        print('{:16s} {:.3f} dB {}'.format(c, delta_dB, res))
    return ret

def main(argv):
    parser = argparse.ArgumentParser("Test fast Splat")
    parser.add_argument('--compare', action='store_true',
                        help="Compare with existing 64-bit reference files")
    parser.add_argument('--verbose', action='store_true',
                        help="Print more information")
    args = parser.parse_args(argv[1:])

    print("Sample precision: {}-bit".format(splat.SAMPLE_WIDTH))
    cases = run_tests(verbose=args.verbose)
    if args.compare is True:
        compare_all(cases)

    return True

if __name__ == '__main__':
    res = main(sys.argv)
    sys.exit(0 if res is True else 1)
