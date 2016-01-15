# Splat - splat/genlib.py
#
# Copyright (C) 2016 Guillaume Tucker <guillaume@mangoz.org>
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

import data
import sources
import gen
import filters
import interpol
from _splat import dB2lin as dB

def qwirky(tempo, rate, wobble=5, brush=5, bells=5):
    """``splat.genlib.qwirky`` is a funny generator for cheesy chords and
    silly melodies.

    wobble
      some kind of tremolo

    brush
      some kind of hiss

    bells
      some kind of chimes
    """
    T = 60.0 / tempo
    sig_len = T * 6
    envelope_spline = interpol.spline([
            (0.0, dB(-24)), (T/50, dB(0), 0.0), (T/10, dB(-6), -0.5),
            (T/3, dB(-9), -0.1), (sig_len, dB(-100.0))])
    envelope_frag = data.Fragment(rate=rate, channels=1, duration=sig_len)
    envelope_frag.offset(1.0)
    envelope_frag.amp(envelope_spline.signal)
    gen_filters = [(filters.linear_fade, (0.005,)),
                   (lambda frag, e: frag.amp(e), (envelope_frag,))]

    def ph(a, k):
        sig_frag = data.Fragment(rate=rate, channels=1, duration=sig_len)
        sources.sine(sig_frag, a/1000.0, k)
        return sig_frag

    g = gen.OvertonesGenerator()
    g.overtones = list()
    if wobble > 0:
        w0 = (wobble - 5) * 0.2
        g.overtones = [
            (1.0, ph(0.7 + w0, T*6), dB(0.0)),
            (2.0, ph(0.4 + w0/2, T*6.2), dB(-12.0 + (wobble - 5) * 6/5)),
            (4.0, ph(0.2 + w0/4, T*12.6), dB(-18.0 + (wobble - 5) * 9/5)),
            (8.0, ph(0.1 + w0/8, T*18.8), dB(-40.0 + (wobble - 5) * 20/5)),
            ]
    if brush > 0:
        dB0 = -24.0 + (brush - 5.0) * 9.0
        dBx = 2.1 - (brush - 5.0) * 0.01
        g.overtones += list(
            ((float(x), 0.0, dB(dB0 - dBx * x))
             for x in range(3, 48, 3)))
    if bells > 0:
        dB0 = -18.0 + (bells - 5) * 1.5
        dBx = 3.0 - (bells - 5) * 0.25
        g.overtones += list(
            ((float(x), 0.0, dB(dB0 - x*dBx))
             for x in range(5, 35, 5)))
    g.filters = gen_filters

    return g
