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

def groovy(tempo, rate, fade=0.03, body=5, buzz=5, rattle=5):
    """``splat.genlib.groovy`` is good for punchy bass lines

    fade
      time used in the ``linear_fade`` filter, to round the corners

    body
      overall thickness

    buzz
      presence in the higher frequencies

    rattle
      texture in the medium frequencies
    """
    T = 60.0 / tempo
    duration = T*6.0
    envelope_spline = interpol.spline([
            (0.0, dB(0.0)), (T/4, dB(2.5), 0.0),
            (T/2, dB(-2.0), 0.0), (duration, dB(-100.0))])
    envelope_frag = data.Fragment(rate=rate, channels=1, duration=duration)
    envelope_frag.offset(1.0)
    envelope_frag.amp(envelope_spline.signal)
    gen_filters = [(filters.linear_fade, (fade,)),
                   (lambda frag, e: frag.amp(e), (envelope_frag,))]

    gen_tri = gen.TriangleGenerator()
    gen_tri.filters = gen_filters
    tri_dB = (-4.0 + (buzz - 5) * 4) if buzz > 0 else None

    gen_ot = gen.OvertonesGenerator()
    gen_ot.overtones = list()
    if body > 0:
        gen_ot.ot_decexp(1.3 + (float(body) - 5) / 6, 24)
    if rattle > 0:
        dB0 = -18.0 + (rattle - 5) * 0.5
        dBx = 1.15 - (rattle - 5) * 0.2
        gen_ot.overtones += \
            ((n * 1.01, 0.0, dB(dB0 - (n * dBx))) for n in range(18, 48))
    gen_ot.filters = gen_filters

    class GroovyGen(gen.Generator):

        def __init__(self, gen_tri, gen_ot, tri_dB, *args, **kw):
            super(GroovyGen, self).__init__(*args, **kw)
            self._gen_tri = gen_tri
            self._gen_ot = gen_ot
            self._tri_dB = tri_dB

        def run(self, start, end, freq, levels):
            self._gen_ot.frag = self.frag
            self._gen_ot.run(start, end, freq, levels=levels)
            self._gen_tri.frag = self.frag
            if self._tri_dB is not None:
                lvl = dB(levels + self._tri_dB)
                self._gen_tri.run(start, end, freq*2, levels=lvl)

    return GroovyGen(gen_tri, gen_ot, tri_dB)
