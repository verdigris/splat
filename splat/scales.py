# Splat - splat/scales.py
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

import math

class Scale(object):
    note_degree = { 'Ab': 11, 'A': 0,  'A#': 1,
                    'Bb': 1,  'B': 2,
                              'C': 3,  'C#': 4,
                    'Db': 4,  'D': 5,  'D#': 6,
                    'Eb': 6,  'E': 7,
                              'F': 8,  'F#': 9,
                    'Gb': 9,  'G': 10, 'G#': 11 }

    def __init__(self, base=2, fund=440.0, key='A', steps=12, *args, **kw):
        self._base = base
        self._fund = fund
        self._key = key
        self._steps = steps
        self._degrees = tuple(sorted(self._init_degrees(*args, **kw)))
        for d in self._degrees:
            assert((d >= 0) and (d <= self.base))

    def _init_degrees(self):
        raise NotImplementedError

    @property
    def base(self):
        return self._base

    @property
    def length(self):
        return len(self._degrees)

    @property
    def degrees(self):
        return self._degrees

    @property
    def fundamental(self):
        return self._fund

    @fundamental.setter
    def fundamental(self, value):
        self._fund = float(value)

    @property
    def key(self):
        return self._key

    def get_freq(self, degree, cycle):
        f0 = self.fundamental * self.degrees[degree]
        return f0 * math.pow(self.base, cycle)

    def get_note(self, note):
        if (len(note) > 1) and (note[1] in ['#', 'b']):
            n = 2
        else:
            n = 1
        degree = Scale.note_degree[note[:n]] - Scale.note_degree[self.key]
        if len(note) > n:
            cycle = int(note[n:])
        else:
            cycle = 0
        if degree < 0:
            degree += 12
        return degree, cycle

    def get_note_freq(self, note):
        degree, cycle = self.get_note(note)
        return self.get_freq(degree, cycle)

    def __getitem__(self, note):
        return self.get_note_freq(note)


class LogScale(Scale):
    def _init_degrees(self):
        degrees = list()
        for n in range(self._steps):
            degrees.append(math.pow(self.base, (float(n) / self._steps)))
        return degrees


class HarmonicScale(Scale):
    def _init_degrees(self, k=3):
        degrees = list()
        for n in range(self._steps):
            f = math.pow(k, n)
            f = math.pow(self.base, math.log(f, self.base) % 1)
            degrees.append(f)
        return degrees
