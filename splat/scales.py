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
    """Musical scale note frequency calculator.

    For a given scale definition, a Scale object will make each note correspond
    to a frequency value which can then be typically used with a sound
    generator.  It can be used as a read-only dictionary producing a frequency
    for a given note name.
    """

    _note_step = { 'Ab': 11, 'A': 0,  'A#': 1,
                   'Bb': 1,  'B': 2,
                             'C': 3,  'C#': 4,
                   'Db': 4,  'D': 5,  'D#': 6,
                   'Eb': 6,  'E': 7,
                             'F': 8,  'F#': 9,
                   'Gb': 9,  'G': 10, 'G#': 11 }

    def __init__(self, key='A', pitch=440.0, steps=12, base=2, *args, **kw):
        """Create a scale with a given ``key`` name which defaults to 'A', a
        ``pitch`` frequency which defaults to 440Hz, a number of ``steps``
        which defaults to 12 for the classic 12 semitones scale, and a ``base``
        which defines the cycle coefficient of the scale, by default 2 for
        classic octaves.
        """
        self._pitch = pitch
        self._key = key
        self._steps = steps
        self._base = base
        self._abs_cycles = False
        self._f0 = None

    @property
    def key(self):
        """Key of the scale (first note name)."""
        return self._key

    @property
    def pitch(self):
        """Pitch frequency."""
        return self._pitch

    @property
    def steps(self):
        """Number of steps in each cycle."""
        return self._steps

    @property
    def base(self):
        """Cycle base coefficient."""
        return self._base

    @property
    def f0(self):
        """Frequency of the first note of the scale.

        This corresponds to the frequency of the first step in the first cycle
        of the scale, for example 'A0' with classic scales in 'A' using
        relative notation or 'A4' using absolute notation."""
        return self._f0

    @property
    def abs_cycles(self):
        """Get or set absolute (True) or relative (False) cycle notation."""
        return self._abs_cycles

    @abs_cycles.setter
    def abs_cycles(self, value):
        self._abs_cycles = value

    def get_note_step(self, note):
        """Get the step number corresponding to a given note name."""
        return self._note_step[note]

    def get_freq(self, step, cycle):
        """Get the frequency for a given step and cycle numbers.

        This is an abstract function which needs to be implemented by each
        concrete scale class."""
        raise NotImplementedError

    def get_note(self, note, key=None):
        """Get the step and cycle values for a given note name."""
        if key is None:
            key = self.key
        if (len(note) > 1) and (note[1] in ['#', 'b']):
            n = 2
        else:
            n = 1
        step = self.get_note_step(note[:n]) - self.get_note_step(key)
        if len(note) > n:
            cycle = int(note[n:])
        elif self.abs_cycles:
            cycle = 4
        else:
            cycle = 0
        if self.abs_cycles:
            cycle -= 4
        elif step < 0:
            step += 12
        return step, cycle

    def get_note_freq(self, note):
        """Get the frequency for a given note name."""
        step, cycle = self.get_note(note)
        return self.get_freq(step, cycle)

    def __getitem__(self, note):
        """Directly get the frequency associated with a given note name."""
        return self.get_note_freq(note)


class LogScale(Scale):
    """Logarithmic scale, also known as equi-tempered scale.

    All steps are equal in order to evenly distribute the dissonance throughout
    the scale.  This is what is typically used on all modern instruments,
    keyboards in particular.
    """
    def __init__(self, *args, **kw):
        super(LogScale, self).__init__(*args, **kw)
        step, cycle = self.get_note(self.key, 'A')
        self._f0 = self.pitch * math.pow(self.base, float(step) / self.steps)

    def get_freq(self, step, cycle):
        step += (cycle * 12)
        return self.f0 * math.pow(self.base, float(step) / self.steps)


class HarmonicScale(Scale):
    """Harmonic scale, also known as diatonic scale.

    This is based on natural harmonics of the fundamental frequency, which
    results in some consonent intervals while others will sound more dissonent.
    It can be found with diatonic instruments, which are usually designed to
    work in only one key and are typically used in folk music.
    """
    def __init__(self, k=3, *args, **kw):
        super(HarmonicScale, self).__init__(*args, **kw)
        # Alternative shorter implementation which creates rounding errors:
        # f = math.pow(k, n)
        # f = math.pow(self.base, math.log(f, self.base) % 1)
        coefs = list()
        m = self._steps / 2
        n = self._steps - m
        for i in range(m):
            r = math.pow(k, i)
            while r > self.base:
                r /= self.base
            coefs.append(r)
        for i in range(1, n + 1):
            r = math.pow(k, -i)
            while r < 1.0:
                r *= self.base
            coefs.append(r)
        self._coefs = tuple(sorted(coefs))
        step, cycle = self.get_note(self.key, 'A')
        self._f0 = self.pitch * self._coefs[step]

    def get_freq(self, step, cycle):
        f1 = self.f0 * self._coefs[step]
        return f1 * math.pow(self.base, cycle)
