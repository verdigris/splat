# Splat - splat/seq.py
#
# Copyright (C) 2014, 2015 Guillaume Tucker <guillaume@mangoz.org>
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

import os
import random
import _splat
from .data import Fragment

class SampleSet(object):
    """Set of samples

    A ``SampleSet`` object contains a set of fragments that can be interchanged
    and picked randomly.  The ``entries`` are a list of 3-tuples:
    (``sample``, ``weight``, ``gain``):

      ``sample``
        Fragment object or string with the name of a file to open
      ``weight``
        Weight when picking randomly a entry in the set
      ``gain``
         Linear gain to be applied to the sample using
         :py:meth:`splat.data.Fragment.amp`
    """

    def __init__(self, entries):
        self._samples = []
        self._dict = dict()
        for sample, weight, gain in entries:
            frag = Fragment.open(sample) if isinstance(sample, str) else sample
            frag.amp(gain)
            if frag.name is None:
                frag.name = str(frag)
            self._dict[frag.name] = frag
            for i in range(weight):
                self._samples.append(frag)

    def pick(self):
        """Randomly pick a sample fragment and return it."""
        return random.choice(self._samples)

    @property
    def names(self):
        """List with the sample names."""
        return self._dict.keys()

    def get(self, name):
        """Get the sample fragment that matches the given ``name``."""
        return self._dict[name]


class Pattern(object):
    """Base abstract pattern

    A ``Pattern`` object can be used by a sequencer object to generate a series
    of sounds.  The ``beats`` gives the number of times the
    :py:meth:`splat.seq.Pattern.play` method needs to be called to complete the
    pattern.
    """

    def __init__(self, beats=4):
        self._beats = beats

    @property
    def beats(self):
        """Number of beats in the pattern."""
        return self._beats

    def play(self, frag, bar, beat, t, T):
        """Abstract method, called to generate each beat in the pattern.

        The ``frag`` argument is the Fragment object into which the sound needs
        to be generated (or mixed).  The ``bar`` and ``beat`` numbers can be
        used to alter the pattern's behaviour depending on the context.

        The start time for the current beat is given by ``t``, which can be
        passed to :py:meth:`splat.seq.Pattern.mix` as the offset value.

        The duration of a beat is given by ``T``.  This is useful when patterns
        need to generate sounds that last a fraction of a beat.

        This method must be implemented in derived classes.
        """
        raise NotImplementedError

    def mix(self, master, frag, t, levels=None, skip=0.0):
        """Mix a pattern fragment into a master fragment.

        This is the preferred way to mix a pattern fragment into a master
        fragment as it makes it possible for derived classes such as
        :py:class:`splat.seq.FuzzyPattern` to alter the behaviour.
        """
        master.mix(frag, offset=t, skip=skip, levels=levels)

    def pick_mix(self, master, samples, *args, **kw):
        """Pick a sample in a set and mix it into a fragment.

        Pick a sample in the :py:class:`splat.seq.SampleSet` object ``samples``
        and mix it into the ``master`` fragment by calling
        :py:meth:`splat.seq.Pattern.mix`.
        """
        self.mix(master, samples.pick(), *args, **kw)


class FuzzyPattern(Pattern):
    """Pattern with random errors

    A FuzzyPattern can be used instead of a Pattern as a bass class to add some
    random gain and timing errors to the sample.  This is done via the
    :py:meth:`splat.seq.FuzzyPattern.mix` method.

    The ``time_error`` and ``gain_error`` parameters respectively set the
    maximum amount of error for the start time in milliseconds and in dB for
    the gain when each fragment is mixed by the pattern's
    :py:meth:`splat.seq.Pattern.play` implementation.
    """

    def __init__(self, time_error=0.02, gain_error=3.0, *args, **kw):
        super(FuzzyPattern, self).__init__(*args, **kw)
        self._te = (-time_error, time_error)
        self._ge = (-gain_error, gain_error)

    def mix(self, frag, sample, t, g=0.0):
        """Alternative mixing method implementation with random errors."""
        t += random.uniform(*self._te)
        l = tuple(g + random.uniform(*self._ge) for i in range(frag.channels))
        l = tuple(_splat.dB2lin(c) for c in l)
        frag.mix(sample, t, 0.0, l)


class Silence(Pattern):
    """Silence pattern

    This simple pattern is a convenient way of not generating anything for the
    duration of a bar.
    """

    def play(self, frag, bar, beat, t, T):
        pass


class Sequencer(object):
    """Sequencer abstract and simplistic base class

    The ``tempo`` argument is the number of beats per minute.
    """

    def __init__(self, tempo):
        self.tempo = tempo

    @property
    def tempo(self):
        """Tempo in beats per minute."""
        return self._tempo

    @tempo.setter
    def tempo(self, value):
        self._tempo = value
        self._period = 60.0 / value

    @property
    def period(self):
        """Period in seconds, duration of a beat."""
        return self._period

    @period.setter
    def period(self, value):
        self._period = value
        self._tempo = 60.0 / value

    def run(self, frag, *args, **kw):
        """Main method to be implemented by concrete classes.

        ``frag`` is the fragment object into which the generated patterns are
        mixed.  The other arguments are left to the be defined by each
        implementation.
        """
        raise NotImplementedError


class PatternSequencer(Sequencer):
    """Pattern sequencer

    A PatternSequencer object takes a list of :py:class:`splat.seq.Pattern`
    objects and iterates over them, calling their
    :py:meth:`splat.seq.Pattern.play` method each time.
    """

    def run(self, frag, patterns):
        """Run the sequencer over the given list of ``patterns``.

        Each entry in the ``patterns`` list (or any iterable sequence)
        corresponds to a bar and contains either a single or a tuple of
        :py:class:`splat.seq.Pattern` objects.  The sequencer iterates over the
        whole list and calls the :py:meth:`splat.seq.Pattern.play` method for
        each beat in each bar on each pattern object found.  Each pattern may
        have a different number of beats, but the tempo is controlled by the
        sequencer and not by the pattern.
        """
        pats = list((p,) if isinstance(p, Pattern) else p for p in patterns)
        total_beats = sum(p[0].beats for p in pats)
        frag.grow(total_beats * self.period)
        n_beat = 0
        for bar, group in enumerate(pats):
            beats = group[0].beats
            for pattern in group:
                bar_beat = n_beat
                for beat in range(beats):
                    t = bar_beat * self.period
                    pattern.play(frag, bar, beat, t, self.period)
                    bar_beat += 1
            n_beat += beats
