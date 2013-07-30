# Splat - splat/generator.py
#
# Copyright (C) 2012, 2013 Guillaume Tucker <guillaume@mangoz.org>
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
import sources
from data import Fragment
from filters import FilterChain
import _splat

class Generator(object):

    """Sound data generator.

    This abstract class provides the basic interface to constitute a sound
    generator.  It creates a :py:class:`splat.data.Fragment` object to store
    the generated and mixed down sound data.  A generator typically runs a
    sound source with a given base frequency and start and end times.  It
    allows arbitrary extra arguments to be passed to the sound source for extra
    flexibility.

    The main purpose is to allow a Splat piece to be run with different
    generators or sound sources without rewriting the code and data that define
    the contents of the piece.

    In order to be used, this class typically needs to be sub-classed to
    implement :py:meth:`splat.gen.Generator._run` with a concrete sound source.
    """

    def __init__(self, frag, filters=None):
        """The ``frag`` argument must be a :py:class:`splat.data.Fragment`
        instance.

        A chain of ``filters`` can also be initialised here with a list of
        filter functions and internally create a
        :py:meth:`splat.filters.FilterChain` object.  This can be altered later
        via :py:attr:`splat.gen.Generator.filters`.
        """
        self._frag = frag
        self._filter_chain = FilterChain(filters)
        self._levels = tuple([0.0 for x in range(self.frag.channels)])
        self._time_stretch = 1.0

    @property
    def levels(self):
        """Sound levels as a tuple in dB."""
        return self._levels

    @levels.setter
    def levels(self, values):
        if len(values) != self.frag.channels:
            raise Exception("Channels mismatch")
        self._levels = values

    @property
    def frag(self):
        """:py:class:`splat.data.Fragment` instance with the generated
        sounds."""
        return self._frag

    @property
    def channels(self):
        """Number of channels."""
        return self._frag.channels

    @property
    def sample_rate(self):
        """Sample rate in Hz."""
        return self._frag.sample_rate

    @property
    def filters(self):
        """The :py:class:`splat.filters.FilterChain` being used."""
        return self._filter_chain

    @filters.setter
    def filters(self, filters):
        self._filter_chain = FilterChain(filters)

    @property
    def time_stretch(self):
        """Time stretch factor

        All ``start`` and ``end`` times are multiplied by this value when
        calling :py:meth:`splat.gen.Generator.run`.
        """
        return self._time_stretch

    @time_stretch.setter
    def time_stretch(self, value):
        self._time_stretch = value

    def _run(self, frag, *args, **kw):
        """Main method, designed to be invoked by sub-classes via
        :py:meth:`splat.gen.Generator.run`
        """
        raise NotImplemented

    def run(self, freq, start, end, levels=None, *args, **kw):
        """Main public method to run the generator

        This method is the main entry point to run the generator and actually
        produce some sound data.  It is designed to be easily overridden by
        various types of generators, and will typically call
        :py:meth:`splat.gen.Generator._run` with a sound source and specific
        arguments.
        """
        if levels is None:
            levels = self._levels
        start *= self.time_stretch
        end *= self.time_stretch
        frag = Fragment(self.channels, self.sample_rate, (end - start))
        self._run(frag, freq, levels, *args, **kw)
        self.filters.run(frag)
        self.frag.mix(frag, start)


class SourceGenerator(Generator):

    def __init__(self, source, *args, **kw):
        super(SourceGenerator, self).__init__(*args, **kw)
        self._source = source

    @property
    def source(self):
        """Sound source."""
        return self._source

    def _run(self, frag, *args, **kw):
        self.source(frag, *args, **kw)


class SineGenerator(SourceGenerator):

    """Sine wave generator.

    This is the simplest generator, based on the :py:func:`splat.sources.sine`
    source to generate pure sine waves.
    """

    def __init__(self, *args, **kw):
        super(SineGenerator, self).__init__(sources.sine, *args, **kw)


class SquareGenerator(SourceGenerator):

    def __init__(self, *args, **kw):
        super(SquareGenerator, self).__init__(sources.square, *args, **kw)


class TriangleGenerator(SourceGenerator):

    def __init__(self, *args, **kw):
        super(TriangleGenerator, self).__init__(sources.triangle, *args, **kw)


class OvertonesGenerator(SourceGenerator):

    """Overtones generator.

    Overtones are defined by an ``overtones`` dictionary.  For a description of
    the overtunes, see the :py:func:`splat.sources.overtones` source which is
    used by this generator.

    Note: The time to generate the signal increases with the number of
    overtones.
    """

    def __init__(self, *args, **kw):
        super(OvertonesGenerator, self).__init__(sources.overtones,
                                                 *args, **kw)
        self.overtones = { 1.0: 0.0 }

    def ot_decexp(self, k=1.0, n=24):
        """Set harmonic overtones levels following a decreasing exponential.

        For a given ratio ``k`` and a harmonic ``i`` from 1 to a total number
        of harmonics ``n``, the amplitude of each overtone is set following
        this function:

        .. math::

           o[i] = exp\\left(\\frac{1 - i}{k}\\right)

        As a result, the fundamental frequency (overtone 1.0) always has an
        amplitude of 1.0 (0 dB).

        A higher ``k`` value means the function will decrease faster causing
        less high-frequency harmonics.
        """
        self.overtones = dict()
        for j in (float(i) for i in range(n)):
            l = _splat.lin2dB(math.exp(-j / k))
            self.overtones[j + 1] = l

    def _run(self, frag, freq, levels, *args, **kw):
        super(OvertonesGenerator, self)._run(frag, freq, levels,
                                             self.overtones, *args, **kw)
