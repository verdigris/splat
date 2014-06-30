# Splat - splat/generator.py
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
import random
import sources
import interpol
from data import Fragment
from filters import FilterChain
import _splat

class Generator(object):

    """Sound data generator

    This abstract class provides the basic interface to constitute a sound
    generator.  It uses a :py:class:`splat.data.Fragment` object to store the
    generated and mixed down sound data.  A generator typically runs a sound
    source with start and end times, and parameters such as a frequency and
    amplitude.  It allows arbitrary additional arguments to be passed on to the
    sound source for extra flexibility.

    The main purpose is to allow a *splat* to be run with different
    generators or sound sources without rewriting the code and data that define
    the contents of the piece.

    In order to be used, this class typically needs to be sub-classed with a
    concrete implementation of :py:meth:`splat.gen.Generator._run`.  Another
    possibility is to override the default :py:meth:`splat.gen.Generator.run`
    to define new behaviours.
    """

    def __init__(self, frag=None, filters=None):
        """The ``frag`` argument must be a :py:class:`splat.data.Fragment`
        instance.  If `None`, a default empty fragment will be automatically
        created (2 channels, 48kHz).

        A chain of ``filters`` can also be initialised here with a list of
        filter functions and internally create a
        :py:meth:`splat.filters.FilterChain` object.  This can be altered later
        via :py:attr:`splat.gen.Generator.filters`.
        """
        self.frag = frag
        self._filter_chain = FilterChain(filters)
        self._levels = tuple([0.0 for x in range(self.frag.channels)])

    @property
    def levels(self):
        """Sound levels as a tuple in dB."""
        return self._levels

    @levels.setter
    def levels(self, value):
        if not isinstance(value, float) and len(value) != self.frag.channels:
            raise ValueError("Invalid levels value")
        self._levels = value

    @property
    def frag(self):
        """:py:class:`splat.data.Fragment` instance with the generated
        audio data."""
        return self._frag

    @frag.setter
    def frag(self, value):
        if value is None:
            self._frag = Fragment()
        elif not isinstance(value, Fragment):
            raise TypeError("Fragment must be a splat.data.Fragment")
        else:
            self._frag = value

    @property
    def channels(self):
        """Number of channels."""
        return self.frag.channels

    @property
    def rate(self):
        """Sample rate in Hz."""
        return self.frag.rate

    @property
    def filters(self):
        """The :py:class:`splat.filters.FilterChain` being used."""
        return self._filter_chain

    @filters.setter
    def filters(self, filters):
        self._filter_chain = FilterChain(filters)

    def _run(self, frag, *args, **kw):
        """Main method, designed to be invoked by sub-classes via
        :py:meth:`splat.gen.Generator.run`
        """
        raise NotImplemented

    def run(self, start, end, *args, **kw):
        """Main public method to run the generator

        This method is the main entry point to run the generator and actually
        produce some sound data.  It is designed to be easily overridden by
        various types of generators, and will typically call
        :py:meth:`splat.gen.Generator._run` with a sound source and specific
        arguments.

        The ``start`` and ``end`` arguments are when the sound should start and
        end in seconds.  The ``levels`` keyword can be used to override the
        default values stored in the Generator.levels and passed to the
        ``_run`` method.

        When ``_run`` has been performed on the new fragment, filters are then
        run on it and it is finally mixed with the main internal fragment.

        It's also possible to pass a ``mix_levels`` keyword argument which is
        used when mixing the newly generated fragment into the main generator
        fragment, after the ``_run`` method has been called.  This is
        especially useful when using gains as signals with time relative to the
        beginning and length of the main fragment.
        """
        levels = kw.pop('levels', self._levels)
        start = float(start)
        end = float(end)
        frag = Fragment(self.channels, self.rate, (end - start))
        self._run(frag, levels, *args, **kw)
        self.filters.run(frag)
        self.frag.mix(frag, start, levels=kw.pop('mix_levels', None))


class SourceGenerator(Generator):

    """Generator using a sound source

    This is a basic class to implement generators using a sound source.  See
    :ref:`sources` for examples.  Several sub-classes are already available to
    make best use of the built-in sound sources.
    """

    def __init__(self, source, *args, **kw):
        """The ``source`` argument is meant to be a source function.  This is
        then used in the :py:meth:`splat.gen.SourceGenerator._run`
        implementation of this class.
        """
        super(SourceGenerator, self).__init__(*args, **kw)
        self._source = source

    @property
    def source(self):
        """Sound source."""
        return self._source

    def _run(self, frag, levels, freq, phase=0.0, *args, **kw):
        self.source(frag, levels, freq, phase, *args)


class SineGenerator(SourceGenerator):

    """Sine wave generator

    This is the simplest generator, based on the :py:func:`splat.sources.sine`
    source to generate pure sine waves.
    """

    def __init__(self, *args, **kw):
        super(SineGenerator, self).__init__(sources.sine, *args, **kw)


class SquareGenerator(SourceGenerator):

    """Square wave generator

    This uses the :py:func:`splat.sources.square` source to generate square
    waves.
    """
    def __init__(self, *args, **kw):
        super(SquareGenerator, self).__init__(sources.square, *args, **kw)


class TriangleGenerator(SourceGenerator):

    """Triangle wave generator

    This uses the :py:func:`splat.sources.triangle` source to generate triangle
    waves.
    """
    def __init__(self, *args, **kw):
        super(TriangleGenerator, self).__init__(sources.triangle, *args, **kw)


class OvertonesGenerator(SourceGenerator):

    """Overtones generator

    Overtones are defined by an ``overtones`` list.  For a description of
    the overtunes, see the :py:func:`splat.sources.overtones` source which is
    used by this generator.

    Note: The time to generate the signal increases with the number of
    overtones.
    """

    def __init__(self, *args, **kw):
        super(OvertonesGenerator, self).__init__(sources.overtones, *args,**kw)
        self.overtones = [(1.0, 0.0, 0.0)]

    def ot_decexp(self, k=1.0, n=24):
        """Set harmonic overtones levels following a decreasing exponential.

        Harmonics are defined here as overtones with positive integer ratios.
        For a given parameter ``k`` and a harmonic ``i`` from 1 to a total
        number of harmonics ``n``, the amplitude of each overtone is set
        following this function:

        .. math::

           o[i] = exp\\left(\\frac{1 - i}{k}\\right)

        As a result, the fundamental frequency (overtone 1.0) always has an
        amplitude of 1.0 (0 dB).

        A higher ``k`` value means the function will decrease faster causing
        less high-frequency harmonics.
        """
        self.overtones = list()
        for j in (float(i) for i in range(n)):
            l = _splat.lin2dB(math.exp(-j / k))
            self.overtones.append(((j + 1), 0.0, l))

    def _run(self, frag, levels, freq, phase=0.0, *args, **kw):
        super(OvertonesGenerator, self)._run(frag, levels, freq, phase,
                                             self.overtones, *args, **kw)


class Particle(object):

    """Sound particle

    A sound particle is a basic element to generate complex sounds via a
    :py:class:`splat.gen.ParticleGenerator` object.  It only contains the
    parameters to be used in combination with a sound source to create a
    usually short piece of sound.  Sound particles are usually grouped in a
    :py:class:`splat.gen.ParticlePool` object.
    """

    def __init__(self, start, end, f_log):
        """The ``start`` and ``end`` arguments set the time limits of the
        particle in seconds.  The ``f_log`` argument is a frequency on a
        logarithmic scale, which is then typically used by the sound source to
        generate the actual piece of sound."""
        self._start = start
        self._end = end
        self.f_log = f_log

    def __repr__(self):
        return u"[{0}, {1}] {2}".format(self.start, self.end, self.freq)

    @property
    def start(self):
        """Start time in seconds."""
        return self._start

    @property
    def end(self):
        """End time in seconds."""
        return self._end

    @property
    def freq(self):
        """Frequency as a linear value in Hz."""
        return _splat.dB2lin(self.f_log)

    @property
    def length(self):
        """Length in seconds."""
        return self.end - self.start


class ParticlePool(object):

    """Pool of sound particles

    The idea behind sound particles is to use statistics over a large quantity
    of discrete particles to create a continuous sound.  To achieve this, the
    ParticlePool can be used to create a population of particles with an
    envelope of varying density with some randomness in each particle
    parameters.
    """

    def __init__(self, min_f_log, max_f_log, min_len, max_len, envelope,
                 n_slices, density):
        """The ``min_f_log`` and ``max_f_log`` arguments give the boundaries
        for the particle frequencies using logarithmic scale.  The ``min_len``
        and ``max_len`` arguments are for the minimum and maximum length of the
        particle in seconds.  The ``envelope`` is a
        :py:class:`splat.interpol.Spline` object which is then used in slices
        to create more or less particles around each point in time.  The number
        of slices given by ``n_slices`` affects the granularity of the
        envelope.  The ``density`` affects the average number of particle per
        unit of time."""
        self._pts = []
        for y0 in range(n_slices):
            y0 = float(1 + y0) / n_slices
            seg = envelope.slices(y0)

            for x0, x1 in seg:
                if (x1 - x0) < min_len:
                    continue

                slice_max_len = min((x1 - x0), max_len)
                n_events = int((x1 - x0) * density / n_slices)
                for i in range(n_events):
                    length = random.uniform(min_len, slice_max_len)
                    start = random.uniform(x0, (x1 - length))
                    end = start + length
                    f_log = random.uniform(min_f_log, max_f_log)
                    self._pts.append(Particle(start, end, f_log))

        self._start = envelope.start
        self._end = envelope.end


    @property
    def start(self):
        """Start time in seconds."""
        return self._start

    @property
    def end(self):
        """End time in seconds."""
        return self._end

    def count(self, t_start=None, t_end=None):
        """Number of particles left in the pool, either in total or in between
        the ``t_start`` and ``t_end`` times in seconds if they are both
        specified."""
        if (t_start is None) and (t_end is None):
            return len(self._pts)
        n = 0
        for p in self._pts:
            if (p.start > t_start) and (p.end < t_end):
                n += 1
        return n

    def iterate(self, t_start=None, t_end=None, share=1.0):
        """Python generator to iterate through all the particles in the pool
        and remove them at the same time, or only in between the ``t_start``
        and ``t_end`` times in second if they are both specified.  The
        ``share`` argument can be used to only consume a share of the available
        particles.  For example, a value ``share=0.5`` will only use half of
        the available particles."""
        if t_start is None:
            t_start = self.start
        if t_end is None:
            t_end = self.end
        n = int(self.count(t_start, t_end) * share)
        new_pts = []
        for i, p in enumerate(self._pts):
            if (n > 0) and (p.start > t_start) and (p.end < t_end):
                n -= 1
                yield(p)
            else:
                new_pts.append(p)
        self._pts = new_pts


class ParticleGenerator(Generator):

    """Generator with particles

    This generator uses a :py:class:`splat.gen.ParticlePool` with an internal
    :py:class:`splat.gen.Generator` and some parameters to generate a sound.
    """

    def __init__(self, subgen, min_f=20.0, max_f=20000.0, *args, **kw):
        """The ``subgen`` object is a :py:class:`splat.gen.Generator` which is
        called once for each :py:class:`splat.gen.Particle` object used to
        produce the final sound.  Due to the nature of the Particle objects, it
        should take a frequency argument."""
        super(ParticleGenerator, self).__init__(subgen.frag, *args, **kw)
        self._subgen = subgen
        self._start = None
        self._end = None
        self._z = None
        self._eq = None
        self._q = None
        self._min_f_log = _splat.lin2dB(min_f)
        self._max_f_log = _splat.lin2dB(max_f)
        self._gain_fuzz = 0.0
        self._relative_gain_fuzz = 0.0
        self._pool = None
        self.progress_step = 10
        self.do_show_progress = True

    @property
    def subgen(self):
        """Sub-generator object."""
        return self._subgen

    @property
    def start(self):
        """Start time in seconds."""
        return self._start

    @property
    def end(self):
        """End time in seconds."""
        return self._end

    # ToDo: use dB instead of linear (0..1) levels?
    def set_z(self, z_pts, *args, **kw):
        """Create the :py:class:`splat.interpol.Spline` object used for the
        :py:class:`splat.gen.ParticlePool` amplitude envelope with the given
        ``z_pts`` points."""
        self._z = interpol.Spline(z_pts, *args, **kw)

    @property
    def z(self):
        """Get the amplitude envelope Spline object."""
        return self._z

    def set_eq(self, eq_pts, *args, **kw):
        """Set the equalization :py:class:`splat.interpol.Spline` object used
        to alter the amplitude of each particle in function of its frequency to
        create an equalization.  The minimum and maximum frequencies used when
        creating the :py:class:`splat.gen.ParticlePool` is derived from the
        boundaries of the equalization Spline."""
        self._eq = interpol.Spline(eq_pts, *args, **kw)
        self._min_f_log = _splat.lin2dB(self._eq.start)
        self._max_f_log = _splat.lin2dB(self._eq.end)

    @property
    def eq(self):
        """Get the equalization Spline object."""
        return self._eq

    def set_q(self, q_pts, *args, **kw):
        """Set the distribution parameter :py:class:`splat.interpol.Spline`
        object which is used to alter the random function used to pick the
        frequency of each :py:class:`splat.gen.Particle` using the given
        ``q_pts`` points."""
        self._q = interpol.Spline(q_pts, *args, **kw)

    @property
    def q(self):
        """Get the distribution parameter Spline object."""
        return self._q

    @property
    def gain_fuzz(self):
        """Get a 2-tuple with the common and relative gain fuzz values in dB."""
        return (self._gain_fuzz, self._relative_gain_fuzz)

    @gain_fuzz.setter
    def gain_fuzz(self, value):
        """Set the level of gain fuzz which is the amount of randomness added
        to the gain of each :py:class:`splat.gen.Particle` object.  The value
        passed is a 2-tuple with the common gain fuzz applied to all channels
        and the relative gain fuzz used to create a difference in between each
        channel.  These values are in dB."""
        self._gain_fuzz, self._relative_gain_fuzz = value

    def make_pool(self, min_len=0.05, max_len=0.1, n_slices=20, density=100):
        """Build the pool of particles using all the parameters currently
        defined.  This discards any existing pool, and can also be used to
        create a new pool with the same parameters when the previous one has
        been exhausted."""
        self._pool = ParticlePool(
            self._min_f_log, self._max_f_log, min_len, max_len, self._z,
            n_slices, density)

    @property
    def pool(self):
        """Get the :py:class:`splat.gen.ParticlePool` object."""
        return self._pool

    def run(self, start, end, freq, share=1.0, levels=None, *args, **kw):
        """Run the generator using the standard interface, with the extra
        ``share`` argument which is passed to the internal
        :py:meth:`splat.gen.ParticlePool.iterate` object."""
        if self.start is None:
            self._start = start
            self._end = end
        else:
            self._start = min(self.start, start)
            self._end = max(self.end, end)

        freq = _splat.Signal(self.frag, freq, (self.end - self.start))
        n_events = self._pool.count()
        step = n_events / self.progress_step
        progress = 0
        for i, p in enumerate(self._pool.iterate(start, end, share)):
            if self.do_show_progress and (i % step) == 0:
                self.show_progress(progress)
                progress += self.progress_step

            if self._eq:
                g = self._eq.value(p.freq)
            else:
                g = 0.0

            if self._gain_fuzz:
                g += random.uniform(-self._gain_fuzz, self._gain_fuzz)

            if self._relative_gain_fuzz:
                levels = tuple(g + random.uniform(-self._relative_gain_fuzz,
                                                   self._relative_gain_fuzz)
                               for g in range(self.frag.channels))
            else:
                levels = tuple(g for i in range(self.frag.channels))

            if self._q:
                p_freq = self.curve(freq[self.frag.s2n(p.start)][0], p.freq,
                                    self._q.value(p.start))
            else:
                p_freq = p.freq

            self.subgen.run(p.start, p.end, p_freq, levels=levels, *args, **kw)

    def curve(self, freq, p_freq, q):
        """This method implements the stastistical distribution used to pick
        the frequency for each :py:class:`splat.gen.Particle` object.  It
        returns a frequency ``freq`` in linear scale altered around the target
        frequency ``p_freq`` with the parameter ``q``. The default
        implementation produces something a little bit like a normal
        distribution but is not a standard function."""
        f_log = _splat.lin2dB(freq)
        s = _splat.lin2dB(p_freq)
        f_curve = s - 0.5 * ((2 * s) - (2 * f_log)) * math.exp(
            -((s - f_log) * (s - f_log)) / _splat.dB2lin(q))
        return _splat.dB2lin(f_curve)

    def show_progress(self, progress):
        """Print the ``progress`` in percentage of the
        :py:meth:`splat.gen.ParticleGenerator.run` method call."""
        print("Progress: {0}%".format(progress))
