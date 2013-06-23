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

    """Sound data generator.

    This abstract class provides the basic interface to constitute a sound
    generator.  It uses a :py:class:`splat.data.Fragment` object to store the
    generated and mixed down sound data.  A generator typically runs a sound
    source with start and end times, and parameters such as a frequency and
    amplitudes.  It allows arbitrary extra arguments to be passed on to the
    sound source for extra flexibility.

    The main purpose is to allow a Splat piece to be run with different
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
        created.

        A chain of ``filters`` can also be initialised here with a list of
        filter functions and internally create a
        :py:meth:`splat.filters.FilterChain` object.  This can be altered later
        via :py:attr:`splat.gen.Generator.filters`.
        """
        if frag is None:
            frag = Fragment()
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

    def run(self, start, end, *args, **kw):
        """Main public method to run the generator

        This method is the main entry point to run the generator and actually
        produce some sound data.  It is designed to be easily overridden by
        various types of generators, and will typically call
        :py:meth:`splat.gen.Generator._run` with a sound source and specific
        arguments.

        The ``start`` and ``end`` arguments are when the sound should start and
        end in seconds.  The rest of the arguments are passed to the ``_run``
        method, with the ``levels`` keywords always in second place after the
        internal fragment.

        When ``_run`` has been run on the new fragment, filters are run on it
        and it is then mixed with the main internal fragment.
        """
        levels = kw.pop('levels', self._levels)
        start *= self.time_stretch
        end *= self.time_stretch
        frag = Fragment(self.channels, self.sample_rate, (end - start))
        self._run(frag, levels, *args, **kw)
        self.filters.run(frag)
        self.frag.mix(frag, start)


class SourceGenerator(Generator):

    """Generator using a sound source

    This is a basic class to implement generators using a sound source.  See
    :ref:`sources` for examples.  Several sub-classes are already available to
    make best use of the built-in sound sources.
    """

    def __init__(self, source, *args, **kw):
        """The `source` argument is meant to be a source function.  This is
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
        self.source(frag, levels, freq, phase, *args, **kw)


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
        self.overtones = list()
        for j in (float(i) for i in range(n)):
            l = _splat.lin2dB(math.exp(-j / k))
            self.overtones.append(((j + 1), 0.0, l))

    def _run(self, frag, levels, freq, phase=0.0, *args, **kw):
        super(OvertonesGenerator, self)._run(frag, levels, freq, phase,
                                             self.overtones, *args, **kw)


class Particle(object):

    def __init__(self, start, end, f_log):
        self.start = start
        self.end = end
        self.f_log = f_log

    def __repr__(self):
        return u"[{0}, {1}] {2}".format(self.start, self.end, self.freq)

    @property
    def freq(self):
        return _splat.dB2lin(self.f_log)


class ParticlePool(object):

    def __init__(self, min_f_log, max_f_log, min_len, max_len, envelope,
                 n_slices, density):
        self._pts = []
        for y0 in range(n_slices):
            y0 = float(1 + y0) / n_slices
            seg = envelope.slices(y0)

            for x0, x1 in seg:
                n_events = int((x1 - x0) * density / n_slices)

                for i in range(n_events):
                    start = random.uniform(x0, x1)
                    end = random.uniform((start + min_len), (start + max_len))
                    f_log = random.uniform(min_f_log, max_f_log)
                    self._pts.append(Particle(start, end, f_log))

        self._start = envelope.start
        self._end = envelope.end


    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    def count(self, t_start=None, t_end=None):
        if (t_start is None) and (t_end is None):
            return len(self._pts)
        n = 0
        for p in self._pts:
            if (p.start > t_start) and (p.end < t_end):
                n += 1
        return n

    def iterate(self, t_start=None, t_end=None, share=1.0):
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

    def __init__(self, subgen, *args, **kw):
        super(ParticleGenerator, self).__init__(subgen.frag, *args, **kw)
        self._subgen = subgen
        self._start = None
        self._end = None
        self._z = None
        self._eq = None
        self._q = None
        self._min_f_log = _splat.lin2dB(20)
        self._max_f_log = _splat.lin2dB(20000)
        self._gain_fuzz = 0.0
        self._stereo_gain_fuzz = 0.0
        self._pool = None
        self.progress_step = 10
        self.do_show_progress = True

    @property
    def subgen(self):
        return self._subgen

    @property
    def start(self):
        return self._start

    @property
    def end(self):
        return self._end

    # ToDo: use dB instead of linear (0..1) levels?
    def set_z(self, z_pts, *args, **kw):
        self._z = interpol.Spline(z_pts, *args, **kw)

    @property
    def z(self):
        return self._z

    def set_eq(self, eq_pts, *args, **kw):
        self._eq = interpol.Spline(eq_pts, *args, **kw)
        self._min_f_log = _splat.lin2dB(self._eq.start)
        self._max_f_log = _splat.lin2dB(self._eq.end)

    @property
    def eq(self):
        return self._eq

    def set_q(self, q_pts, *args, **kw):
        self._q = interpol.Spline(q_pts, *args, **kw)

    @property
    def q(self):
        return self._q

    def set_gain_fuzz(self, gain_fuzz, stereo_gain_fuzz):
        self._gain_fuzz = gain_fuzz
        self._stereo_gain_fuzz = stereo_gain_fuzz

    @property
    def gain_fuzz(self):
        return (self._gain_fuzz, self._stereo_gain_fuzz)

    def make_pool(self, min_len=0.05, max_len=0.1, n_slices=20, density=100):
        self._pool = ParticlePool(
            self._min_f_log, self._max_f_log, min_len, max_len, self._z,
            n_slices, density)

    @property
    def pool(self):
        return self._pool

    def run(self, start, end, freq, share=1.0, levels=(0.0, 0.0), *args, **kw):
        if self.start is None:
            self._start = start
            self._end = end
        else:
            self._start = min(self.start, start)
            self._end = max(self.end, end)

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

            if self._stereo_gain_fuzz:
                levels = tuple(g + random.uniform(-self._stereo_gain_fuzz,
                                                   self._stereo_gain_fuzz)
                               for g in range(self.frag.channels))
            else:
                levels = tuple(g for i in range(self.frag.channels))

            if self._q:
                p_freq = self.curve(freq, p.freq, self._q.value(p.start))
            else:
                p_freq = p.freq

            self.subgen.run(p.start, p.end, p_freq, levels=levels, *args, **kw)

    def curve(self, freq, p_freq, q):
        f_log = _splat.lin2dB(freq)
        s = _splat.lin2dB(p_freq)
        f_curve = s - 0.5 * ((2 * s) - (2 * f_log)) * math.exp(
            -((s - f_log) * (s - f_log)) / _splat.dB2lin(q))
        return _splat.dB2lin(f_curve)

    def show_progress(self, progress):
        print("Progress: {0}%".format(progress))
