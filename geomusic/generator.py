import math
import collections
import sources
from data import Fragment
import _geomusic

class Generator(object):
    """Generator to manage sound sources

    This class needs to be sub-classed to implement
    :py:meth:`geomusic.Generator.run` with a concreate sound source.  It
    creates a :py:class:`geomusic.Fragment` instance to store the generated
    and mixed sounds.
    """
    def __init__(self, frag, levels=None):
        """And this is the constructor

        A :py:class:`geomusic.Fragment` instance needs to be passed by the
        ``frag`` argument.
        """
        self._frag = frag
        if levels is None:
            levels = tuple([1.0 for x in range(self.frag.channels)])
        self.levels = levels
        self.chain = None
        self.time_stretch = 1.0

    @property
    def levels(self):
        """Sound levels as a tuple in dB"""
        return self._levels

    @levels.setter
    def levels(self, values):
        self._check_levels(values)
        self._levels = values

    @property
    def frag(self):
        """:py:class:`geomusic.Fragment` instance with the generated sounds"""
        return self._frag

    @property
    def channels(self):
        """Number of channels"""
        return self._frag.channels

    @property
    def sample_rate(self):
        """Sample rate in Hz"""
        return self._frag.sample_rate

    def _run(self, source, start, stop, *args, **kw):
        """Main method, designed to be invoked by sub-classes from
        :py:meth:`geomusic.Generator.run`

        The ``source`` argument is a sound source function, to which the
        ``*args`` and ``**kw`` arguments are passed on.  The sound is
        supposed to start and stop at the times specified by the ``start``
        and ``stop`` arguments.
        """
        start *= self.time_stretch
        stop *= self.time_stretch
        frag = Fragment(self.channels, self.sample_rate, (stop - start))
        source(frag, *args, **kw)
        if self.chain:
            self.chain.run(frag)
        self.frag.mix(frag, start)

    def run(self, start, stop, *args, **kw):
        """Main method to be implemented by sub-classes

        This method is the main entry point to run the generator and actually
        produce sound.  It will typically call
        :py:meth:`geomusic.Generator._run` with a sound source and specific
        arguments.
        """
        raise NotImplementedError

    def _check_levels(self, levels):
        if len(levels) != self.frag.channels:
            raise Exception("Channels mismatch")


class FilterChain(collections.Sequence):
    def __init__(self, filters=list(), *args, **kw):
        super(FilterChain, self).__init__(*args, **kw)
        self._filters = list()
        for f in filters:
            if isinstance(f, tuple):
                self.append(*f)
            else:
                self.append(f)

    def __getitem__(self, i):
        return self._filters[i]

    def __len__(self):
        return len(self._filters)

    def append(self, filter_func, args=()):
        if not isinstance(args, tuple):
            raise Exception("Invalid filter arguments, must be a tuple")
        self._filters.append((filter_func, args))

    def run(self, frag):
        for f, args in self:
            f(frag, *args)


class SineGenerator(Generator):
    """Sine wave generator.

    This is the simplest generator, based on :py:func:`geomusic.sources.sine`
    to generate pure sine waves.
    """

    def run(self, freq, start, stop, levels=None):
        if levels is None:
            levels = self._levels
        self._run(sources.sine, start, stop, freq, levels)


class OvertonesGenerator(Generator):
    """Overtones generator.

    Overtones are defined by a ``overtones`` dictionary.  This uses the
    :py:func:`geomusic.sources.overtones` source.

    Note: The time to generate the signal increases with the number of
    overtones ``n``.
    """
    def __init__(self, *args, **kw):
        super(OvertonesGenerator, self).__init__(*args, **kw)
        self.overtones = { 1.0: tuple(0.0 for i in range(self.frag.channels)) }

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
            l = _geomusic.lin2dB(math.exp(-j / k))
            self.overtones[j + 1] = tuple(l for i in range(self.frag.channels))

    def run(self, freq, start, stop, levels=None):
        if levels is None:
            levels = self._levels
        self._run(sources.overtones, start, stop, freq, levels, self.overtones)
