import math
import collections
import sources
from data import Fragment
import _geomusic

class Generator(object):
    def __init__(self, frag, levels=None):
        self._frag = frag
        if levels is None:
            levels = tuple([1.0 for x in range(self.frag.channels)])
        self.levels = levels
        self.chain = None
        self.time_stretch = 1.0

    @property
    def levels(self):
        return self._levels

    @levels.setter
    def levels(self, values):
        self._check_levels(values)
        self._levels = values

    @property
    def frag(self):
        return self._frag

    @property
    def channels(self):
        return self._frag.channels

    @property
    def sample_rate(self):
        return self._frag.sample_rate

    def run(self, source, start, stop, *args, **kw):
        start *= self.time_stretch
        stop *= self.time_stretch
        frag = Fragment(self.channels, self.sample_rate, (stop - start))
        source(frag, *args, **kw)
        if self.chain:
            self.chain.run(frag)
        self.frag.mix(frag, start)

    def sine(self, freq, start, stop, levels=None):
        if levels is None:
            levels = self._levels
        self.run(sources.sine, start, stop, freq, levels)

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
