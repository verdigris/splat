import math
from data import Fragment
import _geomusic

class Generator(object):
    def __init__(self, frag, levels=None):
        self.frag = frag
        if levels is None:
            levels = tuple([1.0 for x in range(self.frag.channels)])
        self.levels = levels

    @property
    def levels(self):
        return self._levels

    @levels.setter
    def levels(self, values):
        self._check_levels(values)
        self._levels = values

    def sine(self, freq, start, stop, levels=None):
        if levels is None:
            levels = self._levels
        else:
            self._check_levels(levels)

        frag = Fragment(len(levels), self.frag.sample_rate, (stop - start))
        _geomusic.sine(frag, freq, levels)

        if True:
            fade = min((self.frag.sample_rate * 0.008), (len(frag) / 2))
            for i in xrange(int(fade)):
                l = i / fade
                z = ()
                for channel in frag[i]:
                    z += ((channel * l),)
                frag[i] = z
                z = ()
                for channel in frag[-i]:
                    z += ((channel * l),)
                frag[-i] = z

        self.frag.mix(frag, start)

    def _check_levels(self, levels):
        if len(levels) != self.frag.channels:
            raise Exception("Channels mismatch")
