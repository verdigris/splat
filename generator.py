import math
from spam import sine

class Generator(object):
    def __init__(self, frag):
        self.frag = frag
        self._levels = None

    def set_levels(self, levels):
        self._check_levels(levels)
        self._levels = levels

    def sine(self, freq, start, stop, levels=None):
        if levels is None:
            levels = self._levels
        else:
            self._check_levels(levels)
        frag = sine(self.frag.sample_rate, freq, (stop - start), levels)
        if True:
            fade = min((self.frag.sample_rate * 0.01), (frag.length / 2))
            for i in range(int(fade)):
                l = i / fade
                for channel in frag.data:
                    channel[i] *= l
                    channel[-i] *= l
        self.frag.mix(frag, start)

    def _check_levels(self, levels):
        if len(levels) != self.frag.channels:
            raise Exception("Channels mismatch")
