import math
from data import Fragment
try:
    import _geomusic
    use_c = True
except ImportError:
    use_c = False

def _c_sine(sample_rate, freq, duration, levels):
    frag = Fragment(len(levels), sample_rate)
    length = frag.resize(duration)
    _geomusic.sine(frag.data, freq, sample_rate, length, levels)
    return frag

def _py_sine(sample_rate, freq, duration, levels):
    k = 2 * math.pi * freq
    frag = Fragment(len(levels), sample_rate)
    n = frag.resize(duration)
    for i in range(n):
        s = math.sin(k * i / sample_rate)
        for c, l in enumerate(levels):
            frag.data[c][i] = s * l
    return frag

if use_c:
    sine = _c_sine
else:
    sine = _py_sine


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
