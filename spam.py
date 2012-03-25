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
