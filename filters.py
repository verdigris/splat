def linear_fade(frag, duration=0.01):
    fade = min((frag.sample_rate * duration), (len(frag) / 2))
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

def normalize(frag, level=1.0):
    class Stats(object):
        def __init__(self):
            self.min = 0.0
            self.max = 0.0
            self._total = 0.0
            self._n = 0

        @property
        def average(self):
            return self._total / self._n

        @property
        def peak(self):
            avg = self.average
            peak_pos = self.max + avg
            peak_neg = self.min + avg
            return max(abs(peak_pos), abs(peak_neg))

        def process(self, s):
            self.min = min(s, self.min)
            self.max = max(s, self.max)
            self._total += s
            self._n += 1

    stats = [Stats() for x in xrange(frag.channels)]
    for s in frag:
        for c, c_stats in zip(s, stats):
            c_stats.process(c)

    peak = 0.0
    for c in stats:
        peak = max(peak, c.peak)
    gain = level / peak

    avg = [it.average for it in stats]
    for i, s in enumerate(frag):
        z = ()
        for c, c_avg in zip(s, avg):
            z += ((c - c_avg) * gain,)
        frag[i] = z
