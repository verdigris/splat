import os
import wave
import _geomusic

# ToDo: make Fragment from file contents
def open_audio_file(file_name, mode):
    ext = file_name.rpartition('.')[2]
    if ext == 'wav':
        return wave.open(file_name, mode)
    raise Exception('Unsupported file format: {0}'.format(ext))


class Fragment(_geomusic.Fragment):
    def normalize(self, level=1.0):
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

        stats = [Stats() for x in xrange(self.channels)]
        for s in self:
            for c, c_stats in zip(s, stats):
                c_stats.process(c)

        peak = 0.0
        for c in stats:
            peak = max(peak, c.peak)
        gain = level / peak

        avg = [it.average for it in stats]
        for i, s in enumerate(self):
            z = ()
            for c, c_avg in zip(s, avg):
                z += ((c - c_avg) * gain,)
            self[i] = z

    def resize(self, duration):
        n = int(duration * self.sample_rate)
        self._resize(n)
        return n

    def save_to_file(self, file_name, sample_width):
        print("sample width: {0}, channels: {1}, length: {2}".format(
                sample_width, self.channels, len(self)))
        f = wave.open(file_name, 'w')
        f.setnchannels(self.channels)
        f.setsampwidth(sample_width)
        f.setframerate(self.sample_rate)
        f.setnframes(len(self))
        f.writeframes(self.as_bytes(sample_width))
        f.close()
