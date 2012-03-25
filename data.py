import os
import wave

def open_audio_file(file_name, mode):
    ext = file_name.rpartition('.')[2]
    if ext == 'wav':
        return wave.open(file_name, mode)
    raise Exception('Unsupported file format: {0}'.format(ext))

class Fragment(object):
    def __init__(self, channels, sample_rate):
        self._sample_rate = sample_rate
        self.data = [list() for i in range(channels)]

    @property
    def sample_rate(self):
        return self._sample_rate

    @property
    def length(self):
        return len(self.data[0])

    @property
    def duration(self):
        return self.length / self.sample_rate

    @property
    def channels(self):
        return len(self.data)

    def resize(self, duration):
        n = int(duration * self.sample_rate)
        self._resize(n)
        return n

    def _resize(self, n):
        if n > self.length:
            padding = [0.0 for i in range(n - self.length)]
            for channel in self.data:
                channel += padding

    def mix(self, frag, start):
        if frag.channels != self.channels:
            raise Exception("Channels mismatch")
        start_sample = int(start * self.sample_rate)
        self._resize(start_sample + frag.length)
        for master, mix in zip(self.data, frag.data):
            for i, s in enumerate(mix):
                master[i + start_sample] += s

    def save_to_file(self, file_name, sample_width):
        print("sample width: {0}, channels: {1}, length: {2}".format(
                sample_width, self.channels, self.length))
        f = wave.open(file_name, 'w')
        f.setnchannels(self.channels)
        f.setsampwidth(sample_width)
        f.setframerate(self.sample_rate)
        f.setnframes(self.length)
        f.writeframes(self.get_raw_bytes(sample_width))
        f.close()

    def get_raw_bytes(self, sample_width):
        if sample_width != 2:
            raise Exception("only 16 bits for now...")
        raw_data = list()
        for i in range(self.length):
            for channel in self.data:
                sample = channel[i]
                if abs(sample) > 1.0:
                    sample /= abs(sample)
                raw_sample = int(sample * 32767)
                raw_data.append(raw_sample & 0xFF)
                raw_data.append((raw_sample & 0xFF00) >> 8)
        return bytearray(raw_data)
