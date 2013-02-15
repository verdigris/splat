import os
import collections
import struct
import wave
import _geomusic

class AudioFile(collections.Sequence):
    """Audio file factory"""

    @staticmethod
    def open(file_name, mode='r'):
        ext = file_name.rpartition('.')[2]
        if ext == 'wav':
            return WaveFile(file_name, mode)
        raise Exception('Unsupported file format: {0}'.format(ext))

    def __init__(self, file_name, mode):
        raise NotImplementedError

    def close(self):
        pass

    @property
    def channels(self):
        raise NotImplementedError

    @property
    def duration(self):
        raise NotImplementedError

    @property
    def sample_rate(self):
        raise NotImplementedError

    def get_bytearray(self):
        raise NotImplementedError


class WaveFile(AudioFile):
    def __init__(self, file_name, mode):
        self._file = wave.open(file_name, mode)
        self._sampwidth = self._file.getsampwidth()
        if self._sampwidth != 2:
            raise Exception("Unsupported sample width")
        self._channels = self._file.getnchannels()
        self._length = self._file.getnframes()
        self._fmt = '<' + ('h' * self._channels)

    def close(self):
        self._file.close()

    def __getitem__(self, n):
        if n >= self._length:
            raise IndexError
        self._file.setpos(n)
        data = struct.unpack(self._fmt, self._file.readframes(1))
        return tuple((float(data[i]) / 32768.0 for i in range(self._channels)))

    def __len__(self):
        return self._length

    @property
    def channels(self):
        return self._channels

    @property
    def sample_width(self):
        return self._sampwidth

    @property
    def duration(self):
        return self._file.getnframes() / float(self._file.getframerate())

    @property
    def sample_rate(self):
        return self._file.getframerate()

    # ToDo: pass array to avoid re-allocating each time
    def get_bytearray(self, start=0, length=None):
        if length is None:
            length = len(self)
        self._file.setpos(start)
        return bytearray(self._file.readframes(length))


class Fragment(_geomusic.Fragment):
    @classmethod
    def open(cls, file_name, mode='r'):
        file = AudioFile.open(file_name, mode)
        frag = cls(file.channels, file.sample_rate, file.duration)
        rem = len(frag)
        cur = 0
        chunk_size = (64 * 1024) / (file.sample_width * file.channels)
        while rem > 0:
            n = chunk_size if (rem >= chunk_size) else rem
            bytes = file.get_bytearray(cur, n)
            frag.import_bytes(bytes, cur, file.sample_width,
                              file.sample_rate, file.channels)
            rem -= n
            cur += n
        file.close
        return frag

    def resize(self, duration):
        n = int(duration * self.sample_rate)
        self._resize(n)
        return n

    def n2s(self, n):
        return float(n) / self.sample_rate

    def s2n(self, s):
        return int(s * self.sample_rate)

    def save_to_file(self, file_name, sample_width, start=0, end=None):
        f = wave.open(file_name, 'w')
        f.setnchannels(self.channels)
        f.setsampwidth(sample_width)
        f.setframerate(self.sample_rate)
        f.setnframes(len(self))
        bytes = self.as_bytes(sample_width)
        k = self.channels * 2
        start_n = start * k
        if end is not None:
            end_n = end * k
        else:
            end_n = len(self) * k
        f.writeframes(bytes[start_n:end_n])
        f.close()
