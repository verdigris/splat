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

    """Wave file implementation of :py:class:`geomusic.data.AudioFile`"""

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

    """A fragment of sound data.

    Create an empty sound fragment with the given number of ``channels``,
    sample ``rate`` in Hertz and initial ``duration`` in seconds.  If no
    ``duration`` is specified, the fragment will be empty.  All the samples are
    initialised to 0 (silence).

    All Geomusic sound data is contained in ``Fragment`` objects.  They are
    accessible as a mutable sequence of tuples of floating point values to
    represent the samples of the audio channels.  The length of each sample
    tuple is equal to the number of channels of the fragment, the maximum being
    fixed to 16.
    """

    @classmethod
    def open(cls, file_name, mode='r'):
        """Open a file to create a sound fragment

        Open a sound file specified by ``file_name`` and imports its contents
        into a new ``Fragment`` instance, which is then returned.
        """
        af = AudioFile.open(file_name, mode)
        frag = cls(af.channels, af.sample_rate, af.duration)
        rem = len(frag)
        cur = 0
        chunk_size = (64 * 1024) / (af.sample_width * af.channels)
        while rem > 0:
            n = chunk_size if (rem >= chunk_size) else rem
            raw_bytes = af.get_bytearray(cur, n)
            frag.import_bytes(raw_bytes, cur, af.sample_width,
                              af.sample_rate, af.channels)
            rem -= n
            cur += n
        af.close
        return frag

    def resize(self, duration):
        """Resize the fragment to the specified ``duration`` in seconds."""
        n = int(duration * self.sample_rate)
        self._resize(n)
        return n

    def n2s(self, n):
        """Convert a sample index number ``n`` into a time in seconds."""
        return float(n) / self.sample_rate

    def s2n(self, s):
        """Convert a time in seconds ``s`` to a sample index number."""
        return int(s * self.sample_rate)

    def save_to_file(self, file_name, sample_width=2, start=0, end=None):
        """Save the contents of the audio fragment to an audio file.

        A file is created with the provided name ``file_name``, and the
        contents of the audio fragment are written to it.  The ``sample_width``
        argument contains the resolution in bytes for each channels, which by
        default is 2 (16 bits).

        Note: Only ``WAV`` files are currently supported.

        It's also possible to only save a part of the fragment using the
        ``start`` and ``end`` arguments with times in seconds.
        """
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
