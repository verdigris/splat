# Splat - splat/data.py
#
# Copyright (C) 2012, 2013 Guillaume Tucker <guillaume@mangoz.org>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import struct
import wave
import _splat

def open_wav(file_name):
    if file_name.rpartition('.')[2].lower() != 'wav':
        return None

    w = wave.open(file_name, 'rb')
    channels = w.getnchannels()
    n_frames = w.getnframes()
    sample_rate = w.getframerate()
    duration = n_frames / float(sample_rate)
    sample_width = w.getsampwidth()

    frag = Fragment(channels, sample_rate, duration)
    rem = len(frag)
    cur = 0
    chunk_size = (64 * 1024) / (sample_width * channels)

    while rem > 0:
        n = chunk_size if (rem >= chunk_size) else rem
        raw_bytes = bytearray(w.readframes(n))
        frag.import_bytes(raw_bytes, cur, sample_width, sample_rate, channels)
        rem -= n
        cur += n

    w.close()

    return frag


audio_file_openers = [open_wav,]


class Fragment(_splat.Fragment):

    """A fragment of sound data.

    Create an empty sound fragment with the given number of ``channels``,
    sample ``rate`` in Hertz and initial ``duration`` in seconds.  If no
    ``duration`` is specified, the fragment will be empty.  All the samples are
    initialised to 0 (silence).

    All Splat sound data is contained in ``Fragment`` objects.  They are
    accessible as a mutable sequence of tuples of floating point values to
    represent the samples of the audio channels.  The length of each sample
    tuple is equal to the number of channels of the fragment, the maximum being
    fixed to 16.
    """

    @classmethod
    def open(cls, file_name):
        """Open a file to create a sound fragment

        Open a sound file specified by ``file_name`` and imports its contents
        into a new ``Fragment`` instance, which is then returned.
        """
        for opener in audio_file_openers:
            frag = opener(file_name)
            if frag is not None:
                return frag

        raise Exception("Unsupported file format")

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
