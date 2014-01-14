# Splat - splat/data.py
#
# Copyright (C) 2012, 2013, 2014 Guillaume Tucker <guillaume@mangoz.org>
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

def file_ext(file_name):
    return file_name.rpartition(os.extsep)[2].lower()

# -----------------------------------------------------------------------------
# Audio file openers

def open_wav(file_name):
    if file_ext(file_name) != 'wav':
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

# -----------------------------------------------------------------------------
# Audio file savers

def save_wav(file_name, frag, start, end, sample_width=2):
    w = wave.open(file_name, 'w')
    w.setnchannels(frag.channels)
    w.setsampwidth(sample_width)
    w.setframerate(frag.sample_rate)
    w.setnframes(len(frag))
    raw_bytes = frag.as_bytes(sample_width)
    k = frag.channels * sample_width
    start_n = start * k
    end_n = end * k
    w.writeframes(raw_bytes[start_n:end_n])
    w.close()

audio_file_savers = { 'wav': save_wav, }

# -----------------------------------------------------------------------------
# Data classes

class Fragment(_splat.Fragment):

    """A fragment of sound data.

    Create an empty sound fragment with the given number of ``channels``,
    sample ``rate`` in Hertz and initial ``duration`` in seconds.  The default
    number of channels is 2 and the default sample rate is 48000.  If no
    duration is specified, the fragment will be empty.  All the samples are
    initialised to 0 (silence).

    All Splat sound data is contained in :py:class:`splat.data.Fragment`
    objects.  They are accessible as a mutable sequence of tuples of floating
    point values to represent the samples of the audio channels.  The length of
    each sample tuple is equal to the number of channels of the fragment, the
    maximum being fixed to 16.

    .. note::

       It is rather slow to access and modify all the samples of a Fragment
       using the standard Python sequence interface.  The underlying
       implementation in C can handle all common data manipulations much
       faster using the Fragment methods documented below.
    """

    @classmethod
    def open(cls, file_name):
        """Open a file to create a sound fragment by importing audio data.

        Open a sound file specified by ``file_name`` and import its contents
        into a new :py:class:`splat.data.Fragment` instance, which is then
        returned.  Only a limited set of formats are supported (currently only
        ``wav``).  All the samples are converted to floating point values.
        """
        for opener in audio_file_openers:
            frag = opener(file_name)
            if frag is not None:
                return frag

        raise Exception("Unsupported file format")

    def n2s(self, n):
        """Convert a sample index number ``n`` into a time in seconds."""
        return float(n) / self.sample_rate

    def s2n(self, s):
        """Convert a time in seconds ``s`` into a sample index number."""
        return int(s * self.sample_rate)

    def save(self, file_name, fmt=None, start=0, end=None, *args, **kw):
        """Save the contents of the audio fragment into a file.

        A file called ``file_name`` is created and the contents of the audio
        fragment are written to it.  It is possible to save only a part of the
        fragment using the ``start`` and ``end`` arguments with times in
        seconds.

        The ``fmt`` argument is a string to identify the output file format to
        use.  If ``None``, the file name extension is used.  Currently only
        ``wav`` is supported.  Extra arguments that are specific to the file
        format can be added.

        The ``wav`` format accepts an extra ``sample_width`` argument to
        specify the number of bytes per sample for each channel, which is 2 by
        default (16 bits).
        """
        if fmt is None:
            fmt = file_ext(file_name)
        saver = audio_file_savers.get(fmt, None)
        if saver is None:
            raise Exception("Unsupported file format: {0}".format(fmt))
        if end is None:
            end_n = len(self)
        else:
            end_n = int(end * self.sample_rate)
        start_n = int(start * self.sample_rate)
        saver(file_name, self, start_n, end_n, *args, **kw)
