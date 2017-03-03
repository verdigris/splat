# Splat - splat/data.py
#
# Copyright (C) 2012, 2013, 2014, 2015 Guillaume Tucker <guillaume@mangoz.org>
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
import md5
try:
    # For extra standard audio formats
    import audiotools
    has_audiotools = True
    try:
        from cStringIO import StringIO
    except ImportError:
        from StringIO import StringIO
except ImportError:
    has_audiotools = False
import _splat
import splat

# Used in Splat Audio Fragment files (.saf)
SAF_MAGIC = 'Splat!'
SAF_FORMAT = 1

# -----------------------------------------------------------------------------
# Utilities

def _get_fmt(f, fmt):
    if fmt is None:
        if isinstance(f, str):
            fmt = f.rpartition(os.extsep)[2].lower()
        else:
            raise Exception("Format required with file objects")
    return fmt

def _read_chunks(frag, read_frames, frame_size, sample_type):
    rem = len(frag)
    cur = 0
    chunk_size = 65536 / frame_size # read in blocks of 64K

    while rem > 0:
        n = min(chunk_size, rem)
        data = bytearray(read_frames(n))
        frag.import_bytes(data, frag.rate, frag.channels, sample_type, cur)
        rem -= n
        cur += n

# -----------------------------------------------------------------------------
# Audio file openers

def open_wav(wav_file, fmt=None):
    if _get_fmt(wav_file, fmt) != 'wav':
        return None

    w = wave.open(wav_file, 'rb')
    channels = w.getnchannels()
    n_frames = w.getnframes()
    rate = w.getframerate()
    sample_width_bytes = w.getsampwidth()
    frame_size = sample_width_bytes * channels
    s_width = sample_width_bytes * 8
    frag = Fragment(channels, rate, length=n_frames)
    _read_chunks(frag, w.readframes, frame_size, 'int{:d}'.format(s_width))
    w.close()
    return frag

def open_saf(saf_file, fmt=None):
    if _get_fmt(saf_file, fmt) != 'saf':
        return None

    is_str = isinstance(saf_file, str)
    f = open(saf_file, 'rb') if is_str else saf_file
    if f.readline().strip('\n') != SAF_MAGIC:
        return None
    attr_str = f.readline().strip('\n')
    attr = {k: v for k, v in (kv.split('=') for kv in attr_str.split(' '))}
    opts = ['rate', 'channels', 'length', 'precision', 'format']
    rate, channels, length, precision, fmt = (int(attr[x]) for x in opts)
    if fmt > SAF_FORMAT:
        raise Exception("Format is more recent than this version of Splat")
    frag = Fragment(rate=rate, channels=channels, length=length)
    frame_size = channels * precision / 8
    _read_chunks(frag, (lambda x: f.read(x * frame_size)), frame_size,
                 sample_type='float{:d}'.format(precision))
    if is_str:
        f.close()
    if frag.md5() != attr['md5']:
        raise Exception("Fragment MD5 sum mismatch")
    return frag

audio_file_openers = [open_wav, open_saf,]

if has_audiotools is True:
    def open_audiotools(f_name, fmt=None):
        if not isinstance(f_name, str):
            return None
        try:
            f = audiotools.open(f_name)
        except audiotools.UnsupportedFile:
            return None
        sample_width = f.bits_per_sample()
        pcm = f.to_pcm()
        frag = Fragment(rate=f.sample_rate(), channels=f.channels(),
                        length=f.total_frames())
        frame_size = frag.channels * sample_width / 8
        chunk_size = 65536 / frame_size
        rem = f.total_frames()
        cur = 0
        while rem > 0:
            frames = pcm.read(chunk_size)
            if frames.frames == 0:
                # Some lossy formats like MP3 have slightly varying data length
                break
            raw_bytes = bytearray(frames.to_bytes(False, True))
            frag.import_bytes(raw_bytes, frag.rate, frag.channels,
                              'int{:d}'.format(sample_width), cur)
            rem -= frames.frames
            cur += frames.frames
        return frag
    audio_file_openers.append(open_audiotools)

# -----------------------------------------------------------------------------
# Audio file savers

def save_wav(wav_file, frag, start, end, sample_type='int16'):
    sample_width = splat.sample_types[sample_type]
    w = wave.open(wav_file, 'w')
    w.setnchannels(frag.channels)
    w.setsampwidth(sample_width / 8)
    w.setframerate(frag.rate)
    w.setnframes(len(frag))
    raw_bytes = frag.export_bytes(sample_type, start, end)
    w.writeframes(buffer(raw_bytes))
    w.close()

def save_saf(saf_file, frag, start, end):
    is_str = isinstance(saf_file, str)
    f = open(saf_file, 'w') if is_str else saf_file
    raw_bytes = frag.export_bytes(start=start, end=end)
    frame_size = frag.channels * splat.SAMPLE_WIDTH / 8
    length = len(raw_bytes) / frame_size
    md5sum = md5.new(raw_bytes).hexdigest()
    attrs = {
        'version': splat.VERSION_STR,
        'build': splat.BUILD,
        'format': SAF_FORMAT,
        'channels': frag.channels,
        'rate': frag.rate,
        'precision': splat.SAMPLE_WIDTH,
        'length': length,
        'md5': md5sum,
        }
    h = ' '.join('='.join(str(x) for x in kv) for kv in attrs.iteritems())
    f.write(SAF_MAGIC + '\n')
    f.write(h + '\n')
    f.write(buffer(raw_bytes))
    if is_str:
        f.close()

audio_file_savers = { 'wav': save_wav, 'saf': save_saf, }

if has_audiotools is True:
    channel_masks = [0x4, 0x3, 0x7, 0x33, 0x37, 0x137, 0x637, 0x737, 0xF37,
                     0x7F7, 0xFF7, 0x2FF7, 0x6FF7, 0x7FF7, 0x17FF7, 0x2FFF7]

    def save_audiotools(cls, fname, frag, start, end, sample_type='int16',
                        channel_mask=None, compression=None):
        if not isinstance(fname, str):
            raise TypeError('File name needed when saving with audiotools')
        if channel_mask is None:
            channel_mask = channel_masks[frag.channels]
        data = frag.export_bytes(sample_type, start, end)
        str_io = StringIO(str(data))
        sample_width = splat.sample_types[sample_type]
        reader = audiotools.PCMReader(str_io, frag.rate, frag.channels,
                                      channel_mask, sample_width)
        cls.from_pcm(fname, reader, compression, len(frag))

    fmt_cls = {
        'ogg': audiotools.VorbisAudio,
        'flac': audiotools.FlacAudio,
        'mp3': audiotools.MP3Audio,
        }
    for fmt, cls in fmt_cls.iteritems():
        audio_file_savers[fmt] = lambda *a, **k: save_audiotools(cls, *a, **k)

# -----------------------------------------------------------------------------
# Data classes

class Fragment(_splat.Fragment):
    """A fragment of sound data

    Create an empty sound fragment with the given number of ``channels``,
    sample ``rate`` in Hertz and initial ``duration`` in seconds or ``length``
    in number of samples per channel.  The default number of channels is 2, the
    minimum is 1 and the maximum is currently fixed at 16.  The default sample
    rate is 48000.  If no duration is specified, the fragment will be empty.
    Otherwise, the samples are initialised to 0 (silence).

    The fragment can also be allocated using ``mmap``.  If set to ``True``, the
    fragment will be backed with a temporary mmap file.  It can also be a
    string used as a path prefix to create persistent mmap files, or a list of
    paths with existing files to reuse.  In this case, there must be the same
    number of paths as of ``channels`` and the ``length`` should be set to
    match the mmap file sizes as new fragments are always resized to the given
    length.

    All Splat sound data is contained in :py:class:`splat.data.Fragment`
    objects.  They are accessible as a mutable sequence of tuples of floating
    point values to represent the samples across all the audio channels.  The
    length of each sample tuple is equal to the number of channels of the
    fragment.

    .. note::

       It is rather slow to access and modify all the samples of a Fragment
       using the standard Python sequence interface.  The underlying
       implementation in C can handle all common data manipulations much
       faster using the Fragment methods documented below.
    """

    @classmethod
    def open(cls, in_file, fmt=None, name=None):
        """Open a file to create a sound fragment by importing audio data.

        Open a sound file specified by ``in_file``, which can be either a file
        name or a file-like object, and import its contents into a new
        :py:class:`splat.data.Fragment` instance, which is then returned.

        The format can be specified via the ``fmt`` argument, otherwise it may
        be automatically detected whenever possible.  Only a limited set of
        :ref:`audio_files` formats are supported (currently only ``wav`` and
        ``saf``).  All the samples are converted to floating point values with
        the native precision of the fragment.
        """
        fmt = _get_fmt(in_file, fmt)
        for opener in audio_file_openers:
            frag = opener(in_file, fmt)
            if frag is None:
                continue
            if name is not None:
                frag.name = name
            elif isinstance(in_file, str):
                frag.name = os.path.splitext(os.path.basename(in_file))[0]
            return frag
        raise Exception("Unsupported file format")

    def save(self, out_file, fmt=None, start=None, end=None, normalize=True,
             *args, **kw):
        """Save the contents of the audio fragment into a file.

        If ``out_file`` is a string, a file is created with this name;
        otherwise, it must be a file-like object.  The contents of the audio
        fragment are written to this file.  It is possible to save only a part
        of the fragment using the ``start`` and ``end`` arguments with times in
        seconds.  The contents of the fragment will be automatically normalized
        unless ``normalized`` is set to ``False``.

        The ``fmt`` argument is a string to identify the output file format to
        use.  If ``None``, the file name extension is used.  It may otherwise
        be a file-like object, in which case the format needs to be specified.
        Extra arguments that are specific to the file format may be added.

        The ``wav`` format accepts an extra ``sample_type`` argument to specify
        the sample type, which also defines the number of bytes per sample.
        The default is ``int16`` (16 bits).
        """
        fmt = _get_fmt(out_file, fmt)
        saver = audio_file_savers.get(fmt, None)
        if saver is None:
            raise Exception("Unsupported file format: {0}".format(fmt))
        if normalize is True:
            self.normalize()
        saver(out_file, self, start, end, *args, **kw)

    def dup(self):
        """Duplicate this fragment into a new one and return it."""
        dup_frag = Fragment(channels=self.channels, rate=self.rate)
        dup_frag.mix(self)
        return dup_frag

    def grow(self, duration=None, length=None):
        """Resize if ``duration`` or ``length`` is greater than current."""
        if duration is not None:
            if duration > self.duration:
                self.resize(duration=duration)
        elif length is not None:
            if length > len(self):
                self.resize(length=length)
        else:
            raise ValueError("neither new duration nor length was supplied")

    def md5(self, sample_type=splat.SAMPLE_TYPE, as_md5_obj=False):
        """Get the MD5 checksum of this fragment's data.

        The data is first converted to samples as specified by ``sample_type``
        and ``sample_width`` in bytes.  Then the MD5 checksum is returned as a
        string unless ``as_md5_obj`` is set to True in which case an ``md5``
        object is returned instead."""
        md5sum = md5.new(self.export_bytes(sample_type))
        return md5sum if as_md5_obj is True else md5sum.hexdigest()

    def n2s(self, n):
        """Convert a sample index number ``n`` into a time in seconds."""
        return float(n) / self.rate

    def s2n(self, s):
        """Convert a time in seconds ``s`` into a sample index number."""
        return int(s * self.rate)
