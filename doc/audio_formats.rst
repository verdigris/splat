Audio formats
=============

.. _sample_formats:

Audio sample types
------------------

Audio samples in Splat can be signed integers or floating point numbers of
different widths.  The :py:class:`splat.data.Fragment` objects internally deal
with either 32-bit or 64-bit floating point samples, the precision being
defined when the ``_splat`` C extension is compiled.  Sample types are defined
in Splat with the following names (plain strings), which tell the numeric
format (integer or float) and the width in number of bits:

* ``int8``
* ``int16``
* ``int24``
* ``float32``
* ``float64``

The native sample type is defined by ``splat.SAMPLE_TYPE`` and is usually
``float64`` or ``float32`` in some special cases.  As a convenience, the width
of the native sample type is defined by ``splat.SAMPLE_WIDTH``.  To
automatically determine the width of any given sample type, or to list all the
available sample types, the ``splat.sample_types`` dictionary can be used.

.. _audio_files:

Native audio file formats (SAF and WAV)
---------------------------------------

The following names are used to choose the audio file format in
:py:meth:`splat.data.Fragment.open` and :py:meth:`splat.data.Fragment.save`:

``saf``
  The Splat Audio Fragment file format is primarily made to export and import
  fragment objects without losing any of the original precision.  By default,
  they contain floating point samples with the native sample width.  It's also
  possible to explicitely specify the sample width which may result in a
  conversion between 32-bit and 64-bit types with a potential loss of
  information.  The ``saf`` format is useful when building complex splats to
  avoid having to regenerate everything from scratch each time the code is run,
  or to share some data between different splats.  There is currently no known
  program to play these files directly, although it's quite easy to import one
  into a fragment and directly export it again as WAV.

  When using fragments allocated with mmap and backed by persistent files, a
  variant of the ``saf`` format is used.  It will contain the meta-data
  including the names of separate mmap files for each channel of the fragment.
  This is entirely transparent from a user point of view, only that the
  individual mmap files must not be renamed and need to be kept in the same
  directory as the meta-data ``saf`` file.  The meta-data and mmap files can be
  moved together.

``wav``
  This is for standard WAV audio files.  It uses the standard ``wave`` Python
  module.  Any sample width value supported by the library can be used when
  saving into a WAV file, but only 8-bit and 16-bit integer samples can be used
  when opening a WAV file in Splat.


Extra file formats with ``audiotools``
--------------------------------------

When the `audiotools <http://audiotools.sourceforge.net/>`_ package is
installed, additional formats such as ``flac``, ``ogg`` or ``mp3`` are also
supported to import and export data.  With audiotools, only file names can be
used and anonymous file-like objects are not supported unlike with ``wav`` and
``saf`` formats.

This makes it possible to open and save files using fragments, and mix them all
together as in this example::

  import splat.data
  import splat.gen
  from splat import dB

  # Create a generator and generate a sound
  gen = splat.gen.SineGenerator()
  gen.run(0.0, 1.0, 440.0)

  # Open 2 files for different formats (MP3, SAF) into new fragments
  sample1 = splat.data.Fragment.open('sample1.mp3')
  sample2 = splat.data.Fragment.open('sample2.saf')

  # Create master fragment to mix everything together
  master = splat.data.Fragment()

  # Mix the generated fragment at 2.0s
  master.mix(gen.frag, offset=2.0)

  # Mix the first sample at 1.5s skipping the first 0.5s
  master.mix(sample1, offset=1.5, skip=0.5)

  # Mix the second sample at 0.3s with -0.5 dB gain
  master.mix(sample2, offset=0.3, levels=dB(-0.5))

  # Save result to output file as FLAC
  master.save('example.flac')
