.. _sample_formats:

Audio sample formats
^^^^^^^^^^^^^^^^^^^^

Audio samples in Splat can be either signed integer or floating point numbers.
The :py:class:`splat.data.Fragment` objects internally deal with either 32-bit
or 64-bit samples numbers, the precision being defined when the ``_splat`` C
extension is compiled.  Sample formats are defined in Splat with the following
attributes (names used in functions arguments):

``sample_type``
  Type of value used in samples, can be either:

* ``splat.SAMPLE_INT`` for signed integers
* ``splat.SAMPLE_FLOAT`` for floating point numbers

``sample_width``
  Width as a number of bits in each sample.

The following format combinations are currently supported:

+------------------------+--------------+
| Sample type            | Sample width |
+------------------------+--------------+
| ``splat.SAMPLE_INT``   |            8 |
+------------------------+--------------+
| ``splat.SAMPLE_INT``   |           16 |
+------------------------+--------------+
| ``splat.SAMPLE_INT``   |           24 |
+------------------------+--------------+
| ``splat.SAMPLE_FLOAT`` |           32 |
+------------------------+--------------+
| ``splat.SAMPLE_FLOAT`` |           64 |
+------------------------+--------------+

The native sample format is 64-bit floating point; this is used as the default
format when exporting raw data from a fragment.


.. _audio_files:

Audio file formats
^^^^^^^^^^^^^^^^^^

The following names are used to choose the audio file format in
:py:meth:`splat.data.Fragment.open` and :py:meth:`splat.data.Fragment.save`:

``wav``
  This is for standard WAV audio files, with 8-bit or 16-bit integer samples
  and no compression.  It uses the standard ``wave`` Python module.  Any sample
  width value supported by the library can be used, but only 8-bit and 16-bit
  samples can be used when opening a WAV file in Splat.

``saf``
  The Splat Audio Fragment file format is primarily made to export and import
  Fragment objects without loosing any of the original samples precision.  By
  default, they contain floating point samples with the native sample width.
  It's also possible to explicitely specify the sample width which may result
  in a conversion between 32-bit and 64-bit types with potential loss of
  precision.  The ``saf`` format is useful when building complex splats to
  avoid having to regenerate everything from scratch each time the code is run,
  or to share some data across other splats.  There is currently no known
  programme to play these files directly, although it's quite easy to import
  one into a Fragment and export it again as WAV.


Additional audio file formats
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When the `audiotools <http://audiotools.sourceforge.net/>`_ package is
installed, additional formats such as ``flac``, ``ogg`` or ``mp3`` are also
supported to import and export data.  With audiotools, only file names can be
used and anonymous file-like objects are not supported unlike with ``wav`` and
``saf`` formats.
