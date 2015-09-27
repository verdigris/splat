.. _utilities:

Common utilities
================

This section describes common utility functions and classes that are
mentioned and used in examples throughout this manual.

Conversion functions
--------------------

.. autofunction:: splat.lin2dB
.. autofunction:: splat.dB2lin

There is also a shorthand :py:meth:`splat.dB` alias for
:py:meth:`splat.dB2lin`.


Signal objects
--------------

.. autoclass:: splat.Signal
   :members:


.. _fragments:

Fragment objects
================

.. automodule:: splat.data

.. autoclass:: splat.data.Fragment(channels, rate, duration=0.0, length=0)
   :members:

   .. automethod:: splat.data.Fragment.mix
   .. automethod:: splat.data.Fragment.import_bytes
   .. automethod:: splat.data.Fragment.export_bytes
   .. automethod:: splat.data.Fragment.get_peak
   .. automethod:: splat.data.Fragment.normalize
   .. automethod:: splat.data.Fragment.amp
   .. automethod:: splat.data.Fragment.lin2dB
   .. automethod:: splat.data.Fragment.dB2lin
   .. automethod:: splat.data.Fragment.offset
   .. automethod:: splat.data.Fragment.resize
   .. automethod:: splat.data.Fragment.resample
   .. autoattribute:: splat.data.Fragment.rate
   .. autoattribute:: splat.data.Fragment.duration
   .. autoattribute:: splat.data.Fragment.channels
   .. autoattribute:: splat.data.Fragment.name


.. _filters:

Filters and FilterChain objects
===============================

.. automodule:: splat.filters

A *filter function* takes a :py:class:`splat.data.Fragment` and a tuple of
arguments with its specific parameters.  It is expected to run on the entirety
of the fragment.  Filter functions can be combined into a series via the
:py:class:`splat.filters.FilterChain` class.

.. autoclass:: splat.filters.FilterChain
   :members:

.. autofunction:: splat.filters.linear_fade
.. autofunction:: splat.filters.dec_envelope
.. autofunction:: splat.filters.reverse
.. autofunction:: splat.filters.reverb


.. _generators:

Generator objects
=================

.. automodule:: splat.gen

.. autoclass:: splat.gen.Generator
   :members:
   :private-members:


Source generators
-----------------

.. autoclass:: splat.gen.SourceGenerator
   :members:


.. _sources:

Sound sources
-------------

.. automodule:: splat.sources
.. autofunction:: splat.sources.sine
.. autofunction:: splat.sources.square
.. autofunction:: splat.sources.triangle
.. autofunction:: splat.sources.overtones


Source generator objects
------------------------

.. autoclass:: splat.gen.SineGenerator
   :members:

.. autoclass:: splat.gen.SquareGenerator
   :members:

.. autoclass:: splat.gen.TriangleGenerator
   :members:

.. autoclass:: splat.gen.OvertonesGenerator
   :members:


Particle generators
-------------------

.. autoclass:: splat.gen.Particle
   :members:

.. autoclass:: splat.gen.ParticlePool
   :members:

.. autoclass:: splat.gen.ParticleGenerator
   :members:


.. _scales:

Scales
======

.. automodule:: splat.scales

The main purpose of a :py:class:`splat.scales.Scale` object is to calculate a
frequency for a given note name.  This creates a bridge between musical
composition and sound generators.  The classic equi-tempered scale is
implemented with the :py:class:`splat.scales.LogScale` but endless other
implementations can be made to freely define any scale behaviour via the
:py:meth:`splat.scales.Scale.get_freq` method.  This can even be made dynamic,
for example by associating a varying frequency to each note based on the
sequence of previous notes or the current point in time within a composition.

A scale can then be used like a dictionary, mapping note names with
frequencies.  For example:

.. code-block:: python

    import splat.scales

    s = splat.scales.LogScale()
    a_freq = s['A']

should result in 440.0 Hz for the value of ``a_freq`` with the default scale
parameters.  This standard interface should make it possible to use any scale
implementation without rewriting any of an existing musical composition as long
as the scales use the same :ref:`note_names` (i.e. 'A', 'Bb', 'B', 'C', 'C#'
etc...).  It's also possible to define a different naming convention for the
notes when implementing a new scale.

A Scale object has the following characteristics:

pitch
  Reference frequency of the first note, usually called A, for all scales.
key
  Name of the note that identifies the first note of the scale.
steps
  Number of notes in the scale, for example 12 for the classic semi-tones.
base
  Coefficient to multiply a note frequency to get to the next cycle, 2 for a
  classic octave.
cycle
  General term for what an octave is to a classic scale, all the steps of a
  scale repeat themselves in each cycle.


.. _note_names:

Note names
----------

Classic note names use the letters 'A' to 'G', with '#' and 'b' alterations for
'sharp' and 'flat' respectively.  Then a number can be appended to use a
different octave, or generally speaking a different *cycle*, like 'C4'.  This
notation is used by classic scales like the :py:class:`splat.scales.LogScale`
and :py:class:`splat.scales.HarmonicScale` (via the
:py:meth:`splat.scales.Scale.get_note` method), but can be freely overridden in
other scale implementations.  One important thing to be aware of is the
**relative mode** which makes the octave number relative to the key:

A
  Note A, which frequency is equal to the pitch (i.e. 440Hz by default).
A1
  Note A an octave higher, which corresponds to 880Hz by default when the key
  is A.
Eb
  'E-flat' in the same octave.
G#
  'G-sharp' in the same octave.
C0
  Like C, in the same octave.
C1
  C an octave higher.
F-1
  F an octave lower.
Ab3
  'A-flat' three octaves higher.

It's also possible to use **absolute notation** to make the octave number
independent of the key being used, using the
:py:attr:`splat.scales.Scale.abs_cycles` attribute.  In absolute mode, the A
pitch frequency corresponds to octave number 4, so 'A' is the same as 'A4' and
is an octave lower than 'A5'.  This does not only change the resulting octave
(or cycle) number but also determines which octave a given note belongs to.
The tables below show concrete examples in both modes, using a LogScale and
showing in bold the first note for each key:

.. include:: notes_table.rst

Improved Beep
-------------

And now, here's a small improvement on the :ref:`beep` example using a scale:

.. code-block:: python

    import splat.gen
    import splat.filters
    import splat.scales

    gen = splat.gen.SineGenerator()
    gen.filters = [splat.filters.linear_fade]
    s = splat.scales.LogScale(key='A')
    gen.run(0.0, 1.0, s['A'])
    gen.frag.save("A-note.wav")


.. _scale_objects:

Scale objects
-------------

.. autoclass:: splat.scales.Scale
   :members:

.. autoclass:: splat.scales.LogScale
   :members:

.. autoclass:: splat.scales.HarmonicScale
   :members:


.. _sequencer_objects:

Sequencers
==========

.. automodule:: splat.seq

A sequencer can be used to control when sounds are generated in
time. :py:class:`splat.seq.PatternSequencer` objects periodically run a series
of :py:class:`splat.seq.Pattern` objects, calling them for each beat in each
bar so they can generate sounds at regular intervals and mix them into a main
fragment.  For patterns that use existing sounds from files or other generated
fragments, a :py:class:`splat.seq.SampleSet` object can be used to group
equivalent sounds together and pick one randomly.

.. autoclass:: splat.seq.SampleSet
   :members:

.. autoclass:: splat.seq.Pattern
   :members:

.. autoclass:: splat.seq.FuzzyPattern
   :members:

.. autoclass:: splat.seq.Silence
   :members:

.. autoclass:: splat.seq.Sequencer
   :members:

.. autoclass:: splat.seq.PatternSequencer
   :members:


.. _interpolation:

Interpolation and splines
=========================

.. automodule:: splat.interpol

The :py:mod:`splat.interpol` module contains a set of classes designed to work
together and provide polynomial interpolation functionality.  The principle is
to build a continuous function for a given set of input discrete coordinates.
A :py:class:`splat.interpol.Polynomial` object represents a polynomial function
with a series of coefficients.  It can be calculated by reducing a matrix
containing some data points coordinates using
:py:class:`splat.interpol.PolyMatrix`.  It is usually preferred to create a
continuous list of different polynomials between each input data point using a
:py:class:`splat.interpol.PolyList` oject or using the
:py:func:`splat.interpol.spline` function to directly calculate all the
polynomials and create the list.  It is easier to control the interpolation of
a spline containing many low-degree polynomials than a single high degree
polynomial.


.. autoclass:: splat.interpol.Polynomial
   :members:

.. autoclass:: splat.interpol.PolyMatrix
   :members:

.. autoclass:: splat.interpol.PolyList
   :members:

.. autofunction:: splat.interpol.spline
.. autofunction:: splat.interpol.freqmod

Example using a spline to create a continuous frequency modulation::

    import splat.interpol
    import splat.gen

    pts = [(0.0, 220.0), (1.0, 330.0), (2.0, 110.0)]
    fmod = splat.interpol.freqmod(pts)
    gen = splat.gen.SineGenerator()
    gen.run(fmod.start, fmod.end, fmod.f0, fmod.signal)
    gen.frag.save("freqmod.wav")
