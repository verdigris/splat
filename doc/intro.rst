What is Splat?
==============

Introduction
------------

The basic idea is to apply mathematical concepts to musical composition as well
as sound synthesis.  In practice, Splat lets you create just about any sound
you can imagine and code in software.  A Splat composition, here also called a
*splat*, takes the form of a Python script, and can generate as well as import
existing audio data and run various processing operations.

.. figure:: image/splat-foot.png
   :align: center

   Splat is also the sound produced by a stomping bare foot.

This does not involve any real-time sound generation: everything runs at its
own pace to ultimately create an audio file which can then be played.  This is
not too good for live shows, but it makes it possible to run very complicated
operations that may take much longer to create the audio material than it does
to play it.  The power of the computer running the composition will not change
the actual audio result; only the time to generate it will be shorter with a
faster machine.

It is usually easier to prototype new tools (i.e. generators, filters ...) in
pure Python, and when some code is known to be working it can be implemented
again in C for improved performance (see ``_splat.c`` in source code
distribution).

Of course, a major advantage of Splat is that it can be used in conjunction
with all the other Python packages already available, for example with the
Python Imaging Library to create both sound and images at the same time or with
``numpy`` for advanced calculations.  Splat does not depend on any third-party
code, only on the standard Python library.

.. _beep:

Beep
----

As a first concrete example, here's a very small splat which creates a 440 Hz
sound for 1 second, with fade-in and fade-out, and saves it in a file:

.. code-block:: python

    import splat.gen
    import splat.filters

    gen = splat.gen.SineGenerator()
    gen.filters = [splat.filters.linear_fade]
    gen.run(0.0, 1.0, 440.0)
    gen.frag.save("440-Hz.wav")

Skipping the import statements, this first creates the sound generator
object: a :py:class:`splat.gen.SineGenerator` to generate sine waves.
Then it sets a :py:func:`splat.filters.linear_fade` filter which will
run on the generator output to progressively ramp up and down (fade in
and out) the volume at the beginning and end of the generated audio.
Then the generator is being run, starting at time 0 for 1 second and
at a frequency of 440 Hz (standard concert pitch for the A note...).
Generators have a piece of memory associated with them to keep the
data they produce, so that's where the sine wave went.  Finally the
audio is saved in a file which you can play any time and share with
your friends.


Principles
----------

.. rubric:: Fragments

Splat ultimately produces some sound fragments as
:py:class:`splat.data.Fragment` objects.  They contain sound samples in
floating point numbers in a fixed number of channels.  Fragments can be created
by sound generators or from audio files, mixed together, copied and edited and
finally saved into audio files (see :ref:`audio_files`).  They include basic
transformation methods such as :py:meth:`splat.data.Fragment.normalize` to
normalize the amplitude and :py:meth:`splat.data.Fragment.mix` to mix one
fragment, or a part of it, into another.  The sound data can also be edited
directly as a mutable sequence.

.. rubric:: Sources

Fragments can also be passed to :ref:`sources`, which are essentially functions
that will fill the fragment with some sound content.  Each source takes
arbitrary parameters.  For example :py:func:`splat.sources.sine` takes levels,
frequency and phase.  Other sources like :py:func:`splat.sources.square` take a
ratio parameter to specify the relative lengths of the high and low states
within each period.

.. rubric:: Filters

There are also some :ref:`filters`, which alter existing fragment data.
Filters are functions which take a fragment with some parameters and modify it
in place.  Unlike sources which typically ignore any data currently in a
fragment, filters will take this data as their input, then replace it with
modified data in the same fragment.  Filters can be combined together into an
ordered series using a :py:class:`splat.filters.FilterChain` to easily manage a
chain of effects.  They will all run in turn over the same fragment.

.. rubric:: Generators

Then :py:class:`splat.gen.Generator` objects can be used to manage sources,
filters and fragments together and produce some audio data.  They hold a main
fragment and can run some arbitrary code to call sources with various
arguments, run any extra filters, and mix each newly generated sound fragment
into the main fragment.  Generators are designed with an abstract ``run``
method with the aim to make them interchangeable, as much as possible.  This is
to enable a same piece of code to be run again with different generator
implementations for different souding results.  They are a little bit like
usual musical instruments which all have the same definition of notes but
produce different sounds.

.. rubric:: Scales

In order to use musical notes, there are some :ref:`scales` which will
calculate the frequency of a given note based on its name and scale
implementation.  Generators and sources typically work with frequencies in
Hertz, so :ref:`scale_objects` are very useful when generating some musical
notes.  They can be used to create conventional scales (tempered, diatonic) but
also to experiment with new ideas such as arbitrary ratios, microtones or
dynamic scales.

.. rubric:: Sequencers

:ref:`sequencer_objects` can be used to control when notes and sounds get
generated and mixed in time.  The basic :py:class:`splat.seq.PatternSequencer`
can create simple repetitive patterns, but it also serves as a basis for more
creative implementations with irrational rhythms, dynamic tempo, random number
of beats in each bar...

.. rubric:: Interpolation

More on the mathematical side of things, there are also some
:ref:`interpolation` tools.  The main feature is the
:py:func:`splat.interpol.spline` function which creates a polynomial spline
curve as a :py:class:`splat.interpol.PolyList` object with a given list of
points and optional slopes to constrain it.  This can then be typically used to
modulate sources but in the end it just returns interpolated numbers.  There
are endless ways to use it within Splat or any other application.

.. rubric:: Signals

Finally, there is a recurring Splat concept called **signals**.  They refer to
objects which for a given value return an other one, and can be either of the
following types:

Python *callable* (function, object method...)
  They must accept a single floating point argument and return one.  Any
  arbitrary code can be used to create the signal behaviour.
constant *floating point value*
  The same value is always used.  This is much faster than a function returning
  always the same value.
:py:class:`splat.data.Fragment` mono object
  A Fragment with a single channel can be used with the input value being a
  time in seconds to look up a sample value as the signal output.  This is
  especially useful when the signal behaviour can't be implemented by a
  function.  It may also be used as a data cache to avoid repeated expensive
  computations.
:py:attr:`splat.interpol.PolyList.signal` object
  This property provides a special object type with an optimised C
  implementation to use splines as signals.  The end result is the same as
  using the :py:meth:`splat.interpol.PolyList.value` method except that the
  optimised signal is a lot faster.

Splat signals can be used in many places to provide modulations and other kinds
of dynamic behaviour.  For example, sources can be called with a spline signal
as their amplitude parameter in order to generate an arbitrary envelope (or
amplitude modulation, or tremolo when it's periodic).

This is handled under the hood in the **C extension module** ``_splat`` by
looking at the type of signal parameters.  Functions usually automatically
detect when all their arguments are plain floating point values and run a
faster optimised implementation.  Below is an example with a phase modulation
to produce a vibrato effect using a lambda function as a signal (sounds more
like a siren)::

    import math
    import splat.gen

    gen = splat.gen.TriangleGenerator()
    vibrato = lambda t: 0.01 * math.sin(15.0 * t)
    gen.run(0.0, 1.0, 880.0, phase=vibrato)
    gen.frag.save('vibrato.wav')

It's also possible to create a :py:class:`splat.Signal` object to use the
functionality at the Python level.  This basically takes a signal argument
(float, callable or Fragment) and creates a sequence object which can then be
indexed.

Typical Splat
-------------

Your Splat program will most likely do some or all of the following things:

#. Create a :py:class:`splat.gen.Genarator` which contains a fragment and a
   sound source.
#. Optionally, create filters and build a filter chain into the generator.
#. Run the generator to generate some sound in its internal fragment.
#. Create other generators and sequencers to generate more sound fragments.
#. Import existing audio files into fragments (no generators involved here).
#. Run any filters or edit the fragments.
#. Finally mix all the fragments together and save the result into a file.

.. raw:: latex

   \newpage

.. rubric:: And now, in stereo

Below is a more advanced example, also this `example.py on Github <https://github.com/verdigris/splat/blob/master/example.py>`_:

.. literalinclude:: ../example.py
   :language: python
