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
not good for live shows, but it makes it possible to run very complicated
operations that may take much longer to create the audio material than it does
to play it.  The power of the computer running the composition will not change
the actual audio result; only the time to generate it will be shorter with a
faster machine.

It is usually easier to prototype new tools (i.e. generators, filters ...) in
pure Python, and when some code is known to be working it can be implemented
again in C for improved performance (see ``_splat.c`` in source code
distribution).

Note: Only ``WAV`` audio file format is currently supported for file imports
and exports.

Principles
----------

.. rubric:: Fragments

Splat ultimately produces some sound fragments, as
:py:class:`splat.data.Fragment` objects.  They contain sound samples in
floating point numbers for a fixed number of channels.  Fragments can be mixed
together, copied and edited and finally saved into audio files.  They also
include basic transformation methods such as
:py:meth:`splat.data.Fragment.normalize` to normalize the amplitude and
:py:meth:`splat.data.Fragment.mix` to mix one fragment, or a part of it, into
another.

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
in place.  Unlike sources which typically overwrite any data currently in a
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
to enable a same composition to be run again with different generator
implementations for different results.  They are a little bit like musical
instruments, which typically all have the same definition of notes but produce
different sounds.

.. rubric:: Scales

Then to use actual musical notes, there are some **scales** which will
calculate the frequency of a given note based on its name and scale
implementation.  Generators and sources typically work with frequencies in
Hertz, so scale objects are very useful when writing some musical notes.  They
can be used to create conventional scales but also to experiment with new
ideas using arbitrary ratios.

.. rubric:: Interpolation

Other than that, there are some :ref:`interpolation` tools.  The main feature
is the :py:class:`splat.interpol.Spline` class which can create a polynomial
spline curve with a given list of points and optional slopes to constrain it.
This can then be typically used to modulate sources but in the end it just
returns interpolated numbers.  There are endless ways to use it within Splat.

.. rubric:: Signals

Finally, there is a recurring Splat concept called **signals**.  They refer to
objects which can be either of the following types:

* a constant *floating point value*,
* a Python *callable* (function, object method...),
* a :py:class:`splat.data.Fragment` object.

Splat signals can be used in many places to provide modulations and other
dynamic behaviour.  For example, sources can be called with a
:py:meth:`splat.interpol.Spline.value` method as their amplitude parameter in
order to generate an arbitrary envelope (or amplitude modulation).

This is handled under the hood in the **C extension module** ``_splat`` by
looking at the type of signal parameters.  Functions usually automatically
detect when all their arguments are plain floating point values and run a
faster optimised implementation.  Below is an example with a phase modulation
to produce a vibrato effect (sounds more like a siren)::

    import math
    import splat.gen

    gen = splat.gen.TriangleGenerator()
    vibrato = lambda t: 0.01 * math.sin(15.0 * t)
    gen.run(0.0, 1.0, 880.0, phase=vibrato)
    gen.frag.normalize()
    gen.frag.save("vibrato.wav")

It's also possible to create a :py:class:`splat.Signal` object to use the
functionality at the Python level.  This basically takes a signal argument
(float, callable or Fragment) and creates a sequence object which can then be
indexed.

Typical Splat programme
-----------------------

Your Splat programme will most likely do some or all of the following things:

#. Create a :py:class:`splat.gen.Genarator` which contains a fragment and a
   sound source.
#. Optionally, create filters and build a filter chain into the generator.
#. Run the generator to generate some sound in its internal fragment.
#. Create other generators to create other sound fragments.
#. Import existing audio files into fragments (no generators involved here).
#. Run any filters or edit the fragments.
#. Finally mix all the fragments together and save the result into a file.

.. rubric:: And now, an example in stereo (see ``example.py``)

.. literalinclude:: ../example.py
   :language: python
