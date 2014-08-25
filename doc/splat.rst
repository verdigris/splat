Fragment objects
================

.. automodule:: splat.data

.. autoclass:: splat.data.Fragment(channels, rate, duration=0.0, length=0)
   :members:

   .. automethod:: splat.data.Fragment.mix
   .. automethod:: splat.data.Fragment.import_bytes
   .. automethod:: splat.data.Fragment.export_bytes
   .. automethod:: splat.data.Fragment.normalize
   .. automethod:: splat.data.Fragment.amp
   .. automethod:: splat.data.Fragment.offset
   .. automethod:: splat.data.Fragment.resize
   .. autoattribute:: splat.data.Fragment.rate
   .. autoattribute:: splat.data.Fragment.duration
   .. autoattribute:: splat.data.Fragment.channels


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


.. _interpolation:

Interpolation and splines
=========================

.. automodule:: splat.interpol

The :py:mod:`splat.interpol` module contains a set of classes designed to work
together and provide polynomial interpolation functionality.  The principle is
to build a continuous function for a given set of input discrete coordinates.
A :py:class:`splat.interpol.Polynomial` object represents a polynomial function
with a series of coefficients.  It can be calculated by reducing a matrix
containing some coordinates using :py:class:`splat.interpol.PolyMatrix`.  It is
usually preferred to use :py:func:`splat.interpol.spline` to create a long
function composed of a list of different polynomials between each pair of
points.  It is easier to control the interpolation of a spline containing many
low-degree polynomials than a single high degree polynomial.

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
    gen.run(fmod.start, fmod.end, fmod.f0, fmod.value)
    gen.frag.save("freqmod.wav")


.. _utilities:

Utilities
=========


Conversion functions
--------------------

.. autofunction:: splat.lin2dB
.. autofunction:: splat.dB2lin


Signal objects
--------------

.. autoclass:: splat.Signal
   :members:


.. _web:

On the web
==========

There are some concrete examples as `Gists on Github
<https://gist.github.com/gctucker>`_.

Some experimental splats can also be heard on `SoundCloud
<http://soundcloud.com/verdigris-mu>`_.

You can find some splats as well as other music, software and electronics
things on `verdigris.mu <http://verdigris.mu>`_.

Please make your creations or reactions be heard by sending them to
`info@verdigris.mu <mailto:info@verdigris.mu>`_.
