``splat``
=========

Fragment objects
----------------

.. autoclass:: splat.data.Fragment(channels, rate, duration)
   :members:

   .. automethod:: splat.data.Fragment.mix
   .. automethod:: splat.data.Fragment.normalize
   .. automethod:: splat.data.Fragment.amp
   .. autoattribute:: splat.data.Fragment.sample_rate
   .. autoattribute:: splat.data.Fragment.duration
   .. autoattribute:: splat.data.Fragment.channels


Generator objects
-----------------

.. autoclass:: splat.gen.Generator
   :members:
   :private-members:


.. autoclass:: splat.gen.SineGenerator
   :members:

.. autoclass:: splat.gen.OvertonesGenerator
   :members:


.. _sources:

Sound sources
-------------

.. autofunction:: splat.sources.sine
.. autofunction:: splat.sources.overtones


.. _filters:

Filter functions and FilterChain objects
----------------------------------------

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

Interpolation and Spline objects
--------------------------------

This module contains a set of classes designed to work together and provide
polynomial interpolation functionality.  The principle is to build a continuous
function for a given set of input discrete coordinates.  A
:py:class:`splat.interpol.Polynomial` object represents a polynomial function
with a series of coefficients.  It can be calculated by reducing a matrix
containing some coordinates using :py:class:`splat.interpol.PolyMatrix`.  It is
usually preferred to use a :py:class:`splat.interpol.Spline` object to create a
long function composed of a list of different polynomials between each pair of
points.  It is easier to control the interpolation of a spline containing many
low-degree polynomials tha a single high degree polynomial.

.. autoclass:: splat.interpol.Polynomial
   :members:

.. autoclass:: splat.interpol.PolyMatrix
   :members:

.. autoclass:: splat.interpol.Spline
   :members:


.. _general_purpose_functions:

General purpose functions
-------------------------

.. autofunction:: splat.lin2dB
.. autofunction:: splat.dB2lin
