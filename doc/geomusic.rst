``geomusic``
============

Fragment objects
----------------

.. autoclass:: geomusic.Fragment(channels, rate, duration)
   :members:

   .. automethod:: geomusic.Fragment.mix
   .. automethod:: geomusic.Fragment.normalize


Generator objects
-----------------

.. autoclass:: geomusic.Generator
   :members:
   :private-members:


.. autoclass:: geomusic.SineGenerator
   :members:

.. autoclass:: geomusic.OvertonesGenerator
   :members:


.. _sources:

Sound sources
-------------

.. autofunction:: geomusic.sources.sine
.. autofunction:: geomusic.sources.overtones


.. _filters:

Filter functions and FilterChain objects
----------------------------------------

A *filter function* takes a :py:class:`geomusic.Fragment` and a tuple of
arguments with its specific parameters.  It is expected to run on the entirety
of the fragment.  Filter functions can be combined into a series via the
:py:class:`geomusic.FilterChain` class.

.. autoclass:: geomusic.FilterChain
   :members:

.. autofunction:: geomusic.filters.linear_fade
.. autofunction:: geomusic.filters.dec_envelope


General purpose functions
-------------------------

.. autofunction:: geomusic.lin2dB
.. autofunction:: geomusic.dB2lin
