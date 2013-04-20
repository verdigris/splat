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


General purpose functions
-------------------------

.. autofunction:: splat.lin2dB
.. autofunction:: splat.dB2lin
