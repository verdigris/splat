``splat``
=========

Fragment objects
----------------

.. autoclass:: splat.Fragment(channels, rate, duration)
   :members:

   .. automethod:: splat.Fragment.mix
   .. automethod:: splat.Fragment.normalize
   .. automethod:: splat.Fragment.amp


Generator objects
-----------------

.. autoclass:: splat.Generator
   :members:
   :private-members:


.. autoclass:: splat.SineGenerator
   :members:

.. autoclass:: splat.OvertonesGenerator
   :members:


.. _sources:

Sound sources
-------------

.. autofunction:: splat.sources.sine
.. autofunction:: splat.sources.overtones


.. _filters:

Filter functions and FilterChain objects
----------------------------------------

A *filter function* takes a :py:class:`splat.Fragment` and a tuple of arguments
with its specific parameters.  It is expected to run on the entirety of the
fragment.  Filter functions can be combined into a series via the
:py:class:`splat.FilterChain` class.

.. autoclass:: splat.FilterChain
   :members:

.. autofunction:: splat.filters.linear_fade
.. autofunction:: splat.filters.dec_envelope
.. autofunction:: splat.filters.reverse
.. autofunction:: splat.filters.reverb


General purpose functions
-------------------------

.. autofunction:: splat.lin2dB
.. autofunction:: splat.dB2lin
