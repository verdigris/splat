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


Sound Sources
-------------

.. autofunction:: geomusic.sources.sine
.. autofunction:: geomusic.sources.overtones
