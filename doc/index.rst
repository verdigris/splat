.. Splat documentation master file, created by
   sphinx-quickstart on Tue Mar  5 22:42:04 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. only:: html

   Splat!
   ======

   **Splat is a programme** to generate some audio data which you may call
   music.  It's written in Python to make it easy to use, and all the crucial
   processing parts are implemented in C for fast number crunshing.  It's
   distributed under the terms of the **GNU LGPL v3** so you remain a free
   individual when you use it.

   **Splat is not a programme** to generate **live** audio in real time, at
   least not at the moment.  So it's better suited to the studio than to the
   stage.

   Here's a very small **Splat** example which creates a beep 440Hz sound:

   .. code-block:: python

       import splat

       gen = splat.SineGenerator(splat.Fragment(), [splat.filters.linear_fade])
       gen.run(440.0, 0.0, 1.0)
       gen.frag.save("A440.wav")


Splat manual
============

.. toctree::
   :maxdepth: 2

   intro
   audio_formats
   splat
   fast
   web
