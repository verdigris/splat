.. _splat_fast:

Fast mode & SIMD
================

Splat provides a fast mode which uses hardware specific features to perform
faster.  In standard mode, the sample precision is 64-bit.  Fast mode uses
32-bit per sample which allows even faster computation at the cost of accuracy.
The standard tests may not pass in fast mode due to inaccuracy, so there is a
separate set of tests to measure the amplitude of errors compared to standard
mode.  Some functions such as sine wave generation on ARM NEON will run up to
10 times faster in fast mode and the error amplitude is kept under -40dB in all
cases.  The main idea behind fast mode is to allow fast prototyping and shorter
development cycles with a quality that is still good enough to have a clear
idea of the end result.  It may also enable the use of Splat on embedded
systems with limited resources.  Standard mode should still be used for final
runs with high quality audio.

Fast mode is enabled at compilation time using the ``SPLAT_FAST`` macro and by
enabling compiler flags corresponding to the CPU hardware.  Currently, some
functions have been implemented using ARM NEON and x86 SSE intrinsics.  There
are convenience commands ``build-neon-v7`` and ``build-sse`` to build Splat
with these flags turned on.

Full testing of fast mode is automated with the ``rebuild-fast-test.sh`` script
which will rebuild in standard mode, run standard tests, run fast tests in
standard mode to generate reference data, rebuild in fast mode (NEON or SSE),
run the fast tests again, compare the results with standard mode and print the
peak error amplitudes in dB.
