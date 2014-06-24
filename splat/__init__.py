# Splat - splat/__init__.py
#
# Copyright (C) 2012, 2013, 2014 Guillaume Tucker <guillaume@mangoz.org>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import _splat
from _splat import lin2dB, dB2lin
from _splat import SAMPLE_INT, SAMPLE_FLOAT
from _splat import NATIVE_SAMPLE_TYPE, NATIVE_SAMPLE_WIDTH

__all__ = ['gen', 'data', 'filters', 'sources', 'scales', 'interpol']

VERSION = (1, 2)
VERSION_STR = '.'.join([str(v) for v in VERSION])
BUILD = 8

__version__ = VERSION_STR

def check_version(ver):
    if ver != VERSION:
        raise Exception("Version mismatch: {0}, required: {1}".format(
                VERSION, ver))

audio_formats = [(SAMPLE_FLOAT, 64), (SAMPLE_FLOAT, 32),
                 (SAMPLE_INT, 16), (SAMPLE_INT, 8)]

class Signal(_splat.Signal):

    """A general purpose signal.

    This class provides a thin Python wrapper around the signal functionality
    implemented in the C ``_splat`` extension.  It takes a signal and provides
    a sequence object which can be indexed and iterated.  Signals are limited
    in time and have a fixed duration.
    """

    def __init__(self, frag, sig_obj, duration=None):
        """The ``frag`` argument is a :py:class:`splat.data.Fragment` object
        with which the Signal will be compatible.  This is used to determine
        the sample rate and the duration of the signal.  Then ``sig_obj`` is
        the signal, which can be either a floating point value, a callable or a
        Fragment object.  The ``duration`` can be optionally used to provide an
        arbitrary Signal duration in seconds, different from the ``frag``
        duration which is used by default.

        When using a **callable** signal, it must accept a single floating
        point argument and return a floating point value.  As a minimal
        example, it may be a lambda function (that's a tremolo)::

          sig = splat.Signal(frag, lambda x: 0.8 + 0.2 * sin(x * 1000.0))

        The input value is often a time value, but it may be any dimension
        depending on the usage of the signal.  This is compatible with the
        :py:meth:`splat.interpol.Spline.value` method, so a Spline object can
        be used in conjunction with a Signal object.
        """
        args = (frag, sig_obj)
        if duration is not None:
            args += (duration,)
        super(Signal, self).__init__(*args)
