# Splat - splat/filters.py
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

import collections
from _splat import dec_envelope, reverse, reverb

class FilterChain(collections.Sequence):

    """Chain of filters to process existing data

    This class is used by the :py:class:`splat.gen.Generator` classes to define
    a chain of filter functions that are run on each newly created
    :py:class:`splat.data.Fragment` instance from a sound source
    (:ref:`sources`).
    """

    def __init__(self, filters=None, *args, **kw):
        """This class implements the :py:class:`collections.Sequence` interface
        to hold a chain of filter functions.  The initial filter functions are
        provided via the ``filters`` list.  The idea is to be able to go
        through them in a specific order over a same audio fragment so the
        output of one filter is the input of the next.
        """
        super(FilterChain, self).__init__(*args, **kw)
        self._filters = list()
        if filters is not None:
            for f in filters:
                if isinstance(f, tuple):
                    self.append(*f)
                else:
                    self.append(f)

    def __getitem__(self, i):
        return self._filters[i]

    def __len__(self):
        return len(self._filters)

    def append(self, filter_func, args=()):
        """Add a filter function to the chain.

        The filter function ``filter_func`` is added to the end of the chain,
        and the provided ``args`` tuple is associated with it to provide
        specific parameters when invoking it by
        :py:meth:`splat.filters.FilterChain.run`.
        """
        if not isinstance(args, tuple):
            raise Exception("Invalid filter arguments, must be a tuple")
        self._filters.append((filter_func, args))

    def run(self, frag):
        """Run all the filter functions on a sound fragment.

        All the filter functions in the chain are run with their associated
        arguments on the :py:class:`splat.data.Fragment` argument ``frag``.
        """
        for f, args in self:
            f(frag, *args)


def linear_fade(frag, duration=0.01):
    """Apply a linear fade-in and fade-out on the fragment.

    For the given ``duration`` in seconds, apply a gain that varies linearly
    from 0.0 to 1.0 and then from 1.0 to 0.0 respectively at the beginning and
    end of the fragment.  This duration is adjusted evenly if the fragment is
    too short.
    """
    fade = min((frag.rate * duration), (len(frag) / 2))
    for i in range(int(fade)):
        l = i / fade
        for j in (i, -i):
            frag[j] = tuple((s * l) for s in frag[j])

def reverb_delays(level=-6.0, length=2.0, dec=0.25, rate=400):
    n2 = int(length * rate)
    n1 = n2 / 80
    r2 = rate
    r1 = r2 / 2
    d = [((float(t) / r1), level - float(t * 2) * dec) for t in range(n1 / 2)]
    d += [((float(t) / r2), level - float(t) * dec) for t in range(n1, n2)]
    return d
