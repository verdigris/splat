# Splat - splat/interpol.py
#
# Copyright (C) 2013, 2014 Guillaume Tucker <guillaume@mangoz.org>
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

import copy
import _splat

class Polynomial(object):

    """Polynomial object

    A series of coefficients are used to create a polynomial function.  The
    index of each coefficient corresponds to the power of the input variable,
    starting with 0 and increasing by 1 each time.
    """

    def __init__(self, coefs):
        """The ``coefs`` are a sequence of floats with the coefficients of the
        polynomial function.  The degree of the polynomial is equal to the
        length of ``coefs``.
        """
        self._coefs = tuple(coefs)

    def __repr__(self):
        return repr(self.coefs)

    @property
    def coefs(self):
        """A tuple with the polynomial coefficients."""
        return self._coefs

    def derivative(self):
        """Create a :py:class:`splat.interpol.Polynomial` object with the
        derivative of this polynomial."""
        coefs = tuple((k * (i + 1)) for i, k in enumerate(self.coefs[1:]))
        return Polynomial(coefs)

    def integral(self, y0=0.0):
        """Create a :py:class:`splat.interpol.Polynomial` object with an
        integral of this polynomial.  The ``y0`` value is the arbitrary
        initial value."""
        coefs = (y0,) + tuple((k / (i + 1)) for i, k in enumerate(self.coefs))
        return Polynomial(coefs)

    def value(self, x, dB2lin=False):
        """Return the value of the polynomial for the given ``x`` input
        value."""
        return _splat.poly_value(self.coefs, x, dB2lin)


class PolyMatrix(object):

    """Matrix to calculate polynomial coefficients

    For a given list of input points, polynomial coefficients are calculated so
    that the function passes through all of these points.  The order of the
    polynomial function is equal to the number of input points minus one.

    The calculation is achieved by building a matrix from the input coordinates
    and then reducing it with linear operations.  It uses what is now known as
    the *Gauss-Jordan elimination*.  A :py:class:`splat.interpol.Polynomial`
    object is created as a result.  It is also possible to append user-defined
    rows to the matrix in addition to the ones automatically created from the
    input points.
    """

    def __init__(self, pts):
        """The ``pts`` are a list of 2-tuples with ``(x, y)`` coordinates."""
        self._pts = pts
        self._m = []
        self._r = None
        self._poly = None

    @property
    def m(self):
        """Get the matrix numbers as ``n`` lists, each of them containing a row
        of ``(n + 1)`` elements."""
        if not self._m:
            self._build()
        return self._m

    @property
    def reduced(self):
        """Get the reduced matrix, which contains the identity matrix and the
        polynomial coefficients at the end of each row."""
        if self._r is None:
            self._reduce()
        return self._r

    @property
    def poly(self):
        """Get the :py:class:`splat.interpol.Polynomial` object derived from
        the reduced matrix."""
        if self._poly is None:
            n = len(self.reduced)
            self._poly = Polynomial(tuple((l[n] for l in self.reduced)))
        return self._poly

    def add_row(self, row):
        """Append a row to the matrix.  It must be a list of floats with the
        same length as the others (``n + 1``).  This resets the reduced matrix
        and the polynomial."""
        self.m.append(row)
        self._r = None
        self._poly = None

    def _build(self):
        n = len(self._pts)
        for pt in self._pts:
            x, y = pt[0], pt[1]
            self._m.append([(x ** i) for i in range(n)] + [y])

    def _reduce(self):
        self._r = copy.deepcopy(self.m)
        n = len(self._r)
        for i in range(n):
            a = self._r[i][i]
            for k in range(n + 1):
                self._r[i][k] /= a
            for j in range(n):
                if i != j:
                    r = self._r[j][i]
                    for k in range(i, n + 1):
                        self._r[j][k] -= (self._r[i][k] * r)


class PolyList(object):

    """List of 3-tuples containing a range where values are defined and a
    :py:class:`splat.interpol.Polynomial` object to compute the values.  Each
    item in the list corresponds to a segment between 2 input points.  The
    segments do not need to be contiguous, although looking for a value outside
    of the defined ranges will cause an error."""

    def __init__(self, pols, scale=1.0, dB2lin=False):
        """The list of polynomial is built using the sequence of 3-tuples
        ``pols``.  The ``scale`` argument is a convenient way to use the same
        input points and independently adjust the scale of the values, which
        get multiplied by this scale coefficient (1.0 by default)."""
        self._pols = pols
        self._pols_coefs = list((x0, x1, pol.coefs) for x0, x1, pol in pols)
        self._signal = _splat.Spline(self._pols_coefs, scale, dB2lin)
        self._scale = scale
        self._dB2lin = dB2lin

    def __getitem__(self, i):
        return self._pols.__getitem__(i)

    @property
    def pols(self):
        return copy.copy(self._pols)

    @property
    def scale(self):
        """Scale factor for all returned values."""
        return self._scale

    @scale.setter
    def scale(self, value):
        self._scale = value

    @property
    def start(self):
        """First ``x`` value from the list of points.  The result is undefined
        for values smaller than this."""
        return self._pols[0][0]

    @property
    def end(self):
        """Last ``x`` value from the list of points.  The result is undefined
        for values greater than this."""
        return self._pols[-1][1]

    @property
    def signal(self):
        """Object with optimised signal implementation."""
        return self._signal

    def points(self):
        xlist = list(p[0] for p in self)
        xlist.append(self[-1][1])
        return list((x, self.value(x)) for x in xlist)

    def integral(self, y0=0.0):
        """Get a new PolyList object with the integral of this one.

        For each polynomial, its integral is calculated based on its
        coefficients and the last value of the previous polynomial to maintain
        continuity.  The initial value for the first polynomial can be set with
        ``y0``.  The scale is preserved from the current object.
        """
        pols = []
        for x0, x1, p in self:
            p2 = p.integral(0.0)
            y1 = p2.value(x0)
            p2._coefs = (y0 - y1,) + p2._coefs[1:]
            pols.append((x0, x1, p2))
            y0 = p2.value(x1)
        return PolyList(pols, self.scale)

    def derivative(self):
        """Get a new PolyList object with the derivative of this one.

        For each polynomial, its derivative is calculated based on its
        coefficients.  The scale is preserved from the current object.
        """
        return PolyList(list((x0, x1, p.derivative()) for x0, x1, p in self),
                        self.scale)

    def value(self, x):
        """Return the spline value for a given ``x`` input value, or ``None``
        if undefined."""
        y = _splat.spline_value(self._pols_coefs, x, self._dB2lin)
        if y is None:
            return None
        return y * self.scale

    def slices(self, y0, xmin=None, xmax=None, xstep=0.001):
        """Get slices of the a given ``y0`` value.

        Return a list of 2-tuples with ``x`` ranges for which the spline is
        greater or equal to the given ``y0`` value.  Optional arguments
        ``xmin`` and ``xmax`` can be used to restrict the interval instead of
        the whole spline definition range.  The ``xstep`` argument defines the
        granularity used in generating the range coordinates.  Smaller values
        make more accurate results but also increase the processing time.

        This can for example be used to create arbitrary random distributions
        by adding evenly distributed random events within slices while
        iterating over a ``y`` range.  It can also be used to calculate an
        approximation of the area of the spline by adding the area of
        rectangular slices between a set of ``y`` values.
        """
        if xmin is None:
            xmin = self.start
        if xmax is None:
            xmax = self.end
        y1 = self.value(xmin)
        if y1 is None:
            return None
        slices = []
        x0 = None
        s1 = bool(y1 < y0)
        for i in range(1, int((xmax - xmin) / xstep)):
            x = xmin + float(i) * xstep
            y = self.value(x)
            if y is None:
                break
            s = bool(y < y0)
            if s1 and not s:
                x0 = x
            elif (x is not None) and (s and not s1):
                slices.append((x0, x))
            y1 = y
            s1 = s
        return slices


def spline(pts, scale=1.0, n=2, dB2lin=False):
    """Build a spline from a series of points and slope coordinates.

    For a given list of ``(x, y)`` or ``(x, y, d)`` coordinates, where ``d`` is
    a slope value, create and return a :py:class:`splat.interpol.PolyList` list
    of polynomials of degree ``n`` or ``n + 1`` if the slope is specified.
    This can then be used as a continuous ``f(x) = y`` function.  Each
    polynomial is used to interpolate between 2 input points, except at the end
    where 3 points may share a same polynomial.  They are calculated so that
    they go through all of the given ``(x, y)`` 2-tuple coordinates.  The slope
    (or derivative value) at each point is either determined by the previous
    interpolation polynomial or constrained in the last value in 3-tuple input
    points.  The calculation is done using a
    :py:class:`splat.interpol.PolyMatrix` object.
    """
    pts = list(tuple(float(x) for x in pt) for pt in sorted(pts))
    m = PolyMatrix(pts[:(n + 1)])
    dpol = m.poly.derivative()
    pols = []

    for i in range(len(pts) - (n - 1)):
        seg = pts[i:(i + n)]
        m = PolyMatrix([])
        m0 = None
        x0, x1 = seg[0][0], seg[1][0]
        if len(seg[1]) > 2:
            d1 = seg[1][2]
            q = n + 2
            m0 = [0.0] + [i * (x1 ** (i - 1)) for i in range(1, (q))]+ [d1]
        else:
            q = n + 1
        for pt in seg:
            x, y = pt[0], pt[1]
            m.add_row([(x ** p) for p in range(q)] + [y])
        m.add_row([0.0] + [i * (x0 ** (i - 1)) for i in range(1, (q))]
                  + [dpol.value(x0)])
        if m0 is not None:
            m.add_row(m0)
        dpol = m.poly.derivative()
        pols.append((x0, x1, m.poly))

    if n > 2:
        pols.append((x1, pts[-1][0], m.poly))

    return PolyList(pols, scale, dB2lin)


def freqmod(pts):
    """Create a spline for frequency modulation.

    Return a :py:class:`splat.interpol.PolyList` object which can be used as a
    phase modulation to create the desired frequency curve defined by the
    ``pts`` when passed to a sound source or generator.  The returned object
    also has an extra attribute ``f0`` with the initial frequency."""
    f0 = pts[0][1]
    norm_pts = list((pt[0], (pt[1] - f0)) + pt[2:] for pt in pts)
    freqmod = spline(norm_pts, (1 / f0)).integral()
    freqmod.f0 = f0
    return freqmod
