import copy

class Polynomial(object):

    def __init__(self, coefs):
        self._coefs = coefs
        self._deriv = None

    @property
    def pol(self):
        return self._coefs

    @property
    def derivative(self):
        if self._deriv is None:
            coefs = tuple((k * (i + 1)) for i, k in enumerate(self._coefs[1:]))
            self._deriv = Polynomial(coefs)
        return self._deriv

    def value(self, x):
        res = 0
        for p, k in enumerate(self._coefs):
            res += k * (x ** p)
        return res


class PolyMatrix(object):

    def __init__(self, pts):
        self._pts = pts
        self._m = []
        self._r = None
        self._poly = None

    @property
    def m(self):
        if not self._m:
            self._build()
        return self._m

    @property
    def reduced(self):
        if self._r is None:
            self._reduce()
        return self._r

    @property
    def poly(self):
        if self._poly is None:
            self._reduce()
        return self._poly

    def add_pts(self, pts):
        self._pts += pts

    def add_row(self, row):
        self._m.append(row)

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
        self._poly = Polynomial(tuple((l[n] for l in self._r)))


class Spline(object):

    def __init__(self, pts, n=2):
        self._pts = pts
        self._n = 2
        self._pols = []
        self._build()

    @property
    def polynomials(self):
        return self._pols

    def value(self, x):
        for x0, x1, pol in self._pols:
            if (x0 <= x) and (x <= x1):
                return pol.value(x)
        return None

    def _build(self):
        m = PolyMatrix(self._pts[:(self._n + 1)])
        dpol = m.poly.derivative

        for i in range(len(self._pts) - (self._n - 1)):
            seg = self._pts[i:(i + self._n)]
            m = PolyMatrix([])
            m0 = None
            x0, x1 = seg[0][0], seg[1][0]
            if len(seg[1]) > 2:
                d1 = seg[1][2]
                q = self._n + 2
                m0 = [0.0] + [i * (x1 ** (i - 1)) for i in range(1, (q))]+ [d1]
            else:
                q = self._n + 1
            for pt in seg:
                x, y = pt[0], pt[1]
                m.add_row([(x ** p) for p in range(q)] + [y])
            m.add_row([0.0] + [i * (x0 ** (i - 1)) for i in range(1, (q))]
                      + [dpol.value(x0)])
            if m0 is not None:
                m.add_row(m0)
            dpol = m.poly.derivative
            self._pols.append((x0, x1, m.poly))

        if self._n > 2:
            self._pols.append((x1, pts[-1][0], m.poly))
