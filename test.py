import sys
import md5
import math
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
import splat.data
import splat.gen
import splat.interpol

g_id = 1

# -----------------------------------------------------------------------------
# utilities

def set_id(test, name):
    global g_id
    setattr(test, 'test_name', name)
    setattr(test, 'test_id', g_id)
    g_id += 1

def check_md5(frag, hexdigest):
    md5sum = md5.new(frag.as_bytes(2))
    if md5sum.hexdigest() != hexdigest:
        print("MD5 mismatch: {0} {1}".format(md5sum.hexdigest(), hexdigest))
        return False
    else:
        return True

def check_samples(frag, samples):
    for n, s in samples.iteritems():
        for c, d in zip(frag[n], s):
            if not floatcmp(c, d):
                print("Sample mismatch {0}: {1} {2}".format(n, frag[n], s))
                return False
    return True

def floatcmp(f1, f2):
    return abs(f1 - f2) < (10 ** -6)

# -----------------------------------------------------------------------------
# test functions

def test_frag():
    frag = splat.data.Fragment(duration=1.0)
    return (check_samples(frag, {int(len(frag) / 2): (0.0, 0.0)}) and
            check_md5(frag, 'fe384f668da282694c29a84ebd33481d'))
set_id(test_frag, 'Fragment')

def test_gen_frag():
    gen = splat.gen.SineGenerator()
    return (isinstance(gen.frag, splat.data.Fragment) and
            gen.frag.duration == 0.0)
set_id(test_gen_frag, "Generator Fragment")

def test_sine():
    gen = splat.gen.SineGenerator()
    f = 1000.0
    gen.run(0.0, 1.0, f)
    n = int(0.1234 * gen.frag.duration * gen.frag.sample_rate)
    s = math.sin(2 * math.pi * f * float(n) / gen.frag.sample_rate)
    return (check_samples(gen.frag, {n: (s, s)}) and
            check_md5(gen.frag, 'ec18389e198ee868d61c9439343a3337'))
set_id(test_sine, "SineGenerator")

def test_square():
    gen = splat.gen.SquareGenerator()
    f = 1000.0
    gen.run(0.0, 1.0, f)
    nf = gen.frag.sample_rate / f
    samples = {int(nf * 0.1): (1.0, 1.0), int(nf * 0.9): (-1.0, -1.0)}
    return (check_samples(gen.frag, samples) and
            check_md5(gen.frag, '0ca047e998f512280800012b05107c63'))
set_id(test_square, "SquareGenerator")

def test_triangle():
    gen = splat.gen.TriangleGenerator()
    f = 1000.0
    ratio = 0.567
    gen.run(0.0, 1.0, f, 0.0, ratio, levels=(0.0, 0.0))
    nf = gen.frag.sample_rate / f
    x1 = 0.25
    t1 = int(nf * ratio * x1)
    s1 = (t1 * 2.0 / (ratio * nf)) - 1.0
    x2 = 0.75
    ratio2 = 1 - ratio
    t2 = int(nf * (ratio + (ratio2 * x2)))
    a2 = -2.0 / (ratio2 * nf)
    b2 = 1.0 - (a2 * ratio * nf)
    s2 = (t2 * a2) + b2
    samples = {t1: (s1, s1), t2: (s2, s2)}
    return (check_samples(gen.frag, samples) and
            check_md5(gen.frag, 'b6d9eb000b328134cd500173b24f1c88'))
set_id(test_triangle, "TriangleGenerator")

def test_overtones():
    gen = splat.gen.OvertonesGenerator()
    gen.ot_decexp(1.0)
    f = 1000.0
    gen.run(0.0, 1.0, f)
    return check_md5(gen.frag, 'ee045e012673ff7ed4ab9bd590b57368')
set_id(test_overtones, "OvertonesGenerator")

def test_spline():
    pts = [(1.23, 4.56), (4.32, 2.54, 1.25), (5.458, -4.247)]
    s = splat.interpol.Spline(pts)
    for p in pts:
        x, y = p[0], p[1]
        y1 = s.value(x)
        if not floatcmp(y, y1):
            print("Spline error: s({0}) = {1} instead of {2}".format(x, y1, y))
            return False
    return True
set_id(test_spline, "Spline")

# -----------------------------------------------------------------------------
# main function

def main(argv):
    tests = []
    for name, value in globals().iteritems():
        if name.startswith('test_'):
            tests.append(value)
    failures = []
    n = len(tests)
    for t in sorted(tests, cmp=lambda a, b: cmp(a.test_id, b.test_id)):
        res = t()
        print("{:03d}/{:03d} {:4s} - {}".format(
                t.test_id, n, 'OK' if res is True else 'FAIL', t.test_name))
        if res is False:
            failures.append(t)

    print("--------------------------------------------------")
    print("Results: {0}/{1} passed".format((n - len(failures)), n))
    if failures:
        print("SOME TESTS FAILED:")
        for f in failures:
            print("    {:03d} {}".format(f.test_id, f.test_name))
    else:
        print("All good.")
    return (len(failures) == 0)

if __name__ == '__main__':
    ret = main(sys.argv)
    sys.exit(0 if ret is True else 1)
