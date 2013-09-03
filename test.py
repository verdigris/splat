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

# -----------------------------------------------------------------------------
# utilities

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
test_frag.test_name = 'Fragment'

def test_gen_frag():
    gen = splat.gen.SineGenerator()
    return (isinstance(gen.frag, splat.data.Fragment) and
            gen.frag.duration == 0.0)
test_gen_frag.test_name = "Generator Fragment"

def test_sine():
    gen = splat.gen.SineGenerator()
    f = 1000.0
    gen.run(f, 0.0, 1.0)
    n = int(0.1234 * gen.frag.duration * gen.frag.sample_rate)
    s = math.sin(2 * math.pi * f * float(n) / gen.frag.sample_rate)
    check_samples(gen.frag, {n: (s, s)})
    return check_md5(gen.frag, 'ec18389e198ee868d61c9439343a3337')
test_sine.test_name = "SineGenerator"

def test_overtones():
    gen = splat.gen.OvertonesGenerator()
    gen.ot_decexp(1.0)
    f = 1000.0
    gen.run(f, 0.0, 1.0)
    return check_md5(gen.frag, 'ee045e012673ff7ed4ab9bd590b57368')
test_overtones.test_name = "OvertonesGenerator"

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
test_spline.test_name = "Spline"

# -----------------------------------------------------------------------------
# main function

def main(argv):
    tests = []
    for name, value in globals().iteritems():
        if name.startswith('test_'):
            tests.append(value)
    failures = []
    n = len(tests)
    for i, t in enumerate(tests):
        res = t()
        print("{0:03d}/{1:03d} {2:4s} - {3}".format(
                i + 1, n, 'OK' if res is True else 'FAIL', t.test_name))
        if res is False:
            failures.append(t.test_name)

    print("--------------------------------------------------")
    print("Results: {0}/{1} passed".format((n - len(failures)), n))
    if failures:
        print("SOME TESTS FAILED:")
        for f in failures:
            print("    {0}".format(f))
    else:
        print("All good.")
    return (len(failures) == 0)

if __name__ == '__main__':
    ret = main(sys.argv)
    sys.exit(0 if ret is True else 1)
