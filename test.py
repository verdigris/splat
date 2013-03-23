import sys
import md5
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
import splat

# -----------------------------------------------------------------------------
# utilities

def check_md5(frag, hexdigest):
    md5sum = md5.new(frag.as_bytes(2))
    if md5sum.hexdigest() != hexdigest:
        print("MD5 mismatch: {0} {1}".format(md5sum.hexdigest(), hexdigest))
        return False
    else:
        return True

# -----------------------------------------------------------------------------
# test functions

def test_frag():
    frag = splat.Fragment(duration=1.0)
    return check_md5(frag, 'fe384f668da282694c29a84ebd33481d')
test_frag.test_name = 'Fragment'

def test_sine():
    gen = splat.SineGenerator(splat.Fragment(duration=1.0))
    gen.run(1000, 0.0, 1.0)
    return check_md5(gen.frag, 'ec18389e198ee868d61c9439343a3337')
test_sine.test_name = "SineGenerator"

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
