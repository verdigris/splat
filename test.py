import sys
import md5
import math
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
import splat
import splat.data
import splat.gen
import splat.sources
import splat.interpol

splat.check_version((1, 1))

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

def check_multiple_md5(frags, hexdigest):
    for frag in frags:
        if check_md5(frag, hexdigest) is False:
            return False
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

def test_frag_offset():
    offset = 1234.5
    frag = splat.data.Fragment(duration=1.0)
    frag.offset(offset)
    frag_twice = splat.data.Fragment(duration=1.0)
    frag_twice.offset(offset)
    offset2 = -57.25
    offset_check = offset + offset2
    frag_twice.offset(offset2)
    frag_signal = splat.data.Fragment(duration=1.0)
    frag_signal.offset(lambda x: offset)
    n = int(0.5678 * frag.duration * frag.sample_rate)
    return (check_samples(frag, { n: (offset, offset) }) and
            check_samples(frag_signal, { n: (offset, offset) }) and
            check_samples(frag_twice, { n: (offset_check, offset_check) }) and
            check_multiple_md5([frag, frag_signal],
                               '248070c79f99014cf800d05ea81e0679'))
set_id(test_frag_offset, "Fragment offset")

def test_gen_frag():
    gen = splat.gen.SineGenerator()
    return (isinstance(gen.frag, splat.data.Fragment) and
            gen.frag.duration == 0.0)
set_id(test_gen_frag, "Generator Fragment")

def test_signal():
    duration = 0.0123
    x1 = int(duration * 0.234)
    x2 = int(duration * 0.789)
    frag = splat.data.Fragment(duration=duration)
    float_value = 1.234
    float_tuple = (float_value,)
    sig_float = splat.Signal(frag, float_value)
    if len(sig_float) != len(frag):
        print("Signal and Fragment lengths mismatch: {} {}".format(
                len(sig_float), len(frag)))
        return False
    for i, s in enumerate(splat.Signal(frag, float_value)):
        if s != float_tuple:
            print("Incorrect float signal value[{}]: {} {}".format(
                    i, s, float_tuple))
            return False
    if sig_float[x1] != sig_float[x2] != float_tuple:
        print("Incorrect float signal indexed values")
        return False
    func = lambda x: x * 0.1469
    (y1, y2) = (func(x) for x in (x1, x2))
    for i, (y,) in enumerate(splat.Signal(frag, func)):
        x = func(frag.n2s(i))
        if not floatcmp(x, y):
            print("Incorrect function signal value[{}]: {} {}".format(i, x, y))
            return False
    frag2 = splat.data.Fragment(duration=duration, channels=1)
    splat.sources.sine(frag2, 456.789, 0.0, 0.0)
    for x, y in zip(frag2, splat.Signal(frag, frag2)):
        if x != y:
            print("Incorrect fragment signal value: {} {}".format(x, y))
            return False
    for i, (y, z) in enumerate(splat.Signal(frag, (func, float_value))):
        x = func(frag.n2s(i))
        if not floatcmp(x, y) or z != float_value:
            print("Incorrect mixed signal value: {} {}".format(
                    (x, float_value), (y, z)))
            return False
    return True
set_id(test_signal, "Signal")

def test_sine():
    freq = 1237.9
    frag_float = splat.data.Fragment(duration=1.0)
    splat.sources.sine(frag_float, -0.5, freq, 0.0)
    frag_signal = splat.data.Fragment(duration=1.0)
    splat.sources.sine(frag_signal, -0.5, freq, lambda x: 0.0)
    frag_freq = splat.data.Fragment(duration=1.0)
    frag_freq.offset(freq)
    frag_frag = splat.data.Fragment(duration=1.0)
    splat.sources.sine(frag_frag, -0.5, frag_freq, 0.0)
    return check_multiple_md5(
        [frag_float, frag_signal, frag_frag],
        '46a8962a759033371f45c4ade9f2bfbd')
set_id(test_sine, "Sine source")

def test_sine_gen():
    gen = splat.gen.SineGenerator()
    f = 1000.0
    gen.run(0.0, 1.0, f)
    n = int(0.1234 * gen.frag.duration * gen.frag.sample_rate)
    s = math.sin(2 * math.pi * f * float(n) / gen.frag.sample_rate)
    return (check_samples(gen.frag, {n: (s, s)}) and
            check_md5(gen.frag, 'ec18389e198ee868d61c9439343a3337'))
set_id(test_sine_gen, "SineGenerator")

def test_square():
    freq = 1237.9
    frag_float = splat.data.Fragment(duration=1.0)
    splat.sources.square(frag_float, -0.5, freq, 0.0)
    frag_signal = splat.data.Fragment(duration=1.0)
    splat.sources.square(frag_signal, -0.5, freq, lambda x: 0.0)
    frag_freq = splat.data.Fragment(duration=1.0)
    frag_freq.offset(freq)
    frag_frag = splat.data.Fragment(duration=1.0)
    splat.sources.square(frag_frag, -0.5, frag_freq, 0.0)
    return check_multiple_md5(
        [frag_float, frag_signal, frag_frag],
        '6a6ab2e991baf48a6fe2c1d18700e40e')
set_id(test_square, "Square source")

def test_square_gen():
    gen = splat.gen.SquareGenerator()
    f = 1000.0
    gen.run(0.0, 1.0, f)
    nf = gen.frag.sample_rate / f
    samples = {int(nf * 0.1): (1.0, 1.0), int(nf * 0.9): (-1.0, -1.0)}
    return (check_samples(gen.frag, samples) and
            check_md5(gen.frag, '0ca047e998f512280800012b05107c63'))
set_id(test_square_gen, "SquareGenerator")

def test_triangle():
    freq = 1237.5
    frag_float = splat.data.Fragment(duration=1.0)
    splat.sources.triangle(frag_float, -0.5, freq, 0.0)
    frag_signal = splat.data.Fragment(duration=1.0)
    splat.sources.triangle(frag_signal, -0.5, freq, lambda x: 0.0)
    frag_freq = splat.data.Fragment(duration=1.0)
    frag_freq.offset(freq)
    frag_frag = splat.data.Fragment(duration=1.0)
    splat.sources.triangle(frag_frag, -0.5, frag_freq, 0.0)
    return check_multiple_md5(
        [frag_float, frag_signal, frag_frag],
        '4bce3885732ba2f5450e79e42155adaa')
set_id(test_triangle, "Triangle source")

def test_triangle_gen():
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
set_id(test_triangle_gen, "TriangleGenerator")

def test_overtones():
    freq = 1237.5
    ot = [(1.3, 0.0, -2.5), (5.7, 10.0, -12.9)]
    frag_float = splat.data.Fragment(duration=1.0)
    splat.sources.overtones(frag_float, -0.5, freq, 0.0, ot)
    frag_mixed = splat.data.Fragment(duration=1.0)
    splat.sources.overtones(frag_mixed, -0.5, freq, lambda x: 0.0, ot)
    frag_signal = splat.data.Fragment(duration=1.0)
    ot_signal = [(1.3, 0.0, -2.5), (5.7, lambda x: 10.0, -12.9)]
    splat.sources.overtones(frag_signal, -0.5, freq, lambda x: 0.0, ot_signal)
    frag_freq = splat.data.Fragment(duration=1.0)
    frag_freq.offset(1237.5)
    frag_frag = splat.data.Fragment(duration=1.0)
    splat.sources.overtones(frag_frag, -0.5, frag_freq, 0.0, ot)
    return check_multiple_md5(
        [frag_float, frag_mixed, frag_signal, frag_frag],
        '8974a1eea0db97af1aa171f531685e9d')
set_id(test_overtones, "Overtones source")

def test_overtones_gen():
    gen = splat.gen.OvertonesGenerator()
    gen.ot_decexp(1.0)
    f = 1000.0
    gen.run(0.0, 1.0, f)
    return check_md5(gen.frag, 'ee045e012673ff7ed4ab9bd590b57368')
set_id(test_overtones_gen, "OvertonesGenerator")

def test_polynomial():
    k0, k1, k2, k3 = coefs = (2.345, 3.6, 6.5, 100)
    p = splat.interpol.Polynomial(coefs)
    if coefs != p.coefs:
        print("Polynomial coefs mismatch")
        return False
    d = p.derivative()
    dcoefs = (k1, (k2 * 2), (k3 * 3))
    if d.coefs != dcoefs:
        print("Derivative error: {} instead of {}".format(d.coefs, dcoefs))
        return False
    i = d.integral(k0)
    if i.coefs != coefs:
        print("Integral error: {} instead of {}".format(i.coefs, coefs))
        return False
    return True
set_id(test_polynomial, "Polynomial")

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

def test_particle():
    start = 1.2
    end = 3.4
    f = 1234.56
    p = splat.gen.Particle(start, end, splat.lin2dB(f))
    if (p.start != start) or (p.end != end):
        print("Start/end times error")
    if not floatcmp(p.freq, f):
        print("Frequency conversion error: {} instead of {}".format(f, p.freq))
    return True
set_id(test_particle, "Particle")

def test_particle_pool():
    min_f = 123.45
    min_f_log = splat.lin2dB(min_f)
    max_f = 678.9
    max_f_log = splat.lin2dB(max_f)
    min_len = 0.1
    max_len = 0.3
    envelope = splat.interpol.Spline([(0.5, 0.0), (1.3, 1.0), (2.1, 0.0)])
    n_slices = 20
    density = 200
    count = 196
    half_count = 103
    pool = splat.gen.ParticlePool(min_f_log, max_f_log, min_len, max_len,
                                  envelope, n_slices, density)
    if (pool.start != envelope.start) or (pool.end != envelope.end):
        print("Invalid start/end times")
        return False
    if pool.count() != count: # plain magic
        print("Invalid number of particles: {}".format(pool.count()))
        return False
    for p in pool.iterate(share=0.5):
        if (p.length < min_len) or (p.length > max_len):
            print("Invalid particle length")
            return False
        if (p.start < envelope.start) or (p.start > envelope.end):
            print("Invalid particle start/end times")
            return False
        if (p.freq < min_f) or (p.freq > max_f):
            print("Invalid particle frequency")
            return False
    if abs(pool.count() - half_count) > (count / 35.0):
        print("Invalid number of particles left: {}".format(pool.count()))
        return False
    return True
set_id(test_particle_pool, "ParticlePool")

def test_particle_gen():
    class TestGenerator(splat.gen.Generator):
        def __init__(self, start, end, f_min, f_max, *args, **kw):
            super(TestGenerator, self).__init__(*args, **kw)
            self._start = start
            self._end = end
            self._f_min = f_min
            self._f_max = f_max

        def run(self, start, end, freq, levels):
            if (start < self._start) or (end > self._end):
                print(start, end, self._start, self._end)
                raise Exception("Invalid subgen start/end run times")
            if (freq < self._f_min) or (freq > self._f_max):
                raise Exception("Invalid subgen frequency")

    z_start = 0.5
    z_mid = 1.8
    z_end = 2.5
    eq_min_f = 100.0
    eq_mid_f = 2345.0
    eq_max_f = 8000.0
    subgen = TestGenerator(z_start, z_end, eq_min_f, eq_max_f)
    pgen = splat.gen.ParticleGenerator(subgen)
    if pgen.subgen is not subgen:
        print("Failed to get ParticleGenerator.subgen")
        return False
    pgen.set_z([(z_start, 0.0), (z_mid, 1.0, 0.0), (z_end, 0.0)])
    if (pgen.z.start != z_start) or (pgen.z.end != z_end):
        print("Incorrect envelope start/end times")
        return False
    pgen.set_eq([(eq_min_f, -30.0), (eq_mid_f, 10.0), (eq_max_f, -50)])
    if (pgen.eq.start != eq_min_f) or (pgen.eq.end != eq_max_f):
        print("Incorrect eq start/end frequencies")
        return False
    gfuzz = (3.0, 1.0)
    pgen.gain_fuzz = gfuzz
    if pgen.gain_fuzz != gfuzz:
        print("Failed to get gain fuzz")
        return False
    pgen.make_pool()
    if pgen.pool.count() != 126: # plain magic again
        print("Incorrect number of particles")
        return False
    if not floatcmp(pgen.curve(123.45, 456.78, 56.78), 154.33118):
        print("ParticleGenerator.curve error")
        return False
    pgen.do_show_progress = False
    pgen.run(0.0, 5.0, 1234.56)
    if pgen.pool.count() != 0:
        print("Not all particles were consumed in the run")
        return False
    return True
set_id(test_particle_gen, "ParticleGenerator")

# -----------------------------------------------------------------------------
# main function

def main(argv):
    print("Sample precision: {} bits".format(splat.sample_precision()))

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
