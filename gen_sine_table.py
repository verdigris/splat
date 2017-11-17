from __future__ import print_function
import sys
import math
import cmath
import argparse
import splat.interpol

def make_spline(n, m):
    step = cmath.pi / n
    sine_pts = [(x, math.sin(x)) for x in (y * step for y in range(n + 1))]
    return splat.interpol.spline(sine_pts, n=m)

def gen_table(pols, out=print):
    out("static const struct splat_sine_poly _sine_table[{}]".format(
            len(pols)))
    out("__attribute__ ((aligned (16))) = {");
    for i, (x0, x1, coefs) in enumerate(pols):
        out("\t{{{{ /* {} */".format(i))
        for coef in coefs.coefs:
            out("\t\t{:.55e},".format(coef))
        out("\t}},")
    out("};")
    out()
    out("const struct splat_sine_poly *splat_sine_table = _sine_table;")
    out("const size_t splat_sine_table_len = {};".format(len(pols)))
    out("const size_t splat_sine_table_mask = 0x{:x};".format(len(pols) - 1))

def main(argv):
    ranges = {
        'full': cmath.pi * 2.0, 'half': cmath.pi, 'quarter': cmath.pi / 2.0
    }
    parser = argparse.ArgumentParser("Sine function interpolation table")
    parser.add_argument('-m', default=3, type=int, help="Polynomials order")
    parser.add_argument('-n', default=64, type=int, help="Number of points")
    parser.add_argument('--range', default='half', choices=ranges,
                        help="Range of the points relative to 2 * PI")
    args = parser.parse_args(argv[1:])

    if math.log(args.n, 2) % 1.0:
        print("The number of points needs to be a power of 2")
        return False

    spline = make_spline(args.n, args.m)
    print("/* Automatically generated file - read the manual */")
    print()
    print("#include \"_splat.h\"")
    print()
    gen_table(spline.pols)

    return True

if __name__ == '__main__':
    res = main(sys.argv)
    sys.exit(0 if res is True else 1)
