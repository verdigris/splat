import sys
import random
import math
import geomusic as geo

# Fixed parameters
gap = 0.1 # silence gap
l = 5.0 # total length in seconds
p = 8.0 # sounds per seconds
base = 110.0 # base frequency in Hz
m = 2 ** 0.5 # frequency span coeff
times = 5 # number of iterations

modes = {
    'lin': (m, lambda b, k, r: b * (1 + (k * r))),
    'log': ((2 ** m), lambda b, k, r: b * math.log((k * r), 2)),
    'exp': (math.log(m, 2), lambda b, k, r: b * 2 ** (k * r)),
}

def main(argv):
    if len(argv) > 1:
        mode = argv[1]
    else:
        mode = 'exp'

    if len(argv) > 2:
        file_name = argv[2]
    else:
        file_name = 'test.wav'

    print("mode: {0}, file name: {1}".format(mode, file_name))

    frag = geo.Fragment(2, 48000, (l + (2 * gap)))
    gen = geo.Generator(frag)
    gen.chain = geo.FilterChain([(geo.filters.linear_fade, (0.008,))])
    k, func = modes[mode]
    random.seed()
    rand = random.random
    end = l - gap
    b = base

    for i in xrange(times):
        t = gap
        while t < end:
            f = func(b, k, rand())
            t += (0.5 + rand()) / p
            start = min(t, end)
            t += (0.5 + rand()) / p
            stop = min(t, end)
            gen.sine(f, start, stop, (rand(), rand()))
        b *= 1.5

    geo.filters.normalize(frag, 0.4)
    frag.save_to_file(file_name, 2)

    return True

if __name__ == '__main__':
    import time
    start = time.time()
    ret = main(sys.argv)
    stop = time.time()
    print("Total time: {0:.3f}".format(stop - start))
    if ret is True:
        sys.exit(0)
    else:
        sys.exit(1)
