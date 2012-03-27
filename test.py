import sys
import random
import math
sys.path.append('..')
import geomusic as geo

# Fixed parameters
l = 3.0 # total length in seconds
p = 8.0 # sounds per seconds
b = 55.0 # base frequency in Hz
m = 2 ** 7 # frequency span coeff

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

    frag = geo.Fragment(2, 48000, l)
    gen = geo.Generator(frag)
    k, func = modes[mode]
    n = int(l * p)
    t = 0.0
    random.seed()

    for i in range(n):
        f = func(b, k, random.random())
        start = t
        stop = t + l / n
        t = stop
        gen.sine(f, start, stop, (random.random(), random.random()))

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
