import sys
import random
import math
sys.path.append('..')
import geomusic as geo

if len(sys.argv) > 1:
    mode = sys.argv[1]
else:
    mode = 'exp'

if len(sys.argv) > 2:
    file_name = sys.argv[2]
else:
    file_name = 'test.wav'

print("mode: {0}, file name: {1}".format(mode, file_name))

random.seed()
frag = geo.Fragment(2, 48000)
gen = geo.Generator(frag)
gen.set_levels((0.7, 0.7))

l = 7.0 # total length in seconds
p = 8.0 # sounds per seconds
b = 55.0 # base frequency in Hz
m = 2 ** 7 # frequency span coeff

modes = {
    'lin': (m, lambda b, k, r: b * (1 + (k * r))),
    'log': ((2 ** m), lambda b, k, r: b * math.log((k * r), 2)),
    'exp': (math.log(m, 2), lambda b, k, r: b * 2 ** (k * r)),
}

k, func = modes[mode]
n = int(l * p)
t = 0.0

for i in range(n):
    f = func(b, k, random.random())
    start = t
    stop = t + l / n
    t = stop
    gen.sine(f, start, stop)

frag.save_to_file(file_name, 2)

sys.exit(0)
