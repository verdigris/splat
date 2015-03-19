import sys
import math
import cmath
import splat.gen
import splat.data
import splat
from splat import dB2lin as dB

def main(argv):
    print("Sample precision: {}-bit".format(splat.SAMPLE_WIDTH))

    freq = 123.45
    duration = 30.0
    phase = 23.456
    pts = [0.0, 0.2, 0.234, 0.456, 0.45602, 0.7124, 0.89, 0.90001]
    overtones = [(1.0, 0.0, 0.45), (2.58, 12.345, 0.45), (200.0, 0.0, 1.0)]

    # -------------------------------------------------------------------------

    frag = splat.data.Fragment(channels=1)
    gen = splat.gen.SineGenerator(frag=frag)
    gen.run(0.0, duration, freq, phase)
    frag.save('sine-{}.wav'.format(splat.SAMPLE_WIDTH), normalize=False)

    print('sine')
    for t in pts:
        t = duration * t
        n = frag.s2n(t)
        t = float(n) / frag.rate
        st = 2 * cmath.pi * freq * (t + phase)
        y1 = frag[n][0]
        y2 = math.sin(st)
        delta = abs(y2 - y1)
        if delta != 0.0:
            delta_dB = "{:.3f}".format(splat.lin2dB(delta))
        else:
            delta_dB = 'infinity'
        print(t, y1, y2, delta_dB)

    # -------------------------------------------------------------------------

    frag = splat.data.Fragment(channels=1)
    gen = splat.gen.OvertonesGenerator(frag=frag)
    gen.overtones = overtones
    gen.run(0.0, duration, freq)
    frag.save('overtones-{}.wav'.format(splat.SAMPLE_WIDTH), normalize=False)

    max_ratio = frag.rate / 2.0 / freq
    ot_clipped = []
    for ot in gen.overtones:
        if ot[0] < max_ratio:
            ot_clipped.append(ot)

    print('overtones')
    for t in pts:
        t = duration * t
        n = frag.s2n(t)
        y1 = frag[n][0]
        y2 = sum(a * math.sin(2 * cmath.pi * freq * r * (t + ph))
                 for r, ph, a in ot_clipped)
        delta = abs(y2 - y1)
        if delta != 0.0:
            delta_dB = "{:.3f}".format(splat.lin2dB(delta))
        else:
            delta_dB = 'infinity'
        print(t, y1, y2, delta_dB)

    return True

if __name__ == '__main__':
    res = main(sys.argv)
    sys.exit(0 if res is True else 1)
