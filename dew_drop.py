import sys
import geomusic as geo

def main(argv):
    gen = geo.Generator(geo.Fragment(2, 48000, 18.0))
    gen.chain = geo.FilterChain([geo.filters.linear_fade])
    gen.time_stretch = 1.8
    s = geo.scales.LogScale(fund=440.0)

    # print all notes frequencies
    for octave in range(-2, 1):
        for note in ['A', 'B', 'C#', 'D', 'E', 'F#', 'G#']:
            note_name = "{0}{1}".format(note, octave)
            print('{0:4s}: {1:.3f}'.format(note_name, s[note_name]))

    # voice 1
    gen.levels = (0.9, 0.9)
    gen.sine(s['A-2'], 0.00, 0.90)
    gen.sine(s['D-2'], 0.90, 1.80)
    gen.sine(s['F#-2'], 1.80, 2.36)
    gen.sine(s['E-2'], 2.36, 3.26)
    gen.sine(s['D-2'], 3.26, 3.81)
    gen.sine(s['F#-2'], 3.81, 4.72)
    gen.sine(s['E-2'], 4.72, 5.28)
    gen.sine(s['B-2'], 5.28, 6.18)
    gen.sine(s['E-2'], 6.18, 7.08)
    gen.sine(s['A-2'], 7.08, 7.64)
    gen.sine(s['D-2'], 7.64, 8.19)
    gen.sine(s['E-2'], 8.19, 9.10)
    gen.sine(s['A-2'], 9.10, 10.00)

    # voice 2
    gen.levels = (1.0, 0.7)
    gen.sine(s['C#-1'], 0.00, 0.56)
    gen.sine(s['E-1'], 0.56, 0.90)
    gen.sine(s['A-1'], 0.90, 1.46)
    gen.sine(s['F#-2'], 1.46, 1.80)
    gen.sine(s['A-1'], 1.80, 2.14)
    gen.sine(s['C#-1'], 2.14, 2.36)
    gen.sine(s['B-1'], 2.36, 2.70)
    gen.sine(s['D-1'], 2.70, 3.26)
    gen.sine(s['B-1'], 3.26, 3.81)
    gen.sine(s['C#-1'], 3.81, 4.36)
    gen.sine(s['A-1'], 4.36, 4.72)
    gen.sine(s['G#-2'], 4.72, 5.28)
    gen.sine(s['D-1'], 5.28, 5.83)
    gen.sine(s['F#-1'], 5.83, 6.18)
    gen.sine(s['B-1'], 6.18, 6.53)
    gen.sine(s['D-1'], 6.53, 7.08)
    gen.sine(s['E-1'], 7.08, 7.64)
    gen.sine(s['F#-1'], 7.64, 7.98)
    gen.sine(s['A-1'], 7.98, 8.54)
    gen.sine(s['D-1'], 8.54, 9.10)
    gen.sine(s['C#-1'], 9.10, 10.00)

    # voice 3
    gen.levels = (0.7, 1.0)
    gen.sine(s['E'], 0.00, 0.34)
    gen.sine(s['D'], 0.34, 0.56)
    gen.sine(s['E'], 0.56, 0.90)
    gen.sine(s['D'], 0.90, 1.24)
    gen.sine(s['C#'], 1.80, 2.14)
    gen.sine(s['A'], 2.14, 2.36)
    gen.sine(s['B'], 2.36, 2.70)
    gen.sine(s['D'], 2.70, 2.92)
    gen.sine(s['E'], 2.92, 3.26)
    gen.sine(s['F#'], 3.26, 3.81)
    gen.sine(s['D'], 4.72, 4.93)
    gen.sine(s['E'], 4.93, 5.28)
    gen.sine(s['F#'], 5.28, 5.83)
    gen.sine(s['A'], 5.83, 6.18)
    gen.sine(s['B'], 6.18, 6.53)
    gen.sine(s['C#'], 6.53, 6.74)
    gen.sine(s['D'], 6.74, 7.08)
    gen.sine(s['C#'], 7.08, 7.64)
    gen.sine(s['B'], 7.98, 8.54)
    gen.sine(s['G#-1'], 8.54, 9.10)
    gen.sine(s['A'], 9.10, 10.00)

    geo.filters.normalize(gen.frag, 0.95)
    gen.frag.save_to_file('dew_drop.wav', 2)

if __name__ == '__main__':
    main(sys.argv)
    sys.exit(0)
