# Geomusic - dew_drop.py
#
# Copyright (C) 2012, 2013 Guillaume Tucker <guillaume@mangoz.org>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import time
import geomusic

def set_fade(gen, duration):
    gen.filters = geomusic.FilterChain(
        [(geomusic.filters.linear_fade, (duration,))])

def main(argv):
    gen = geomusic.OvertonesGenerator(geomusic.Fragment(2, 48000, 18.0))
    gen.time_stretch = 1.8
    s = geomusic.scales.LogScale(fund=440.0)

    # print frequencies of all the notes of the scale over 3 octaves
    for octave in range(-2, 1):
        for note in ['A', 'B', 'C#', 'D', 'E', 'F#', 'G#']:
            note_name = "{0}{1}".format(note, octave)
            print('{0:4s}: {1:.3f}'.format(note_name, s[note_name]))
        print("-------------")

    print("Voice 1")
    gen.levels = (-2.5, -2.5)
    gen.ot_decexp(2.0)
    set_fade(gen, 0.04)
    gen.run(s['A-2'], 0.00, 0.90)
    gen.run(s['D-2'], 0.90, 1.80)
    gen.run(s['F#-2'], 1.80, 2.36)
    gen.run(s['E-2'], 2.36, 3.26)
    gen.run(s['D-2'], 3.26, 3.81)
    gen.run(s['F#-2'], 3.81, 4.72)
    gen.run(s['E-2'], 4.72, 5.28)
    gen.run(s['B-2'], 5.28, 6.18)
    gen.run(s['E-2'], 6.18, 7.08)
    gen.run(s['A-2'], 7.08, 7.64)
    gen.run(s['D-2'], 7.64, 8.19)
    gen.run(s['E-2'], 8.19, 9.10)
    gen.run(s['A-2'], 9.10, 10.00)

    print("Voice 2")
    gen.levels = (0.0, -2.5)
    gen.ot_decexp(1.6)
    set_fade(gen, 0.02)
    gen.run(s['C#-1'], 0.00, 0.56)
    gen.run(s['E-1'], 0.56, 0.90)
    gen.run(s['A-1'], 0.90, 1.46)
    gen.run(s['F#-2'], 1.46, 1.80)
    gen.run(s['A-1'], 1.80, 2.14)
    gen.run(s['C#-1'], 2.14, 2.36)
    gen.run(s['B-1'], 2.36, 2.70)
    gen.run(s['D-1'], 2.70, 3.26)
    gen.run(s['B-1'], 3.26, 3.81)
    gen.run(s['C#-1'], 3.81, 4.36)
    gen.run(s['A-1'], 4.36, 4.72)
    gen.run(s['G#-2'], 4.72, 5.28)
    gen.run(s['D-1'], 5.28, 5.83)
    gen.run(s['F#-1'], 5.83, 6.18)
    gen.run(s['B-1'], 6.18, 6.53)
    gen.run(s['D-1'], 6.53, 7.08)
    gen.run(s['E-1'], 7.08, 7.64)
    gen.run(s['F#-1'], 7.64, 7.98)
    gen.run(s['A-1'], 7.98, 8.54)
    gen.run(s['D-1'], 8.54, 9.10)
    gen.run(s['C#-1'], 9.10, 10.00)

    print("Voice 3")
    gen.levels = (-2.5, 0.0)
    gen.ot_decexp(1.2)
    set_fade(gen, 0.015)
    gen.run(s['E'], 0.00, 0.34)
    gen.run(s['D'], 0.34, 0.56)
    gen.run(s['E'], 0.56, 0.90)
    gen.run(s['D'], 0.90, 1.24)
    gen.run(s['C#'], 1.80, 2.14)
    gen.run(s['A'], 2.14, 2.36)
    gen.run(s['B'], 2.36, 2.70)
    gen.run(s['D'], 2.70, 2.92)
    gen.run(s['E'], 2.92, 3.26)
    gen.run(s['F#'], 3.26, 3.81)
    gen.run(s['D'], 4.72, 4.93)
    gen.run(s['E'], 4.93, 5.28)
    gen.run(s['F#'], 5.28, 5.83)
    gen.run(s['A'], 5.83, 6.18)
    gen.run(s['B'], 6.18, 6.53)
    gen.run(s['C#'], 6.53, 6.74)
    gen.run(s['D'], 6.74, 7.08)
    gen.run(s['C#'], 7.08, 7.64)
    gen.run(s['B'], 7.98, 8.54)
    gen.run(s['G#-1'], 8.54, 9.10)
    gen.run(s['A'], 9.10, 10.00)

    if (len(argv) > 1) and (argv[1] == 'reverb'):
        print("Reverb...")
        d = [((float(t) / 400), -(6 + float(t) * 0.2)) for t in range(800)]
        geomusic.filters.reverb(gen.frag, d)

    print("Saving to file...")
    padded = geomusic.Fragment(2, 48000, (gen.frag.duration + 1.0))
    padded.mix(gen.frag, 0.5)
    padded.normalize(-0.1)
    padded.save_to_file('dew_drop.wav')


if __name__ == '__main__':
    start = time.time()
    main(sys.argv)
    stop = time.time()
    print("Total time: {0:.3f}".format(stop - start))
    sys.exit(0)
