# Splat - dew_drop.py
#
# Copyright (C) 2012, 2013 Guillaume Tucker <guillaume@mangoz.org>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
import time
import splat
import splat.data
import splat.gen
import splat.filters
import splat.scales

def set_fade(gen, duration):
    gen.filters = splat.filters.FilterChain(
        [(splat.filters.linear_fade, (duration,))])

def main(argv):
    gen = splat.gen.OvertonesGenerator(splat.data.Fragment(2, 48000, 18.0))
    gen.time_stretch = 1.8
    s = splat.scales.LogScale(fund=440.0)

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
    gen.run(0.00, 0.90, s['A-2'])
    gen.run(0.90, 1.80, s['D-2'])
    gen.run(1.80, 2.36, s['F#-2'])
    gen.run(2.36, 3.26, s['E-2'])
    gen.run(3.26, 3.81, s['D-2'])
    gen.run(3.81, 4.72, s['F#-2'])
    gen.run(4.72, 5.28, s['E-2'])
    gen.run(5.28, 6.18, s['B-2'])
    gen.run(6.18, 7.08, s['E-2'])
    gen.run(7.08, 7.64, s['A-2'])
    gen.run(7.64, 8.19, s['D-2'])
    gen.run(8.19, 9.10, s['E-2'])
    gen.run(9.10, 10.00, s['A-2'])

    print("Voice 2")
    gen.levels = (0.0, -2.5)
    gen.ot_decexp(1.6)
    set_fade(gen, 0.02)
    gen.run(0.00, 0.56, s['C#-1'])
    gen.run(0.56, 0.90, s['E-1'])
    gen.run(0.90, 1.46, s['A-1'])
    gen.run(1.46, 1.80, s['F#-2'])
    gen.run(1.80, 2.14, s['A-1'])
    gen.run(2.14, 2.36, s['C#-1'])
    gen.run(2.36, 2.70, s['B-1'])
    gen.run(2.70, 3.26, s['D-1'])
    gen.run(3.26, 3.81, s['B-1'])
    gen.run(3.81, 4.36, s['C#-1'])
    gen.run(4.36, 4.72, s['A-1'])
    gen.run(4.72, 5.28, s['G#-2'])
    gen.run(5.28, 5.83, s['D-1'])
    gen.run(5.83, 6.18, s['F#-1'])
    gen.run(6.18, 6.53, s['B-1'])
    gen.run(6.53, 7.08, s['D-1'])
    gen.run(7.08, 7.64, s['E-1'])
    gen.run(7.64, 7.98, s['F#-1'])
    gen.run(7.98, 8.54, s['A-1'])
    gen.run(8.54, 9.10, s['D-1'])
    gen.run(9.10, 10.00, s['C#-1'])

    print("Voice 3")
    gen.levels = (-2.5, 0.0)
    gen.ot_decexp(1.2)
    set_fade(gen, 0.015)
    gen.run(0.00, 0.34, s['E'])
    gen.run(0.34, 0.56, s['D'])
    gen.run(0.56, 0.90, s['E'])
    gen.run(0.90, 1.24, s['D'])
    gen.run(1.80, 2.14, s['C#'])
    gen.run(2.14, 2.36, s['A'])
    gen.run(2.36, 2.70, s['B'])
    gen.run(2.70, 2.92, s['D'])
    gen.run(2.92, 3.26, s['E'])
    gen.run(3.26, 3.81, s['F#'])
    gen.run(4.72, 4.93, s['D'])
    gen.run(4.93, 5.28, s['E'])
    gen.run(5.28, 5.83, s['F#'])
    gen.run(5.83, 6.18, s['A'])
    gen.run(6.18, 6.53, s['B'])
    gen.run(6.53, 6.74, s['C#'])
    gen.run(6.74, 7.08, s['D'])
    gen.run(7.08, 7.64, s['C#'])
    gen.run(7.98, 8.54, s['B'])
    gen.run(8.54, 9.10, s['G#-1'])
    gen.run(9.10, 10.00, s['A'])

    if (len(argv) > 1) and (argv[1] == 'reverb'):
        print("Reverb...")
        d = splat.filters.reverb_delays()
        splat.filters.reverb(gen.frag, d)

    print("Saving to file...")
    padded = splat.data.Fragment(2, 48000, (gen.frag.duration + 1.0))
    padded.mix(gen.frag, 0.5)
    padded.normalize()
    padded.save('dew_drop.wav')


if __name__ == '__main__':
    start = time.time()
    main(sys.argv)
    stop = time.time()
    print("Total time: {0:.3f}".format(stop - start))
    sys.exit(0)
