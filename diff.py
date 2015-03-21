# Splat - dew_drop.py
#
# Copyright (C) 2015 Guillaume Tucker <guillaume@mangoz.org>
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
import argparse
import splat
import splat.data

def frag_diff(f1, f2):
    delta = f2.dup()
    delta.amp(-1)
    delta.mix(f1)
    return delta

def main(argv):
    parser = argparse.ArgumentParser(
        "Compute the difference between sound files")
    parser.add_argument('f1', help="path to first sound file")
    parser.add_argument('f2', help="path to second sound file")
    parser.add_argument('--save', help="path to save the difference")
    args = parser.parse_args(argv[1:])

    f1 = splat.data.Fragment.open(args.f1)
    f2 = splat.data.Fragment.open(args.f2)
    delta = frag_diff(f1, f2)
    print("Peak error: {:.3f} dB".format(
            splat.lin2dB(delta.get_peak()[0]['peak'])))
    if args.save:
        delta.save(args.save, normalize=False)
    return True

if __name__ == '__main__':
    ret = main(sys.argv)
    sys.exit(0 if ret is True else 1)
