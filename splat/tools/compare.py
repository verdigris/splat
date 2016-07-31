# Splat - splat/tools/compare.py
#
# Copyright (C) 2015, 2016 Guillaume Tucker <guillaume@mangoz.org>
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

"""Summary
^^^^^^^

The ``splat compare`` command is a simple utility to compare two sound
files.  It will save the difference created by substracting
corresponding samples from one file to the other and print the peak
amplitude of the difference in dB.

Syntax
^^^^^^

This is the output of ``splat compare --help``::

  usage: Compare sound files and show the peak signal difference
         [-h] [--save SAVE] f1 f2

  positional arguments:
    f1           path to first sound file
    f2           path to second sound file

  optional arguments:
    -h, --help   show this help message and exit
    --save SAVE  path to save the difference

For example, to compare two files, print the peak difference in dB and
save the result of the substracted fragments in ``delta.wav``::

  $ splat compare file-1.wav file-2.wav --save=delta.wav
  Peak delta: -45.803 dB

Functions
^^^^^^^^^

The following function can also be imported to perform the same
comparision between Fragment objects directly in Python instead of via
the command line:

.. autofunction:: splat.tools.compare.frag_delta
"""

import sys
import argparse
import splat
import splat.data

def frag_delta(f1, f2):
    """Return a fragment with the difference between two fragments.

    Fragment ``f2`` gets substracted from ``f1`` and a new fragment
    with the difference is returned.  The input fragments are not
    modified by this operation.
    """
    delta = f2.dup()
    delta.amp(-1)
    delta.mix(f1)
    return delta

def file_peak_delta_dB(f1_path, f2_path, save_path=None):
    f1 = splat.data.Fragment.open(f1_path)
    f2 = splat.data.Fragment.open(f2_path)
    delta = frag_delta(f1, f2)
    if save_path:
        delta.save(save_path, normalize=False)
    return splat.lin2dB(delta.get_peak()[0]['peak'])

def main(argv):
    parser = argparse.ArgumentParser(
        "Compare sound files and show the peak signal difference")
    parser.add_argument('f1', help="path to first sound file")
    parser.add_argument('f2', help="path to second sound file")
    parser.add_argument('--save', help="path to save the difference")
    args = parser.parse_args(argv[1:])

    peak_delta = file_peak_delta_dB(args.f1, args.f2, args.save)
    print("Peak delta: {:.3f} dB".format(peak_delta))
    return True

if __name__ == '__main__':
    ret = main(sys.argv)
    sys.exit(0 if ret is True else 1)
