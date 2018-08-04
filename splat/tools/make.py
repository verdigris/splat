# Splat - splat/tools/make.py
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

The ``splat make`` command can be used to produce a GNU Makefile and
automatically run the standard ``make`` command to generate complex
Splat compositions that are made up of several parts.  Each part is
implemented in a separate Python module to generate a Fragment object
and store it in a file (in Splat Audio Format).  A main mixer module
then mixes all the individual file fragments together.  This saves
regenerating the parts of a composition that do not change while
working on other parts.  It can also speed up the overall process by
generating several fragments in parallel using ``make``'s standard
parallel build features (``-j`` option).  Dependencies between any
input data files, fragment modules and the mixer are also all managed
by ``make``.

The convention is to save the output of each fragment module in a file
with the name of the module and the ``.saf`` extension.  The final
output file is named after the top-level mixer module name.  For
example, if the mixer is ``mixer.py`` and the fragments are
``frag1.py`` and ``frag2.py``, the fragment sound files will be
``frag1.saf`` and ``frag2.saf`` and the final output file will be
``mixer.wav``.  The output audio format is WAV by default but can be
changed with the ``--format`` option.

Syntax
^^^^^^

This is the output of ``splat make --help``::

  usage: Generate Splat makefile and run make [-h] [--format FORMAT] [--print]
                                              [--opt OPT]
                                              main

  positional arguments:
    main                  name of the mixer module

  optional arguments:
    -h, --help            show this help message and exit
    --format FORMAT, -f FORMAT
                          output audio format
    --print               print makefile and do not run make
    --opt OPT, -o OPT     options to pass to make

The ``splat make`` command takes only one positional argument which is
the name of the top-level mixer module.  It will then import this
module and look for special members to gather all the information
needed to generate a makefile.  When called with ``--print``, it will
just print the generated makefile.  Othwerwise, it will pass it on to
``make`` in order to automatically generate the piece.  The ``--opt``
or ``-o`` option can be used to pass extra options to make.

For example, if the top-level mixer module is called ``mixer`` and to
enable parallel runs with 8 threads::

  splat make -o"-j8" mixer

Special module members
^^^^^^^^^^^^^^^^^^^^^^

Starting from the top-level mixer module, ``splat make`` will look for
the following module members:

``SPLAT_FRAGS``

  Dictionary with fragment module names as keys and a 2-tuple with the
  starting time and mixer level as values.

  It is typically used in the top-level mixer module to describe how
  to mix all the fragments together.  For example::

    SPLAT_FRAGS = {
        'beep': (0.0, dB(0)),
        'buzz': (1.0, dB(-9)),
    }

  This will create a dependency on ``beep.saf`` and ``buzz.saf`` to
  ensure they get generated first from ``beep.py`` and ``buzz.py``.
  Then beep will be mixed at the beginning of the output with a level
  of 0 dB and buzz with an offset of 1 second and a level of -9 dB.

``SPLAT_MODS``
  List with the names of other Python modules this module depends on.

  Dependencies will be recursively created within these modules.  It
  is typically used when one module imports another one, for example
  with settings or common code.  Modules which generate a fragment are
  not listed in ``SPLAT_MODS`` but in ``SPLAT_FRAGS`` instead in order
  to manage the dependency with their generated files.

``SPLAT_DEPS``
  List with arbitrary file names this module depends on.

  This is typically used to express the dependency of one module or
  fragment module on some input data files such as sound samples or
  any other kind of media.  Python module files may also be treated as
  plain files to avoid generating dependencies recursively by listing
  them in ``SPLAT_DEPS`` instead of ``SPLAT_MODS``.

A dependency chain may look like this:

  ``mixer.wav`` -> ``mixer.py`` -> ``frag1.saf`` -> ``frag1.py`` -> ``input-sample.wav``

Here, ``mixer.py`` is the top-level mixer which produces the final
``mixer.wav`` file.  It defines a ``SPLAT_FRAGS`` dictionary with
``frag1``, so ``frag1.py`` is run to generate ``frag1.saf``.  There's
an input sound file ``input-sample.wav`` which is not generated by
Splat but listed in ``SPLAT_DEPS`` in ``frag1.py`` as it uses it to
produce ``frag1.saf``.

So, the input sample is used by ``frag1.py`` to generate
``frag1.saf``.  Then ``mixer.py`` takes it and mixes it (typically
with other fragments) and generates ``mixer.wav``.

Main functions
^^^^^^^^^^^^^^

Since each module needs to be run by ``make`` in order to produce a
file, they need to define a main function.  The following standard
main functions for both mixer and fragment modules are provided:

.. autofunction:: splat.tools.make.main_mixer
.. autofunction:: splat.tools.make.main_frag
"""

from __future__ import print_function
import sys
import subprocess
import shlex
import argparse
import splat
import splat.data
import splat.filters

def build_dep_tree(mod, name, dep_tree, clean_list):
    # Generic files this module directly depends on
    deps = getattr(mod, 'SPLAT_DEPS', list())

    # Any Python modules this module depends on
    mods = getattr(mod, 'SPLAT_MODS', None)
    if mods is not None:
        for m in mods:
            target = '.'.join([m, 'py'])
            deps.append(target)
            if target not in dep_tree:
                m_mod = __import__(m)
                build_dep_tree(m_mod, target, dep_tree, clean_list)

    # Splat fragment modules this module depends on
    frags = getattr(mod, 'SPLAT_FRAGS', None)
    if frags is not None:
        if isinstance(frags, dict):
            frag_seq = frags.keys()
        else:
            frag_seq = frags
        for frag in frag_seq:
            target = '.'.join([frag, 'saf'])
            deps.append(target)
            clean_list.append(target)
            if target not in dep_tree:
                frag_mod = __import__(frag)
                build_dep_tree(frag_mod, target, dep_tree, clean_list)

    # Add this module to the tree with all the things it depends on
    dep_tree[name] = deps

def makefile(main, fmt):
    mk = """\
ver_check := $(shell python3 -c 'import splat; splat.check_version(({}, {}))' \\
	|| echo ERROR)

ifeq ($(ver_check),ERROR)
  $(error 'Unmet required Splat version')
endif

""".format(*splat.VERSION)

    mk += """\
%.saf: %.py
	@echo "  GEN     " $@
	@python3 $< $@

%.{}: %.py
	@echo "  MIX     " $@
	@python3 $< $@
""".format(fmt)

    main_mix = "{}.{}".format(main, fmt)
    dep_tree = dict()
    clean_list = [main_mix]
    main_mod = __import__(main)
    build_dep_tree(main_mod, main_mix, dep_tree, clean_list)
    mk += '\nall: {}\n\n'.format(main_mix)
    for target, deps in dep_tree.items():
        if deps:
            mk += target + ':'
            for dep in deps:
                mk += " \\\n\t{}".format(dep)
            mk += '\n\n'

    mk += """\

.PHONY: clean

clean:
"""
    for target in clean_list:
        mk += "\t@echo \"  CLEAN    {}\"\n".format(target)
        mk += "\t@rm -f {}\n".format(target)

    return mk

def run_make(mk, opts=None):
    args = ['make', '-f', '-']
    if opts:
        args += opts
    p = subprocess.Popen(args, stdin=subprocess.PIPE)
    p.communicate(input=mk.encode('utf-8'))
    return p.returncode

def main(argv):
    parser = argparse.ArgumentParser("Generate Splat makefile and run make")
    parser.add_argument('main',
                        help="name of the mixer module")
    parser.add_argument('--format', '-f', default="wav",
                        help="output audio format")
    parser.add_argument('--print', action='store_true', default=False,
                        help="print makefile and do not run make")
    parser.add_argument('--opt', '-o', help="options to pass to make")
    args = parser.parse_args(argv[1:])

    mk = makefile(args.main, args.format)
    if args.print:
        print(mk)
        return True
    else:
        opt = shlex.split(args.opt) if args.opt else None
        ret = run_make(mk, opt)
        return True if (ret == 0) else False

if __name__ == '__main__':
    ret = main(sys.argv)
    sys.exit(0 if ret is True else 1)

# standard main function for a mixer target
def main_mixer(argv, fragments, master=None, filters=None):
    """Main function to be used by a mixer module.

    The ``argv`` are the command line arguments, usually ``sys.argv``.
    Then ``fragments`` is a dictionary, usually ``SPLAT_FRAGS`` or one
    with a similar structure to mix all the fragments together and
    create the final audio file.  Optionally, a ``master`` Fragment
    object can be passed, otherwise the function will create one with
    default settings (stereo 48 kHz).  It's also possible to pass a
    list of ``filters`` to run a FilterChain as a Generator would do
    on the ``master`` fragment.  The function than saves the fragment
    in a file with the name given in ``argv[1]``.
    """
    f_out = argv[1]
    if master is None:
        master = splat.data.Fragment()
    for frag, (ofst, ml) in fragments.items():
        fname = '.'.join([frag, 'saf'])
        master.mix(splat.data.Fragment.open(fname), offset=ofst, levels=ml)
    if filters is not None:
        splat.filters.FilterChain(filters).run(master)
    master.save(f_out)

# standard main function for a fragment target
def main_frag(argv, func):
    """Main function to be used by a fragment module.

    The ``argv`` are the command line arguments, usually ``sys.argv``.
    This function will first run a call-back ``func`` function which
    should return a Fragment object.  It will then save this fragment
    into a file with the name passed in ``argv[1]`` which is usually
    the module name with the ``saf`` format as defined by convention
    in the Makefile.
    """
    f_out = argv[1]
    frag = func(argv[1:])
    frag.save(f_out)
