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
            frag_seq = frags.iterkeys()
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
ver_check := $(shell python -c 'import splat; splat.check_version(({}, {}))' \\
	|| echo ERROR)

ifeq ($(ver_check),ERROR)
  $(error 'Unmet required Splat version')
endif

""".format(*splat.VERSION)

    mk += """\
%.saf: %.py
	@echo "  GEN     " $@
	@python $< $@

%.{}: %.py
	@echo "  MIX     " $@
	@python $< $@
""".format(fmt)

    main_mix = "{}.{}".format(main, fmt)
    dep_tree = dict()
    clean_list = [main_mix]
    main_mod = __import__(main)
    build_dep_tree(main_mod, main_mix, dep_tree, clean_list)
    mk += '\nall: {}\n\n'.format(main_mix)
    for target, deps in dep_tree.iteritems():
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
    p.communicate(input=mk)
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
    f_out = argv[1]
    if master is None:
        master = splat.data.Fragment()
    for frag, (ofst, ml) in fragments.iteritems():
        fname = '.'.join([frag, 'saf'])
        master.mix(splat.data.Fragment.open(fname), offset=ofst, levels=ml)
    if filters is not None:
        splat.filters.FilterChain(filters).run(master)
    master.save(f_out)

# standard main function for a fragment target
def main_frag(argv, func):
    f_out = argv[1]
    frag = func(argv[1:])
    frag.save(f_out)
