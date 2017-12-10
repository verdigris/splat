# Splat - setup.py
#
# Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017
#   Guillaume Tucker <guillaume@mangoz.org>
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
import os.path
try:
    import setuptools
except ImportError:
    import distutils.core as setuptools
from setuptools import setup, Extension

# The manual needs to be generated with Sphinx but is only required when
# running the sdist command.
if os.path.exists('Splat.pdf'):
    data_files = [('.', ['Splat.pdf',]),]
elif len(sys.argv) > 1 and sys.argv[1] == 'sdist':
    raise Exception('Splat.pdf is missing')
else:
    data_files = []

setup(name='verdigris.mu-splat', version='1.6',
      description="Sound generator, synthesizer and editor",
      author="Guillaume Tucker",
      author_email="guillaume@mangoz.org",
      url="https://github.com/verdigris/splat",
      py_modules=['test', 'example', 'dew_drop'],
      ext_modules=[Extension('_splat',
                             sources=['_splat.c', 'signal.c', 'spline.c',
                                      'frag.c', 'source.c', 'filter.c',
                                      'sine_table.c', 'mmap.c'],
                             depends=['_splat.h'])],
      packages=['splat', 'splat.tools'],
      scripts=['tools/splat'],
      data_files=data_files,
      long_description=open('README.rst', 'rb').read(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Other Audience',
        'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: C',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Topic :: Multimedia :: Sound/Audio',
        'Topic :: Multimedia :: Sound/Audio :: Editors',
        'Topic :: Multimedia :: Sound/Audio :: Sound Synthesis',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        ])
