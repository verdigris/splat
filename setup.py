from distutils.core import setup, Extension

setup(name='splat', version='0.0',
      ext_modules=[Extension('_splat', ['_splat.c'])],
      packages=['splat'])
