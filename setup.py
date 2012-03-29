from distutils.core import setup, Extension

setup(name='_geomusic', version='0.0',
      ext_modules=[Extension('_geomusic', ['_geomusic.c'])],
      packages=['geomusic'])
