from distutils.core import setup, Extension

setup(name='splat', version='1.1',
      description="Splat - sound generator",
      author="Guillaume Tucker - http://verdigris.mu",
      author_email="guillaume@mangoz.org",
      ext_modules=[Extension('_splat', ['_splat.c'])],
      packages=['splat'])
