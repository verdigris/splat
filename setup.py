try:
    import setuptools
    from setuptool import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

setup(name='verdigris.mu-splat', version='1.3',
      description="Splat - sound generator",
      author="Guillaume Tucker",
      author_email="guillaume@mangoz.org",
      url="https://github.com/verdigris/splat",
      py_modules=['test', 'example', 'dew_drop'],
      ext_modules=[Extension('_splat', ['_splat.c'])],
      packages=['splat'],
      data_files=[('.', ['Splat.pdf',]),],
      long_description="""
The basic idea is to apply mathematical concepts to
musical composition as well as sound synthesis. In practice,
Splat lets you create just about any sound you can imagine
and code in software.
""",
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
