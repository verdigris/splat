from _geomusic import lin2dB, dB2lin
from generator import Generator, FilterChain, SineGenerator, OvertonesGenerator
from data import Fragment
from filters import FilterChain
import filters
import sources
import scales

VERSION = (1, 0)
VERSION_STR = '.'.join([str(v) for v in VERSION])
BUILD = 1
