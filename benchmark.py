import sys
import time
sys.path.append('..')
import geomusic as geo

def bm(name, call, *args, **kw):
    start = time.time()
    ret = call(*args, **kw)
    stop = time.time()
    time_ms = (stop - start) * 1000
    print('{0}: {1:.3f}ms'.format(name, time_ms))
    return ret

def main(argv):
    time_s = 10.0
    frag = bm('Fragment', geo.Fragment, 2, 48000, time_s)
    gen = bm('Generator', geo.Generator, frag, ((0.7, 0.7)))
    bm('sine', gen.sine, 1000, 0, time_s)
    raw = bm('get_raw_bytes', frag.get_raw_bytes, 2)
    return True

if __name__ == '__main__':
    ret = main(sys.argv)
    if ret is True:
        sys.exit(0)
    else:
        sys.exit(1)
