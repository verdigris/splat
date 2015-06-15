#!/bin/sh

set -e

simd="${1-sse}"
[ -x build-"$simd" ] || {
    echo "Unsupported SIMD mode: $simd"
    exit 1
}

rm -rf build
./build-std
python test.py
python fast_test.py

rm -rf build
./build-"$simd"
python fast_test.py --compare

exit 0
