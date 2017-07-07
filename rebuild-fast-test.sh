#!/bin/sh

set -e

simd="$1"

if [ -z "$simd" ]; then
    arch=$(uname -m)

    case "$arch" in
        "armv7l")
            simd=neon-v7
            ;;
	"aarch64")
	    simd=neon-v8
	    ;;
        "x86_64")
            simd=sse
            ;;
        *)
            echo "Unknown architecture, please specify"
            ;;
    esac
fi

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
