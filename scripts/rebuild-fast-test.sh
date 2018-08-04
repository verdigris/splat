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

fast_cmd=scripts/build-"$simd"

[ -x "$fast_cmd" ] || {
    echo "Unsupported SIMD mode: $simd"
    exit 1
}

rm -rf build
./scripts/build-std
python3 test.py
python3 fast_test.py

rm -rf build
"$fast_cmd"
python3 fast_test.py --compare

exit 0
