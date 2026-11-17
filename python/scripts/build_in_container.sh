#!/bin/sh
# Runs build_wheel.sh inside a manylinux or musllinux container, with the repository mounted at /work:
#   docker run --rm -v "$PWD:/work" -w /work -e EXPECTED_TAG -e ASSET -e VERSION IMAGE sh python/scripts/build_in_container.sh
set -eu

# The images carry CPython builds under /opt/python; pip gives them uv, cmake and ninja without a download script.
toolchain=$(ls -d /opt/python/cp313-*/bin | head -n 1)
"$toolchain/python" -m pip install --quiet uv cmake ninja
PATH="$toolchain:$PATH"
export PATH

# The project needs GCC 15, and neither image ships it (both carry GCC 14). manylinux_2_28 (AlmaLinux 8) has it as
# gcc-toolset-15 in its own repositories; musllinux_1_2 (Alpine 3.22) takes g++ from Alpine 3.23, the first branch
# with GCC 15. Neither changes the libc the library links, so the platform tags stay. build_native.py stops the
# build if the GCC it finds is older than 15.
if command -v dnf >/dev/null; then
    dnf install -y -q gcc-toolset-15-gcc-c++
    . /opt/rh/gcc-toolset-15/enable
fi
if command -v apk >/dev/null; then
    apk add --no-cache bash
    apk add --no-cache --upgrade --repository=https://dl-cdn.alpinelinux.org/alpine/v3.23/main g++ binutils
fi

exec bash python/scripts/build_wheel.sh
