#!/usr/bin/env bash
# Writes the ABI of libportopt.so as libabigail XML: scripts/abi_dump.sh OUT.xml
#
# The release attaches this file for its Linux x64 library (portopt-abi-linux-x64.xml): the baseline an abidiff check
# compares against once the ABI is frozen at 1.0. It is made with GCC 15 on ubuntu-24.04 (scripts/install_abi_tools.sh),
# so a later dump with the same toolchain differs only where the code does. The library is built
# with -g because abidw reads the types from the debug information; --headers-dir keeps only what portopt.h declares.
#
# Environment: CC and CXX name the compiler (default gcc-15 and g++-15).
set -euo pipefail
out=$(realpath -m "${1:?usage: abi_dump.sh OUT.xml}")
cd "$(dirname "$0")/.."

cmake -S . -B build/abi -G Ninja -DCMAKE_BUILD_TYPE=Release -DPO_BUILD_TESTS=OFF \
    "-DCMAKE_C_COMPILER=${CC:-gcc-15}" "-DCMAKE_CXX_COMPILER=${CXX:-g++-15}" \
    -DCMAKE_C_FLAGS=-g -DCMAKE_CXX_FLAGS=-g
cmake --build build/abi --target portopt
mkdir -p "$(dirname "$out")"
abidw --headers-dir include/portopt --drop-private-types --out-file "$out" build/abi/libportopt.so
grep -c 'po_solve' "$out" > /dev/null
echo "wrote $out"
