#!/usr/bin/env bash
# Checks that a downstream project can use the C++ wrapper both through add_subdirectory and through find_package
# after `cmake --install`, and that find_package refuses the next minor version (0.x breaks the ABI between minor
# versions). cpp.yml runs it on every platform.
#
# Environment:
#   CPP_CMAKE_ARGS  extra configure arguments, for example -DCMAKE_CXX_COMPILER=cl
set -euo pipefail
cd "$(dirname "$0")/../.."

version=$(sed -n 's/^ *set(PO_VERSION \([0-9.]*\)).*/\1/p' CMakeLists.txt)
major_minor=${version%.*}
next_minor="${major_minor%%.*}.$((${major_minor#*.} + 1))"

# shellcheck disable=SC2206
args=(${CPP_CMAKE_ARGS:-} -G Ninja -DCMAKE_BUILD_TYPE=Release)

cmake -S . -B build/cpp-install "${args[@]}" -DPO_BUILD_TESTS=OFF
cmake --build build/cpp-install
rm -rf build/cpp-prefix
cmake --install build/cpp-install --prefix "$PWD/build/cpp-prefix"

consumer=(cmake -S cpp/tests/package "${args[@]}")
"${consumer[@]}" -B build/consumer-subdirectory -DPORTOPT_CONSUMER_MODE=add_subdirectory
cmake --build build/consumer-subdirectory
ctest --test-dir build/consumer-subdirectory --output-on-failure

"${consumer[@]}" -B build/consumer-find -DPORTOPT_CONSUMER_MODE=find_package \
    "-DPORTOPT_CONSUMER_VERSION=${major_minor}" "-DCMAKE_PREFIX_PATH=$PWD/build/cpp-prefix"
cmake --build build/consumer-find
ctest --test-dir build/consumer-find --output-on-failure

rm -rf build/consumer-other
if "${consumer[@]}" -B build/consumer-other -DPORTOPT_CONSUMER_MODE=find_package \
    "-DPORTOPT_CONSUMER_VERSION=${next_minor}" "-DCMAKE_PREFIX_PATH=$PWD/build/cpp-prefix" > /dev/null 2>&1; then
    echo "error: find_package(PortOpt ${next_minor}) accepted the installed ${version}"
    exit 1
fi
echo "ok: find_package(PortOpt ${next_minor}) refuses the installed ${version}"
