#!/usr/bin/env bash
# Installs what scripts/abi_dump.sh needs on ubuntu-24.04: GCC 15 (the compiler the project targets, from the toolchain
# PPA), CMake, Ninja and libabigail. Launchpad answers the PPA's signing key with a transient HTTP 500 at times, so
# adding the PPA is retried.
set -euo pipefail
sudo=""
if [ "$(id -u)" -ne 0 ]; then sudo=sudo; fi

$sudo apt-get -o Acquire::Retries=5 update -qq
$sudo apt-get -o Acquire::Retries=5 install -y -qq software-properties-common ca-certificates git cmake ninja-build abigail-tools
for attempt in 1 2 3 4 5; do
    $sudo add-apt-repository -y ppa:ubuntu-toolchain-r/test && break
    if [ "$attempt" -eq 5 ]; then exit 1; fi
    sleep $((attempt * 15))
done
$sudo apt-get -o Acquire::Retries=5 update -qq
$sudo apt-get -o Acquire::Retries=5 install -y -qq gcc-15 g++-15
