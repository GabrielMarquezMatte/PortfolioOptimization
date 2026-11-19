#!/usr/bin/env bash
# Checks a PortOpt package as a user gets it: restores csharp/tests/PortOpt.Consumer from the given package folder
# (and nuget.org), runs it without PORTOPT_LIB, then publishes it with Native AOT for the host RID and runs the
# executable. Usage: check_package.sh <folder with the .nupkg> <version>
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "usage: $0 <package folder> <version>" >&2
    exit 2
fi

# A path dotnet understands, also from Git Bash on Windows.
native_path() {
    (cd "$1" && (pwd -W 2>/dev/null || pwd))
}

packages="$(native_path "$1")"
version="$2"
consumer="$(native_path "$(dirname "$0")/../tests/PortOpt.Consumer")"
work="$(native_path "$(mktemp -d)")"
trap 'rm -rf "$work"' EXIT

# The package must find its own library; a fresh package cache keeps a repacked version from being served stale.
unset PORTOPT_LIB
export NUGET_PACKAGES="$work/packages"
properties=(-p:RestoreAdditionalProjectSources="$packages" -p:PortOptVersion="$version")

dotnet run --project "$consumer" -c Release "${properties[@]}"

rid="$(dotnet --info | sed -n 's/^ *RID: *//p' | head -n 1 | tr -d '\r')"
dotnet publish "$consumer" -c Release -r "$rid" -p:PublishAot=true "${properties[@]}" -o "$work/aot"
"$work/aot/PortOpt.Consumer"
