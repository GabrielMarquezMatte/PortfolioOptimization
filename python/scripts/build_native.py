"""Builds the portopt shared library and copies it into the package's _lib directory.

Usage:  python python/scripts/build_native.py [--build-dir DIR] [extra cmake configure arguments]
"""

import argparse
import os
import platform
import shutil
import subprocess
import sys
from collections.abc import Mapping
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LIB_DIR = ROOT / "python" / "src" / "portopt" / "_lib"
ISA_FLAGS = ("-march", "-mcpu", "-mavx", "-mtune")
LIBRARIES = ("libportopt.so", "libportopt.dylib", "libportopt.dll", "libportopt.dll.a")
MINIMUM_GCC = 15


def check_flags(environ: Mapping[str, str]) -> None:
    """The distributed binaries are generic (no -march), like gcc-bench; a baked-in ISA flag breaks older CPUs."""
    for name in ("CFLAGS", "CXXFLAGS"):
        for flag in environ.get(name, "").split():
            if flag.startswith(ISA_FLAGS):
                sys.exit(f"error: {name} contains {flag}; the wheel must run on the baseline CPU")


def check_gcc(version: str) -> None:
    """The project needs GCC 15 or newer; an older one in a container must stop the build, not ship."""
    if int(version.split(".")[0]) < MINIMUM_GCC:
        sys.exit(f"error: GCC {MINIMUM_GCC} or newer is required, found {version}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", type=Path, default=ROOT / "build" / "wheel-native")
    args, cmake_args = parser.parse_known_args()
    check_flags(os.environ)

    configure = ["cmake", "-S", str(ROOT), "-B", str(args.build_dir), "-G", "Ninja",
                 "-DCMAKE_BUILD_TYPE=Release", "-DPO_BUILD_TESTS=OFF"]
    if platform.system() in ("Linux", "Windows"):
        version = subprocess.run(["gcc", "-dumpfullversion"], check=True, capture_output=True, text=True).stdout
        check_gcc(version.strip())
        configure += ["-DCMAKE_C_COMPILER=gcc", "-DCMAKE_CXX_COMPILER=g++"]
    subprocess.run([*configure, *cmake_args], check=True)
    subprocess.run(["cmake", "--build", str(args.build_dir), "--target", "portopt"], check=True)

    LIB_DIR.mkdir(parents=True, exist_ok=True)
    for old in LIB_DIR.iterdir():
        old.unlink()
    for name in LIBRARIES:
        if (built := args.build_dir / name).exists():
            shutil.copy2(built, LIB_DIR / name)
            print(f"copied {name}")


if __name__ == "__main__":
    main()
