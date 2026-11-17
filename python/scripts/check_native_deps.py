"""Fails if the native library imports a compiler runtime that a clean machine does not have, or exports anything
besides the C interface (po_*): a static libstdc++ must not leak its symbols to the program that loads the library.

Run after building the library:  python python/scripts/check_native_deps.py PATH
"""

import platform
import subprocess
import sys

WINDOWS_FORBIDDEN = ("libstdc++", "libgcc", "libwinpthread")
LINUX_FORBIDDEN = ("libstdc++", "libgcc_s")
MACOS_ALLOWED = ("/usr/lib/", "/System/")
# Symbols the toolchain adds to the dynamic table of a shared library (musl, aarch64): not ours.
# ELF symbol types that name no function or object (aarch64 lists a section symbol such as .data as GLOBAL).
NAMELESS_TYPES = {"SECTION", "FILE"}
TOOLCHAIN_SYMBOLS = {"_init", "_fini", ".init", ".fini", "__bss_start", "_edata", "_end"}


def forbidden(system: str, needed: list[str]) -> list[str]:
    if system == "Windows":
        return [name for name in needed if name.lower().startswith(WINDOWS_FORBIDDEN)]
    if system == "Linux":
        return [name for name in needed if name.startswith(LINUX_FORBIDDEN)]
    return [name for name in needed if not name.startswith(MACOS_ALLOWED)]


def foreign_exports(system: str, lines: list[str]) -> list[str]:
    """The defined dynamic symbols of `readelf --dyn-syms --wide` (Linux) or `nm -gU` (macOS) that are not po_*."""
    names = []
    for line in lines:
        fields = line.split()
        if system == "Linux":
            # Num: Value Size Type Bind Vis Ndx Name; a symbol with Ndx UND is an import.
            if len(fields) < 8 or not fields[0].endswith(":") or fields[6] in ("UND", "Ndx") or fields[3] in NAMELESS_TYPES:
                continue
            names.append(fields[7].split("@")[0])
        elif len(fields) == 3:
            names.append(fields[2].removeprefix("_"))
    return [name for name in names if not name.startswith("po_") and name not in TOOLCHAIN_SYMBOLS]


def exports(system: str, path: str) -> list[str]:
    if system == "Linux":
        return foreign_exports(system, run("readelf", "--dyn-syms", "--wide", path))
    if system == "Darwin":
        return foreign_exports(system, run("nm", "-gU", path))
    return []  # Windows: scripts/abi_contract.py exports compares the export table of the DLL with the contract


def run(*command: str) -> list[str]:
    return subprocess.run(command, check=True, capture_output=True, text=True).stdout.splitlines()


def needed(system: str, path: str) -> list[str]:
    if system == "Windows":
        return [line.split("DLL Name:")[1].strip() for line in run("objdump", "-p", path) if "DLL Name:" in line]
    if system == "Linux":
        return [
            line.split("[")[1].rstrip("]")
            for line in run("readelf", "-d", path)
            if "(NEEDED)" in line
        ]
    # otool -L: the file name, then the library's own id, then what it imports.
    return [line.split()[0] for line in run("otool", "-L", path)[2:]]


def main() -> int:
    system = platform.system()
    dependencies = needed(system, sys.argv[1])
    bad = forbidden(system, dependencies)
    if bad:
        print(f"error: {sys.argv[1]} needs {', '.join(bad)}", file=sys.stderr)
        return 1
    foreign = exports(system, sys.argv[1])
    if foreign:
        print(f"error: {sys.argv[1]} exports {len(foreign)} symbols that are not po_*: {', '.join(foreign[:10])}", file=sys.stderr)
        return 1
    print(f"ok: {', '.join(dependencies)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
