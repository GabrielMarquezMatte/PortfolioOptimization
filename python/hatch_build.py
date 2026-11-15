"""Forces a platform-specific wheel tag: the wheel bundles the native library for one OS and CPU (see _lib/), so it
is not pure Python even though every .py file in it is. The package has no CPython extension (ctypes calls the shared
library), so it is interpreter-independent: only the platform component of the tag is specific.

The platform component comes from the bundled library, not from the build machine: a runner's own tag
(macosx_26_0, manylinux_2_39) would shut out every older system the library loads on."""

import platform
import re
import struct
from pathlib import Path

from hatchling.builders.hooks.plugin.interface import BuildHookInterface
from packaging.tags import platform_tags

LIB_DIR = Path(__file__).parent / "src" / "portopt" / "_lib"
LIBRARIES = {"Linux": "libportopt.so", "Darwin": "libportopt.dylib", "Windows": "libportopt.dll"}
LC_BUILD_VERSION = 0x32


def glibc_floor(library: bytes) -> str:
    """The newest GLIBC_2.N symbol version the library needs, as manylinux's `2_N`."""
    minors = [int(minor) for minor in re.findall(rb"GLIBC_2\.(\d+)", library)]
    return f"2_{max(minors, default=17)}"


def macos_floor(library: bytes) -> str:
    """The minimum macOS a thin 64-bit Mach-O declares in LC_BUILD_VERSION, as `MAJOR_0`."""
    commands = struct.unpack_from("<I", library, 16)[0]
    offset = 32
    for _ in range(commands):
        command, size = struct.unpack_from("<II", library, offset)
        if command == LC_BUILD_VERSION:
            minos = struct.unpack_from("<I", library, offset + 12)[0]
            return f"{minos >> 16}_0"
        offset += size
    raise ValueError("the native library has no LC_BUILD_VERSION load command")


def platform_tag(library: bytes) -> str:
    machine = platform.machine().lower()
    system = platform.system()
    if system == "Linux":
        if b"libc.musl-" in library:
            return next(tag for tag in platform_tags() if tag.startswith("musllinux_"))
        return f"manylinux_{glibc_floor(library)}_{machine}"
    if system == "Darwin":
        return f"macosx_{macos_floor(library)}_{machine}"
    return next(iter(platform_tags()))


class NativeLibraryBuildHook(BuildHookInterface):
    def initialize(self, version, build_data):
        build_data["pure_python"] = False
        build_data["infer_tag"] = True
        if version == "editable":
            # An editable install loads the library from PORTOPT_LIB or _lib at run time; _lib may be empty.
            build_data["tag"] = f"py3-none-{next(iter(platform_tags()))}"
            return
        library = LIB_DIR / LIBRARIES[platform.system()]
        # Only the library goes in the wheel: the import library and .def next to it are for the release assets.
        build_data["force_include"][str(library)] = f"portopt/_lib/{library.name}"
        build_data["tag"] = f"py3-none-{platform_tag(library.read_bytes())}"
