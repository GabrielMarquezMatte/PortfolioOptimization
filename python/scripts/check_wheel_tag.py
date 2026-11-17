"""Fails unless every wheel in the given directory installs on this machine and, with an expected tag, does not
claim a newer system than that tag.

Run after building a wheel:  python python/scripts/check_wheel_tag.py dist [EXPECTED_TAG]

A glibc wheel may carry a lower manylinux minor than expected (the library only uses older symbols, so it runs on
more systems); a higher one, another architecture or another platform fails.
"""

import re
import sys
from pathlib import Path

from packaging.tags import sys_tags
from packaging.utils import parse_wheel_filename

MANYLINUX = re.compile(r"manylinux_2_(\d+)_(.+)")


def tag_allowed(actual: str, expected: str) -> bool:
    if actual == expected:
        return True
    actual_match, expected_match = MANYLINUX.fullmatch(actual), MANYLINUX.fullmatch(expected)
    if actual_match is None or expected_match is None:
        return False
    return actual_match[2] == expected_match[2] and int(actual_match[1]) <= int(expected_match[1])


def main() -> int:
    supported = set(sys_tags())
    expected = sys.argv[2] if len(sys.argv) > 2 else None
    wheels = sorted(Path(sys.argv[1]).glob("*.whl"))
    if not wheels:
        print(f"error: no wheel in {sys.argv[1]}", file=sys.stderr)
        return 1
    for wheel in wheels:
        tags = parse_wheel_filename(wheel.name)[3]
        if tags.isdisjoint(supported):
            print(f"error: {wheel.name} does not install on this machine", file=sys.stderr)
            return 1
        if expected is not None and not all(tag_allowed(tag.platform, expected) for tag in tags):
            print(f"error: {wheel.name} is not tagged {expected}: the bundled library needs a newer system "
                  "than the tag claims", file=sys.stderr)
            return 1
        print(f"ok: {wheel.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
