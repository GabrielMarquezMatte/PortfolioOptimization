"""Fails if a PO_API function of portopt.h is not used by the C++ wrapper headers.

The wrapper is the second consumer of the C interface: a function it never calls is a function no test of it has
exercised. Run:  python cpp/scripts/check_abi_coverage.py [--header PATH] [--wrapper DIR]
"""

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NON_CODE = re.compile(r'//[^\n]*|/\*.*?\*/|"(?:\\.|[^"\\\n])*"', re.DOTALL)
DECLARATION = re.compile(r"\bPO_API\b[^;{]*?\b(po_\w+)\s*\(")
IDENTIFIER = re.compile(r"\bpo_\w+\b")

# Functions the wrapper chooses not to call, each with the reason.
UNUSED_ALLOWED: frozenset[str] = frozenset()


def strip_non_code(text: str) -> str:
    return NON_CODE.sub(" ", text)


def declared(header_text: str) -> set[str]:
    """The names of the PO_API functions."""
    lines = [line for line in strip_non_code(header_text).splitlines() if not line.lstrip().startswith("#")]
    return set(DECLARATION.findall("\n".join(lines)))


def used(wrapper_text: str) -> set[str]:
    """Every po_ identifier outside comments and strings: a call or a function passed as a pointer."""
    return set(IDENTIFIER.findall(strip_non_code(wrapper_text)))


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--header", type=Path, default=ROOT / "include" / "portopt" / "portopt.h")
    parser.add_argument("--wrapper", type=Path, default=ROOT / "cpp" / "include" / "portopt")
    args = parser.parse_args(argv)
    api = declared(args.header.read_text(encoding="utf-8"))
    text = "\n".join(path.read_text(encoding="utf-8") for path in sorted(args.wrapper.glob("*.hpp")))
    missing = sorted(api - used(text) - UNUSED_ALLOWED)
    for name in missing:
        print(f"error: {name} is declared with PO_API and not used by the wrapper")
    if missing:
        return 1
    print(f"ok: the wrapper uses {len(api) - len(UNUSED_ALLOWED)} of {len(api)} PO_API functions")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
