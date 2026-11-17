"""The ABI contract: one normalized line per declaration of include/portopt/portopt.h. Standard library only.

    python scripts/abi_contract.py write                         regenerate include/portopt/portopt.abi
    python scripts/abi_contract.py check [--base FILE]           fail if the header and the contract differ;
                                                                 with --base, also if a line of FILE is gone
    python scripts/abi_contract.py base-ref                      the latest stable tag v1.x or later, or nothing
    python scripts/abi_contract.py def --library NAME OUT        the Windows export list of the contract
    python scripts/abi_contract.py exports LISTING               fail if the exports of a DLL (`objdump -p DLL`)
                                                                 differ from the contract

A line is `define NAME value`, `typedef ...` or `function RETURN name(PARAMETER TYPES)`. Comments, parameter names and
the PO_API platform macro are not part of it, so renaming a parameter is not an ABI change. Until 1.0 a line may
change (PO_ABI_VERSION is 0); from 1.0 the contract only grows: removing or changing a line is what `check --base`
rejects.
"""

import argparse
import difflib
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
HEADER = ROOT / "include" / "portopt" / "portopt.h"
CONTRACT = ROOT / "include" / "portopt" / "portopt.abi"

COMMENTS = re.compile(r"/\*.*?\*/|//[^\n]*", re.DOTALL)
DEFINE = re.compile(r"^#define\s+(PO_\w+)[ \t]+(\S.*)$", re.MULTILINE)
TYPEDEF = re.compile(r"\btypedef\b([^;{}]*);")
FUNCTION = re.compile(r"\bPO_API\b([^;{}]*?)\b(po_\w+)\s*\(([^)]*)\)\s*;")
TOKEN = re.compile(r"\w+|\*")
PLATFORM_MACROS = {"PO_API", "PO_BUILDING"}
NOT_A_NAME = {"const", "void", "struct"}


def join_tokens(tokens: list[str]) -> str:
    """Tokens as C text: a pointer star sticks to the type before it."""
    return " ".join(tokens).replace(" *", "*")


def parameter_type(parameter: str) -> str:
    """The type of one parameter, without its name."""
    tokens = TOKEN.findall(parameter)
    if len(tokens) >= 2 and tokens[-1] not in NOT_A_NAME and tokens[-1] != "*" and tokens[-2] != "struct":
        tokens = tokens[:-1]
    return join_tokens(tokens)


def function_line(prefix: str, name: str, parameters: str) -> str:
    returned = join_tokens(TOKEN.findall(prefix))
    types = [parameter_type(part) for part in parameters.split(",")]
    return f"function {returned} {name}({', '.join(types)})"


def refuse_unreadable(text: str) -> None:
    """Raises for anything the contract format cannot hold: it must never be silently left out of the contract."""
    code = "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("#"))
    for statement in re.findall(r"\bPO_API\b[^;{}]*;", code):
        if not FUNCTION.fullmatch(statement.strip()):
            raise ValueError(f"cannot read this declaration into the contract: {' '.join(statement.split())}")
    if re.search(r"\benum\b", code):
        raise ValueError("enum is not part of the contract format: use #define values")
    if re.search(r"\btypedef\b[^;]*\{", code):
        raise ValueError("a typedef struct with a body is not part of the contract format: use an opaque struct")


def contract_lines(header_text: str) -> list[str]:
    """The contract of a header, in header order."""
    text = COMMENTS.sub(" ", header_text)
    refuse_unreadable(text)
    found: list[tuple[int, str]] = []
    for match in DEFINE.finditer(text):
        if match.group(1) not in PLATFORM_MACROS:
            found.append((match.start(), f"define {match.group(1)} {' '.join(match.group(2).split())}"))
    for match in TYPEDEF.finditer(text):
        found.append((match.start(), f"typedef {' '.join(match.group(1).split())}"))
    for match in FUNCTION.finditer(text):
        found.append((match.start(), function_line(match.group(1), match.group(2), match.group(3))))
    return [line for _, line in sorted(found)]


def removed_lines(base: list[str], current: list[str]) -> list[str]:
    """The lines of `base` that `current` no longer has: what a release may not take away."""
    kept = set(current)
    return [line for line in base if line not in kept]


def function_names(lines: list[str]) -> list[str]:
    names = []
    for line in lines:
        match = re.match(r"function .*?\b(po_\w+)\(", line)
        if match:
            names.append(match.group(1))
    return names


def def_file(lines: list[str], library: str) -> str:
    return "\n".join(["LIBRARY " + library, "EXPORTS", *function_names(lines)]) + "\n"


# A row of the "[Ordinal/Name Pointer] Table" of `objdump -p`: "[   0] +base[   1]  0000 po_abi_version".
EXPORT_ROW = re.compile(r"^\s*\[\s*\d+\]\s+\+base\[\s*\d+\]\s+[0-9a-fA-F]+\s+(\S+)\s*$", re.MULTILINE)


def export_differences(lines: list[str], objdump: str) -> list[str]:
    """How the export table of a DLL (`objdump -p`) differs from the contract's functions."""
    contract = set(function_names(lines))
    exported = set(EXPORT_ROW.findall(objdump))
    problems = [f"exported but not in the contract: {name}" for name in sorted(exported - contract)]
    problems += [f"in the contract but not exported: {name}" for name in sorted(contract - exported)]
    return problems


def latest_stable(tags: list[str]) -> str | None:
    """The highest tag vN.M.P with N >= 1 and no pre-release suffix."""
    versions = []
    for tag in tags:
        match = re.fullmatch(r"v(\d+)\.(\d+)\.(\d+)", tag)
        if match and int(match.group(1)) >= 1:
            versions.append((tuple(int(part) for part in match.groups()), tag))
    if not versions:
        return None
    return max(versions)[1]


def read_lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").splitlines()


def check(header: Path, contract: Path, base: Path | None) -> int:
    expected = contract_lines(header.read_text(encoding="utf-8"))
    committed = read_lines(contract)
    if expected != committed:
        diff = difflib.unified_diff(committed, expected, "portopt.abi", "portopt.h", lineterm="")
        print("\n".join(diff))
        print("error: portopt.h and portopt.abi differ; if the change is intended run: "
              "python scripts/abi_contract.py write")
        return 1
    gone = removed_lines(read_lines(base), committed) if base else []
    for line in gone:
        print(f"error: the published contract has a line this one does not: {line}")
    return 1 if gone else 0


def stable_tag() -> str | None:
    tags = subprocess.run(["git", "tag", "--list", "v*"], cwd=ROOT, capture_output=True, text=True, check=True)
    return latest_stable(tags.stdout.split())


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["write", "check", "base-ref", "def", "exports"])
    parser.add_argument("target", nargs="?", type=Path, help="def: the output file; exports: the output of objdump -p")
    parser.add_argument("--header", type=Path, default=HEADER)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    parser.add_argument("--base", type=Path)
    parser.add_argument("--library", default="libportopt.dll")
    args = parser.parse_args(argv)
    if args.command == "write":
        lines = contract_lines(args.header.read_text(encoding="utf-8"))
        args.contract.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
        return 0
    if args.command == "check":
        return check(args.header, args.contract, args.base)
    if args.command == "base-ref":
        print(stable_tag() or "")
        return 0
    if args.command == "def":
        args.target.write_text(def_file(read_lines(args.contract), args.library), encoding="utf-8", newline="\n")
        return 0
    problems = export_differences(read_lines(args.contract), args.target.read_text(encoding="utf-8"))
    print("\n".join(problems))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
