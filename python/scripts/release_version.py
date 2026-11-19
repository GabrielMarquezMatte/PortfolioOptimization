"""The versions of a release tag, and the check that the repository has one version.

    python python/scripts/release_version.py TAG          the PEP 440 version of the wheel, then the x.y.z of CMake
    python python/scripts/release_version.py check [TAG]  fail unless PO_VERSION of CMakeLists.txt equals the version
                                                           of python/pyproject.toml; with TAG, also unless the tag is
                                                           that version and CHANGELOG.md has its dated section
    python python/scripts/release_version.py prerelease TAG  "true" for a prerelease tag (v0.1.0rc1 too), else "false"
    python python/scripts/release_version.py rehearsal       a prerelease tag of the repository's version (v0.1.0-rc.0),
                                                           which a dry run builds with so it takes the path of a tag

A release candidate is built from the tag alone (v0.1.0-rc.1 gives the wheel 0.1.0rc1 and CMake 0.1.0), so no commit
changes a version for it.
"""

import re
import sys
from pathlib import Path

from packaging.version import Version

ROOT = Path(__file__).resolve().parents[2]
CMAKE_VERSION = re.compile(r"^\s*set\(PO_VERSION ([^\s)]+)\)", re.MULTILINE)
PYPROJECT_VERSION = re.compile(r'^version = "([^"]+)"', re.MULTILINE)


def parse(tag: str) -> Version:
    if not tag.startswith("v"):
        raise ValueError(f"tag {tag!r} must start with 'v'")
    return Version(tag[1:])


def python_version(tag: str) -> str:
    return str(parse(tag))


def cmake_version(tag: str) -> str:
    return parse(tag).base_version


def is_prerelease(tag: str) -> bool:
    return parse(tag).is_prerelease


def rehearsal_tag(cmake: str) -> str:
    return f"v{found(CMAKE_VERSION, cmake, 'CMakeLists.txt')}-rc.0"


def found(pattern: re.Pattern[str], text: str, where: str) -> str:
    match = pattern.search(text)
    if match is None:
        raise ValueError(f"no version in {where}")
    return match.group(1)


def mismatches(cmake: str, pyproject: str, changelog: str, tag: str | None) -> list[str]:
    """How the texts of CMakeLists.txt, python/pyproject.toml and CHANGELOG.md disagree with each other and the tag."""
    cmake_text = found(CMAKE_VERSION, cmake, "CMakeLists.txt")
    pyproject_text = found(PYPROJECT_VERSION, pyproject, "python/pyproject.toml")
    problems = []
    if cmake_text != pyproject_text:
        problems.append(
            f"CMakeLists.txt has PO_VERSION {cmake_text} but python/pyproject.toml has version {pyproject_text}"
        )
    if tag is None:
        return problems
    if cmake_version(tag) != cmake_text:
        problems.append(f"tag {tag} is version {cmake_version(tag)} but CMakeLists.txt has PO_VERSION {cmake_text}")
    section = re.compile(rf"^## \[{re.escape(cmake_text)}\] - \d{{4}}-\d{{2}}-\d{{2}}$", re.MULTILINE)
    if section.search(changelog) is None:
        problems.append(f"CHANGELOG.md has no section '## [{cmake_text}] - YYYY-MM-DD'")
    return problems


def main(argv: list[str]) -> int:
    if argv[:1] == ["prerelease"]:
        print("true" if is_prerelease(argv[1]) else "false")
        return 0
    if argv[:1] == ["rehearsal"]:
        print(rehearsal_tag((ROOT / "CMakeLists.txt").read_text(encoding="utf-8")))
        return 0
    if argv[:1] != ["check"]:
        print(python_version(argv[0]))
        print(cmake_version(argv[0]))
        return 0
    tag = argv[1] if len(argv) > 1 else None
    files = ("CMakeLists.txt", "python/pyproject.toml", "CHANGELOG.md")
    cmake, pyproject, changelog = ((ROOT / name).read_text(encoding="utf-8") for name in files)
    problems = mismatches(cmake, pyproject, changelog, tag)
    for problem in problems:
        print(f"error: {problem}", file=sys.stderr)
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
