"""python/scripts/release_version.py: the versions of a release tag and the check that the repository agrees."""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "python" / "scripts" / "release_version.py"
spec = importlib.util.spec_from_file_location("release_version", SCRIPT)
assert spec is not None and spec.loader is not None
release_version = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release_version)

CMAKE = 'set(PO_VERSION 0.1.0 CACHE STRING "PortOpt version (x.y.z)")\n'
PYPROJECT = '[project]\nname = "pyportopt"\nversion = "0.1.0"\n'
UNRELEASED = "# Changelog\n\n## [Unreleased]\n\n### Added\n"
RELEASED = "# Changelog\n\n## [Unreleased]\n\n## [0.1.0] - 2026-11-20\n\n### Added\n"


def test_plain_tag() -> None:
    assert release_version.python_version("v0.1.0") == "0.1.0"
    assert release_version.cmake_version("v0.1.0") == "0.1.0"


def test_prerelease_tag() -> None:
    for tag in ("v0.1.0-rc.1", "v0.1.0rc1"):
        assert release_version.python_version(tag) == "0.1.0rc1"
        assert release_version.cmake_version(tag) == "0.1.0"


def test_tag_without_v_is_rejected() -> None:
    with pytest.raises(ValueError):
        release_version.python_version("0.1.0")


def test_cmake_and_pyproject_must_agree() -> None:
    assert release_version.mismatches(CMAKE, PYPROJECT, UNRELEASED, None) == []
    problems = release_version.mismatches(CMAKE, PYPROJECT.replace("0.1.0", "0.2.0"), UNRELEASED, None)
    assert problems == ["CMakeLists.txt has PO_VERSION 0.1.0 but python/pyproject.toml has version 0.2.0"]


def test_tag_must_match_the_cmake_version() -> None:
    problems = release_version.mismatches(CMAKE, PYPROJECT, RELEASED, "v0.2.0")
    assert problems == ["tag v0.2.0 is version 0.2.0 but CMakeLists.txt has PO_VERSION 0.1.0"]


def test_tag_needs_a_dated_changelog_section() -> None:
    problems = release_version.mismatches(CMAKE, PYPROJECT, UNRELEASED, "v0.1.0-rc.1")
    assert problems == ["CHANGELOG.md has no section '## [0.1.0] - YYYY-MM-DD'"]
    assert release_version.mismatches(CMAKE, PYPROJECT, RELEASED, "v0.1.0-rc.1") == []
    assert release_version.mismatches(CMAKE, PYPROJECT, RELEASED, "v0.1.0") == []


def test_real_repository_agrees() -> None:
    assert release_version.main(["check"]) == 0
