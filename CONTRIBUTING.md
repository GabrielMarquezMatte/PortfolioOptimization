# Contributing

Small, focused pull requests are preferred over large rewrites. Read [ARCHITECTURE.md](ARCHITECTURE.md) for the
shape of the codebase and [STYLEGUIDE.md](STYLEGUIDE.md) for the code style.

## Ground rules

- The C API in `include/portopt/portopt.h` is the only stable contract. Everything under `src/` can change
  freely. The stability rules are in [ARCHITECTURE.md](ARCHITECTURE.md#stability-rules-of-the-c-api).
- Warnings are errors (`-Wall -Wextra -Wpedantic -Werror`, plus `-Wnrvo` on GCC 14+; `/W4 /WX` on MSVC).
- Commit messages carry no `Co-Authored-By` trailer.
- Do not copy or translate `quadprog` code (GPL). The Goldfarb-Idnani method comes from the paper (Goldfarb and
  Idnani, *Mathematical Programming* 27, 1983). Reading `quadprog` to understand a convention is fine.

## Build and test (C++)

Presets are in `CMakePresets.json`:

```
cmake --preset gcc-release
cmake --build --preset gcc-release
ctest --preset gcc-release
```

`gcc-debug` and `msvc-release` work the same way. `native-release` uses the default compiler (`CC` and `CXX`)
and is what CI uses on Linux and macOS. GCC 15 or newer is the supported GCC.
The Linux wheels are built with GCC 15 too: `gcc-toolset-15` (15.2.1) from AlmaLinux 8 in `manylinux_2_28` and
the `g++` of Alpine 3.23 in `musllinux_1_2` (`python/scripts/build_in_container.sh`); the macOS wheels declare
macOS 11.0.

`ctest` runs the GoogleTest suite, a C11 compile of the public header and the `style` test.

### Windows and MinGW

The MinGW build links `libstdc++`, `libgcc` and `winpthread` statically, so a Windows Python can load
`libportopt.dll` without the MSYS2 `bin` directory on `PATH`. Keep it that way when touching the link options.

## Python tests

The Python package is managed with [uv](https://docs.astral.sh/uv/). `uv sync --project python` creates
`python/.venv` with the `dev` group (pytest, mypy).

Point the package at the library you built with `PORTOPT_LIB`, then run pytest and mypy:

```
PORTOPT_LIB=$PWD/build/gcc-release/libportopt.dll uv run --project python pytest python/tests -q
PORTOPT_LIB=$PWD/build/gcc-release/libportopt.dll uv run --project python mypy --config-file python/pyproject.toml python/src python/tests
```

mypy runs in strict mode (`[tool.mypy]` in `python/pyproject.toml`); from the repository root it finds that file only
through `--config-file`.

The library is `libportopt.dll` (MinGW), `Release/portopt.dll` (MSVC), `libportopt.so` or `libportopt.dylib`.
Without `PORTOPT_LIB` the package looks in `python/src/portopt/_lib`. When it finds no library the tests are
skipped and pytest exits with status 5.

## Style checks

- `python tests/style/check_style.py` runs on its own and as the `style` ctest. It fails on `else`, on folders
  with more than 8 modules and, for Python, on functions above 60 lines, more than 5 parameters or nesting above 3.
- C++ function size, parameter count and nesting are checked by clang-tidy with the repository `.clang-tidy`.
  CI runs it over `src/`; locally, after configuring with `-DCMAKE_EXPORT_COMPILE_COMMANDS=ON`:

  ```
  clang-tidy -p build/gcc-release src/solve/solve.cpp
  ```

## Making a release

The version lives in `CMakeLists.txt` (`PO_VERSION`) and `python/pyproject.toml`; `python
python/scripts/release_version.py check` fails when they differ, and `check vX.Y.Z` also when the tag is another
version or `CHANGELOG.md` has no dated `## [X.Y.Z] - YYYY-MM-DD` section.

1. Rename `## [Unreleased]` in `CHANGELOG.md` to the version with today's date, keep an empty `## [Unreleased]`
   above it, set the version in both files, and push; `ci.yml`, `cpp.yml` and `wheels.yml` must be green.
2. Optionally run `release.yml` by hand (`dry_run` on by default): it builds and tests everything and publishes
   nothing.
3. Push a release candidate tag, `vX.Y.Z-rc.N`: the GitHub Release is marked as a prerelease and the wheels go to
   PyPI as the prerelease `X.Y.ZrcN` (derived from the tag; no commit changes the version for a candidate), which
   pip installs only when asked for (`pip install --pre pyportopt` or `pyportopt==X.Y.ZrcN`). Install it in a clean
   environment to check it.
4. Push `vX.Y.Z`: the GitHub Release and PyPI.

PyPI publishes by trusted publishing: the project names this repository, the workflow `release.yml` and the
environment `pypi`, so no token is stored.

## Documentation

Repository documents are in English. Measurements and justifications go in `docs/`, never in code comments
([STYLEGUIDE.md](STYLEGUIDE.md) says why). A numerical decision, such as a tolerance, goes in
[docs/guide/validation.md](docs/guide/validation.md) with the test that pins it. Add an entry under
`## [Unreleased]` in [CHANGELOG.md](CHANGELOG.md) for any user-visible change.
