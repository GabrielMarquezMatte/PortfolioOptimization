# Versioning and compatibility

What a new version may change, and how to write code that keeps compiling and running across minor versions. PortOpt
follows [Semantic Versioning](https://semver.org/) from 1.0: a major version (2.0) may break your code, while a minor
(1.1) or a patch (1.0.1) does not. The C library, the Python package, the C++ wrapper and the .NET package share one
version number.

## What a minor version may add

- New functions of the C API, among them new setters (`po_spec_set_*`, `po_data_set_*`). New configuration is always
  a new function, never a new parameter of an existing one.
- New values of `PO_OBJECTIVE_*` and `PO_CONSTRAINT_*`, and the binding enums and result fields that go with them.
- New fields of the binding results (Python `Result`, C++ `po::Result`, .NET `Result`) and of the binding specs
  (Python `Spec`, C++ `po::Spec`, .NET `Spec`), appended after the existing ones. A Python `Spec` field is
  keyword-only; a C++ `po::Spec` is written with designated initializers, which a new member does not disturb.
- New properties of .NET `MarketInputs` and `Spec`, each defaulting to the earlier behaviour.
- Bug fixes. A problem that used to be solved wrong, or reported with the wrong status, can come out differently after
  a fix.

A minor version does not remove or rename anything, change a signature or the meaning of an existing function or
value, or renumber a value.

## Values the library returns

`po_status` and `PO_SOLVE_*` are returned to you, so a new value is something an older caller could meet. Two rules
keep that from happening by surprise:

- A new `PO_SOLVE_*` value is returned only to a request that uses something added with it (a setter or an objective
  that an older caller cannot have set). A status that any existing request could receive is a major version.
- A new `po_status` code is a failure like the others: every binding turns an unknown code into its generic error
  (Python `PortOptError`, C++ `po::Error`, .NET `PortOptException`), with the native message.

A `switch` over an enum of PortOpt should still have a default branch.

## Error messages

`po_last_error()` and the messages of the binding exceptions are for people. They may change in any version; branch
on the status (`po_status`, the Python exception type, `po::Error::status`, `PortOptException.Status`), not on the
text.

## The ABI

The C ABI is versioned by `PO_ABI_VERSION`, which `po_abi_version()` reports. It changes only with a major version,
when a signature or a meaning does; a new function leaves it as it is. The Python package, the C++ wrapper and the
.NET package check the exact ABI version when they load the library, so a library of another major version is refused
at load, with a message, rather than failing later. A binding of 1.y needs a library of at least 1.y: the functions it
calls may be newer than an older 1.x library has, and loading such a library fails on the missing export.

The CMake package `PortOpt` of a 1.y release satisfies `find_package(PortOpt 1.x)` for every x up to y.

## How it is enforced

- `include/portopt/portopt.abi` is the contract written out, one line per declaration of `portopt.h`. `ctest` fails
  when the header and the contract differ, and from 1.0 `scripts/abi_contract.py check --base` fails when a line of
  the last stable release is gone or changed.
- Each release publishes the libabigail dump of the Linux library, the baseline of a binary comparison.
- The .NET package tracks its public surface in `PublicAPI.Shipped.txt` and `PublicAPI.Unshipped.txt`, and .NET
  package validation compares it with the last release.
