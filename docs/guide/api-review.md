# API review before 0.1.0

The C interface (`include/portopt/portopt.h`), the Python package `portopt` and the planned C++ wrapper
(`cpp/include/portopt/portopt.hpp`), reviewed before the first release. 0.x may still break (`PO_ABI_VERSION` is
0), so this is the cheap moment to correct names and errors; it freezes nothing.

Criteria: the stability rules of the C API (ARCHITECTURE.md: setters on opaque handles, `int32_t` enums with
`#define` values, `po_status` with a thread-local `po_last_error()`, independent handles, only `int32_t`,
`size_t`, `double` and pointers at the boundary); no surface without a use; the same name for the same thing in C,
Python and C++; error messages that name the argument; every function documented in the header; the uses of the
design (`Data`, `Spec`, `solve`, the metrics and the frontier loop).

Decisions: **keep**, **change**, **remove**, **add**. Each item has an id; a change lists its migration as
before → after. None of them touches `src/` beyond `src/capi`, so none changes a bit of a solution.

## C interface

| Id | Item | Finding | Decision | Breaks | Migration |
|---|---|---|---|---|---|
| C1 | `PO_ABI_VERSION`, `po_abi_version`, `po_version`, `po_last_error` | As the rules require; `po_version` is the library's `x.y.z`. | keep | — | — |
| C2 | `po_status`, `PO_OK` … `PO_INTERNAL_ERROR` | `PO_UNSUPPORTED` no longer occurs for a valid objective; kept as a code, since codes only grow and a future objective may need it. | keep | — | — |
| C3 | `po_spec`, `po_data`, `po_result` | Opaque, independent; `po_solve` copies spec and data. | keep | — | — |
| C4 | `PO_OBJECTIVE_*`, `PO_SOLVE_*`, `PO_CONSTRAINT_*` | Dense values from 0, each documented. | keep | — | — |
| C5 | `po_spec_create`, `po_spec_free` | — | keep | — | — |
| C6 | `po_spec_set_bounds`, `po_spec_set_asset_bounds` | `asset_bounds` takes one pair per *weight* (cash included), not per asset; the header says so. A rename to `weight_bounds` would be more exact but buys little: the Python and C++ layers choose the function from the shape of `bounds`. | keep | — | — |
| C7 | `po_spec_set_target_return`, `_risk_aversion`, `_importance`, `_max_volatility`, `_target_volatility`, `_risk_budget`, `_leverage` | One setter per field, each validated where it is set and documented with its domain; the field names match the Python `Spec` fields. There is no setter that clears a field: a new spec is the way back to a default, which every layer already does. | keep | — | — |
| C8 | `po_data_create`, `po_data_set_risk_free`, `po_data_set_returns`, `po_data_free` | — | keep | — | — |
| C9 | `po_data_asset_count` | The only function without a comment in the header. | change (doc) | no | add "The number of assets n given to po_data_create (cash not counted)." |
| C10 | `po_data_weight_count`, `po_result_weight_count` | — | keep | — | — |
| C11 | `po_solve` | Checks setter applicability, required data and the size of per-asset bounds, in that order; messages name the setter and the objective. | keep | — | — |
| C12 | `po_result_status`, `po_result_weights`, `po_result_multiplier`, `po_result_bound_multipliers`, `po_result_residual`, `po_result_free` | `po_result_status` writes an `int32_t` (a `PO_SOLVE_*`), as the enum rule requires. | keep | — | — |
| C13 | `po_result_iterations` | `int32_t` where counts are `size_t`; it is a diagnostic count, not a size, and `int32_t` is one of the boundary types. | keep | — | — |
| C14 | `po_portfolio_return`, `po_portfolio_volatility`, `po_utility`, `po_tracking_error`, `po_sharpe_ratio`, `po_diversification_ratio`, `po_risk_contributions` | Two carry `portfolio_`, five do not. A uniform prefix would rename five functions for no reader's benefit: `po_return` and `po_volatility` would be the ambiguous ones, and they are the two that carry it. | keep | — | — |
| C15 | Missing functions | Every use of the design is covered (data, spec, solve, the metrics of a weight vector, the frontier as a loop). Batch solving is `po_solve_batch`, out of 0.x. | add nothing | — | — |

## Python

| Id | Item | Finding | Decision | Breaks | Migration |
|---|---|---|---|---|---|
| P1 | `Objective`, `SolveStatus` | `IntEnum`s with the C values. | keep | — | — |
| P2 | `Spec` | Frozen dataclass whose keyword fields are the C setters; compares and hashes by value. | keep | — | — |
| P3 | `Data(mu, cov, *, risk_free, returns, benchmark)` | `mu` and `cov` are positional and required even when the objective needs neither: tracking error is written `Data(None, None, returns=..., benchmark=...)`. | change | no | `Data(None, None, returns=r, benchmark=b)` → also `Data(returns=r, benchmark=b)` (`mu=None`, `cov=None` by default) |
| P4 | `Data.portfolio_return`, `.volatility`, `.utility`, `.tracking_error`, `.sharpe_ratio`, `.diversification_ratio`, `.risk_contributions`, `.asset_count`, `.weight_count` | `volatility` is `po_portfolio_volatility`; the design names it `data.volatility` and it reads well next to `sharpe_ratio`. | keep | — | — |
| P5 | `Result` | Same fields as the C accessors; arrays read-only. | keep | — | — |
| P6 | `solve` | — | keep | — | — |
| P7 | `PortOptError` and the exception mapping | `PO_INVALID_ARGUMENT` is `ValueError` and `PO_UNSUPPORTED` is `NotImplementedError`, but `PO_OUT_OF_MEMORY` is a `PortOptError` where Python has `MemoryError`. | change | yes (an `except PortOptError` no longer catches an allocation failure) | `PortOptError(status=4)` → `MemoryError` with the native message |
| P8 | `portopt.__version__` | Missing; a user cannot tell which release is installed without `importlib.metadata`. | add | no | `portopt.__version__`, the version of the installed distribution (`pyportopt`), the one `pip` shows |
| P9 | `Spec.pointer`, `Data.pointer` | Public raw handles, used by `python/benchmarks/bench_overhead.py` to time the bare C calls. Harmless and useful for interop. | keep | — | — |

## C++ wrapper (planned)

The wrapper takes the Python names in `snake_case`, so one vocabulary serves both.

| Id | Item | Finding | Decision | Breaks | Migration |
|---|---|---|---|---|---|
| X1 | `po::Status`, `po::Error`, `po::Expected<T>` | `Expected`, not `Result`: `Result` is already the solve result. | keep | — | — |
| X2 | `po::Objective`, `po::SolveStatus` | Values equal to the C `#define`s. | keep | — | — |
| X3 | `po::Spec` (value, designated initializers), `po::Bounds`, `po::AssetBounds` | Mirrors the Python `Spec`; the native spec handle lives only inside `solve`. | keep | — | — |
| X4 | `po::Data` (move-only RAII handle) and its metrics | Same method names as Python `Data`. | keep | — | — |
| X5 | `po::Result`, `po::solve`, `po::version` | Same fields as Python `Result`. | keep | — | — |

## Status

Approved by the author on 2026-11-17 as proposed; every other item is kept. Done: C9, P3, P7, P8.

# API review before the 1.0 freeze

The four surfaces (C, Python, C++ and the .NET binding, unreleased), reviewed for what the 1.0 freeze makes
permanent. From 1.0 the contract only grows (`abi_contract.py check --base`, the libabigail dump, .NET package
validation against the last release), so a name or a shape that is wrong at 1.0 stays wrong. The question for each
item is not "does it work" but "can it grow without a break, and is its name true".

Already right, and kept: no public struct (configuration is setters on opaque handles, so a new option is a new
function); `int32_t` enums with `#define` values; one status type and a thread-local message; independent handles; an
immutable result; the enforcement above. Every future feature in sight fits additively: a tolerance, an iteration
limit, a warm start or a new objective is a new setter or a new `PO_OBJECTIVE_*`; a new scalar constraint is a new
`PO_CONSTRAINT_*` of `po_result_multiplier`; vector constraints (linear or group rows) would need a new accessor, which
is a new function.

## Growth rules (missing)

| Id | Item | Finding | Decision | Breaks | Migration |
|---|---|---|---|---|---|
| G1 | When the freeze starts | ARCHITECTURE.md says "freezes at 1.0"; nothing says which release is 1.0. Until then `check --base` has no base and the .NET package has no validation baseline. | change | — | the release after this review is `1.0.0`, with `PO_ABI_VERSION 1` |
| G2 | What a minor version may add | ExcelReader has `docs/guide/versioning.md`; here a user cannot tell what a 1.x may change. | add | no | a versioning guide: a minor version adds functions, setters, `#define` values and result fields (appended); it never removes, renames or changes a signature or a meaning; a new export does not change `PO_ABI_VERSION`; the bindings check the exact ABI version |
| G3 | Enum values the library *returns* (`po_status`, `PO_SOLVE_*`) | `PO_OBJECTIVE_*` and `PO_CONSTRAINT_*` grow safely: a caller only meets a value it chose. A new `PO_SOLVE_*` value is different: an old caller receives it. Python's `SolveStatus(value)` raises `ValueError` on an unknown value, so an old wheel loaded on a newer library (through `PORTOPT_LIB`) would fail a solve that succeeded. | add (rule) | no | rule: a new `PO_SOLVE_*` value is returned only to a request that uses something new (a setter or an objective an older caller cannot have set); any other new status is an ABI version change. A new `po_status` code is safe: every binding maps an unknown one to its generic error |
| G4 | Error messages | The header quotes messages ("`<setter> does not apply to <OBJECTIVE>`"), which reads as a promise. A message that is a contract cannot be improved. | change (doc) | no | one sentence in the header: messages are for people and may change in any version; the status is the contract |

## C interface

| Id | Item | Finding | Decision | Breaks | Migration |
|---|---|---|---|---|---|
| C16 | `po_spec_set_asset_bounds` (C6 kept it) | The same word names two counts: `po_data_asset_count` and `po_spec_set_risk_budget` count assets without cash, while `asset_bounds` takes one pair per *weight*, cash included. C6 judged the rename not worth it while the name could still change; at the freeze it becomes permanent. Python is not affected (it picks the function from the shape of `bounds`). | change | yes (C, C++, C#) | `po_spec_set_asset_bounds` → `po_spec_set_weight_bounds`; C++ `po::AssetBounds` → `po::WeightBounds`; C# `AssetBounds` → `WeightBounds` |
| C17 | `PO_SOLVE_MAX_ITER` | The header gives it a wider meaning than the name: "did not end at a point it can vouch for", including an unreachable target volatility and a gross exposure above 4.5e3. The name suggests a remedy (more iterations) the API does not offer: there is no iteration-limit setter. | change | yes (names only; the value 3 stays) | `PO_SOLVE_MAX_ITER` → `PO_SOLVE_NOT_CONVERGED` (Python `SolveStatus.NOT_CONVERGED`, C++ `not_converged`, C# `NotConverged`) |

## Python

| Id | Item | Finding | Decision | Breaks | Migration |
|---|---|---|---|---|---|
| P10 | `Data.asset_count`, `Data.weight_count` | Plain instance attributes: `data.weight_count = 3` is accepted and then lies. `Spec` and `Result` are frozen. | change | no | read-only properties |

## C++ wrapper

| Id | Item | Finding | Decision | Breaks | Migration |
|---|---|---|---|---|---|
| X6 | `po::Spec`, `po::Result` (aggregates) | Grow by a new member, which breaks no designated initializer as long as it is appended and the existing order stays. | keep (rule in G2) | — | — |
| X7 | `Data::create(assets, mu, cov)` plus `set_risk_free`, `set_returns` | Unlike Python and C#, inputs after the first two come through setters; that is also how a new input arrives without a break. | keep | — | — |

## .NET binding (unreleased: every change here is free)

| Id | Item | Finding | Decision | Breaks | Migration |
|---|---|---|---|---|---|
| N1 | `MarketData(double[]? mu = null, double[,]? covariance = null, double? riskFree = null, double[,]? returns = null, double[]? benchmark = null)` | Optional parameters are compiled into the caller. A sixth input changes the signature, so a consumer built against 1.0 fails with `MissingMethodException`; a second overload with optional parameters makes calls ambiguous. The constructor cannot grow. `Spec` already grows the right way (init properties). | change | no (unreleased) | `new MarketData(mu, cov, riskFree: r)` → `new MarketData(new MarketInputs { Mu = mu, Covariance = cov, RiskFree = r })`; `MarketInputs` is a record with init properties, and a new input is a new property |
| N2 | `Result.Weights`, `Result.BoundMultipliers` as `double[]` | Mutable arrays on an immutable result, but each `Result` owns a fresh copy, and an array passes straight to the `ReadOnlySpan<double>` metrics (`data.Volatility(result.Weights)`); `IReadOnlyList` would lose that. | keep | — | — |
| N3 | `Spec.Bounds` and `Spec.AssetBounds`, both settable | Both set is an `InvalidOperationException` at `Validate` or `Solve`. A closed hierarchy would make it unrepresentable at the cost of two more public types. | keep (renamed by C16) | — | — |
| N4 | `MarketData` where Python and C++ say `Data` | `Data` is too generic next to `System.Data` in .NET. | keep | — | — |

## Status

Approved by the author on 2026-11-21: the next release is 1.0.0 (G1), C16 and C17 renamed, N1 as proposed; every
other item as proposed. Done: G1 (`PO_ABI_VERSION` 1), G2 ([versioning.md](versioning.md)), G3, G4, C16, C17, P10, N1.
The review also found that `check --base` ran nowhere; the `abi-base` job of `ci.yml` runs it from the first v1.x
tag.
