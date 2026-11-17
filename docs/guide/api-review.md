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

Proposed on 2026-11-17; the author approves or overrides each decision before any change is made.
