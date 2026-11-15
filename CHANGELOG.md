# Changelog

All notable changes to this project are documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
Until 1.0 the C ABI may break between 0.x releases (`PO_ABI_VERSION` is 0 and the Python package refuses a
library whose ABI version differs).

## [Unreleased]

### Added

- The C API in `include/portopt/portopt.h`: opaque `po_spec`, `po_data` and `po_result` handles, the complete
  `PO_OBJECTIVE_*` enum, `po_status` with a thread-local `po_last_error()`, and one error boundary through which no
  exception crosses.
- `PO_OBJECTIVE_MIN_VARIANCE`: minimum variance with a scalar box (`po_spec_set_bounds`), per-asset bounds
  (`po_spec_set_asset_bounds`), a target return (`po_spec_set_target_return`) and a risk-free asset
  (`po_data_set_risk_free`, cash last). The other objectives return `PO_UNSUPPORTED`; an infeasible problem is a
  result with status `PO_SOLVE_INFEASIBLE`, not an error.
- Multipliers of the budget and target-return constraints and of the bounds, under one documented sign convention.
- A dense Goldfarb-Idnani dual active-set quadratic programming solver, written from the paper.
- Portfolio metrics: `po_portfolio_return` and `po_portfolio_volatility`.
- The Python package `portopt` (`Spec`, `Data`, `solve`, `Result`), a ctypes layer over the C API; the library is
  found through `PORTOPT_LIB` or `python/src/portopt/_lib`.
- Tests: KKT conditions of random and ill-conditioned problems, closed forms (unconstrained minimum, frontier,
  capital market line) and the C boundary; see [docs/guide/validation.md](docs/guide/validation.md).
- CI on Windows (MSVC), Linux (GCC 15) and macOS (Clang), with `clang-tidy` over `src/` on Linux.
- `PO_OBJECTIVE_UTILITY`, `PO_OBJECTIVE_IMPORTANCE` and `PO_OBJECTIVE_KELLY`, built with `PO_OBJECTIVE_MIN_VARIANCE` by
  one mean-variance builder, and `PO_OBJECTIVE_TRACKING_ERROR`, from a return series and a benchmark. New functions:
  `po_spec_set_risk_aversion`, `po_spec_set_importance`, `po_spec_set_max_volatility`, `po_data_set_returns`,
  `po_utility` and `po_tracking_error`.
- `PO_OBJECTIVE_MAX_RETURN`: maximum return under a volatility cap, solved along a parametric family of QPs, with the
  volatility multiplier `PO_CONSTRAINT_VOLATILITY` (`-s / t*`; `-inf` at a vertical frontier, documented in the
  header).
- A check of the objective x setter matrix in `po_solve`: a setter outside the objective's row or a missing input
  is `PO_INVALID_ARGUMENT` with `"<setter> does not apply to <OBJECTIVE>"` or `"<OBJECTIVE> requires <item>"`.
- Python: `Spec` is a frozen dataclass with the new keyword fields (use `dataclasses.replace` to vary one),
  `Data(..., returns=, benchmark=)`, `Data.utility`, `Data.tracking_error` and `Result.volatility_multiplier`.
- Tests: consistency between objectives (`IMPORTANCE` and `KELLY` against `UTILITY`), `MAX_RETURN` against the
  frontier and the capital market line, and its volatility multiplier against finite differences; see
  [docs/guide/validation.md](docs/guide/validation.md).

### Changed

- Python: `Spec` is no longer a plain class; it compares by value and cannot be assigned to after construction.
