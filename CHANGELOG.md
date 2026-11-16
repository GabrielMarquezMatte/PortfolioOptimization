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
- The leverage, a limit on the gross exposure `sum_i |w_i| <= L` of the risky weights (cash not counted):
  `po_spec_set_leverage`, the constraint `PO_CONSTRAINT_LEVERAGE` of `po_result_multiplier` (`lambda_lev <= 0`, NaN
  without a leverage, 0 when inactive) and a term `lambda_lev g` in the stationarity of the multiplier convention. It
  applies to every solved objective, `MAX_RETURN` included, and is solved as an implicit family of `2^n` rows that
  the QP generates as it needs them (`qp::GrossExposure`), never listed.
- `PO_OBJECTIVE_MAX_SHARPE` and `PO_OBJECTIVE_MAX_DIVERSIFICATION`, solved by homogenization (a QP in `y = kappa w`
  with `kappa = 1'y`). With cash both fix the cash weight at 0, and the budget term is carried by `beta_0`, with the
  budget multiplier NaN.
- `PO_SOLVE_UNBOUNDED`: the supremum of the ratio exists but is not attained (for example `MAX_SHARPE` with unbounded
  short selling and `1' Sigma^-1 mu~ <= 0`); the weights and multipliers are NaN. On an unbounded feasible set, so is
  a problem whose `1' Sigma^-1 v` is within the rounding of the data, where `kappa` is at most 16 times its noise model
  and its sign is not decided (docs/guide/validation.md). A bounded feasible set (a leverage, or no pair of distinct
  weights with an infinite upper and an infinite lower bound) is never `PO_SOLVE_UNBOUNDED`. `PO_SOLVE_INFEASIBLE`
  also covers the problem with no feasible portfolio of positive numerator, and cash bounds that exclude 0.
- `po_sharpe_ratio` and `po_diversification_ratio`.
- `PO_OBJECTIVE_RISK_PARITY`: minimizes `sum_i (c_i - b_i)^2` over the budget row, the bounds and, optionally, a target
  return (`po_spec_set_target_return`, needs `mu`) and a target volatility. New functions `po_spec_set_risk_budget` (the
  budget `b`, default `1 / n`), `po_spec_set_target_volatility` and `po_risk_contributions` (`c`, cash entry 0, NaN when
  the volatility is 0), with `po_result_residual` as the value of `f`. Without a target and with the Spinu point inside
  the box the answer is that point (Newton on the convex Spinu form); otherwise a Gauss-Newton sequence of QPs over the
  polyhedron of the linear rows, an augmented Lagrangian only for `vol(w) = s`, a Newton polish on the KKT system and a
  KKT and feasibility gate relative to the rounding noise of the problem. No cash and no leverage
  (`"risk_free does not apply to RISK_PARITY"`, `"leverage does not apply to RISK_PARITY"`); `target_volatility` and
  `risk_budget` apply to `RISK_PARITY` only. `PO_CONSTRAINT_VOLATILITY` is `df*/ds` (NaN without a target volatility;
  `-inf` when `s` equals the minimum volatility, where the minimum-variance portfolio is the only feasible point and the
  other multipliers are NaN). `PO_SOLVE_INFEASIBLE` only when certified. See
  [docs/guide/math.md](docs/guide/math.md#risk_parity).
- Python: `Spec.target_volatility`, `Spec.risk_budget` (stored as a tuple, so a `Spec` still compares and hashes by
  value), `Data.risk_contributions` and `Result.residual`.
- Every objective is now solved: `po_solve` no longer returns `PO_UNSUPPORTED` for a valid objective (the status and the
  Python `NotImplementedError` mapping remain).
- Python: `Spec.leverage`, `Data.sharpe_ratio`, `Data.diversification_ratio` and `Result.leverage_multiplier`; `Spec`
  with a leverage and array bounds still compares and hashes by value.
- The fuzz battery (`tests/unit/fuzz`): problem shapes over the condition number of the covariance (1e2 to 1e8), the
  daily and annual scales, equal and nearly equal means, a risk-free rate above every mean, tight boxes and cash, with
  references computed by another route (a bisection on the frontier, a golden-section search on the ratio, a brute
  force over the sign patterns of the leverage); see [docs/guide/validation.md](docs/guide/validation.md).

### Changed

- `PO_SOLVE_MAX_ITER` carries no point for every objective: weights, multipliers and `po_result_residual` are NaN (the
  header had left the case of `RISK_PARITY` open). For `RISK_PARITY` it is the status of a target volatility the
  iteration cannot reach, of a target return with no bounds (no minimizer) and of a solution whose gross exposure
  (above about 4.5e3) leaves the budget row unresolvable in double.
- The stationarity formula of the header has a volatility term for `RISK_PARITY` with a target volatility, as it had
  for `MAX_RETURN`.
- `MAX_RETURN` solves `Q(t)` shifted along a face of the rows (the budget, a cash bound or a row of the leverage), and
  an inactive cap is decided by a limit certificate (the KKT conditions of the linear program over the rows) rather
  than by comparing weights between two points of the path. The QP solver checks every row again at the final point.
  Ill-conditioned problems that returned `MAX_ITER` or a wrong point now solve; what is left is documented in
  [docs/guide/validation.md](docs/guide/validation.md#fuzz-battery).

- Python: `Spec` is no longer a plain class; it compares by value and cannot be assigned to after construction.

### Known limitations

- `RISK_PARITY` is a local method: `OPTIMAL` is a KKT point reached from the Spinu solution, not the global minimum
  (there are up to `2^(n-1)` normalized solutions of `c = b`; with negative lower bounds one in another orthant can be
  missed).
- With a target volatility and few or no active bounds at `cond(Sigma)` around 1e8, the weights are determined only to
  the rounding of the KKT residual amplified by the conditioning: up to about 1e-4 relative (measured up to 3.7e-5)
  from the exact KKT point, while the objective value is right to its band. Such a solve at `n = 50` can take about
  10 s; at `n <= 20` at most about 0.15 s.
- A bounded problem with the budget row only (Capped box, equal budget, `n = 20`, `cond(Sigma) = 1e8`, seed 31696) ends
  as `MAX_ITER`, without a point. Documented, not fixed.
- A target return without bounds has no minimizer in general and is `MAX_ITER`.
