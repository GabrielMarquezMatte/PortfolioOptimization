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
  (`po_data_set_risk_free`, cash last). At that step the other objectives returned `PO_UNSUPPORTED` (every
  objective is solved now, below); an infeasible problem is a result with status `PO_SOLVE_INFEASIBLE`, not an
  error.
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

- Performance tools: a Google Benchmark suite over every objective, variant and size (`BM_Solve`, `BM_Kernel`, 105
  cases of 10, 100 and 1000 assets), `portopt_profile` (one pass over the cases and a dump of the solutions in
  hex-float), `tests/benchmarks/quick.py` (a quick suite, A/B by medians and `compare-solutions`), `BM_CApi` and
  `python/benchmarks/bench_overhead.py` for the cost of the C and Python layers, all opt-in (`PO_BUILD_BENCHMARKS`);
  see [docs/performance/measuring.md](docs/performance/measuring.md).
- Performance reports: [docs/performance/benchmarks.md](docs/performance/benchmarks.md) (the before and after of every
  case), `baseline.md`, `profile-1.md`, `profile-2.md` and
  [docs/performance/python-vs-arch.md](docs/performance/python-vs-arch.md) (the C and Python layers against the solve).

### Changed

- Performance, with every result bit for bit identical (70 cases at 10 and 100 assets and all 35 at 1000 assets,
  `compare-solutions`): the initial basis of the Goldfarb-Idnani solver by rows (8 times at n = 1000) and kept along
  the `MAX_RETURN` path, the Cholesky factorization right-looking by column pairs, Householder reflections by rows,
  unit rows (bounds) read at their one entry, the symmetric block product of `gauss_newton`, Dekker's product in the
  compensated dot (where `std::fma` is a library call) and the allocations removed from the bound and gross exposure
  rows. Measured against the commit before them on the same machine and build (generic x86-64, GCC 16.2): the geometric
  mean of the solve time over the 35 cases is 0.71 at 10 assets, 0.29 at 100 and 0.23 at 1000 (the sum of the 35
  solves at 1000 assets went from 290 s to 79 s; `risk_parity/box` from 86.8 s to 21.5 s, `max_return/leverage` from
  32.9 s to 13.2 s). One case is slower, `risk_parity/budget` at 10 assets (5.5 to 5.9 µs). See
  [docs/performance/benchmarks.md](docs/performance/benchmarks.md#plan-6-final-before-and-after).
- `qp::solve` takes an optional `qp::BasisCache`; `Path` of `MAX_RETURN` holds one. Internal; no change of the C API.
- `-ffp-contract=off` on the C API objects as on the core and the tests (no accidental FMA in a `-march=native` build
  of the library).
- The test battery builds with MSVC (`/W4 /WX`) and passes there: `_CRT_SECURE_NO_WARNINGS` for the measurement tests,
  and five oracle self-tests that need an extended `long double` skip where `long double` is `double`.
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
- `RISK_PARITY` reports the multipliers of one more Newton solve of the KKT system at the returned weights (the
  budget, target-return and volatility multipliers and the bound multipliers of the held weights), kept when a
  safeguard accepts it, instead of the least-squares ones at the rounded weights. The weights, the status and the KKT
  gate do not change. The correction is made only at a point that passed the gate, and the rank test of the QR does
  not veto it (at `cond(Sigma) = 1e8` with `N >= 12` and no box the KKT matrix fails that test on about half the
  shapes while the correction is right). With a target volatility at least 1e-3 relative above the minimum
  volatility the measured relative error of the multipliers went from up to 1.5e-8, 1.6e-5 and 1.3 to up to 2.7e-13,
  5.2e-11 and 5.3e-10 at `cond(Sigma)` 1e4, 1e6 and 1e8 (`N` from 5 to 50; bound multipliers within 1.5e-10 at 1e8,
  `N` = 5 and 8); measured, not derived, and pinned by ceilings 10 times the largest measured on the cases of each
  test. Nearer the minimum volatility, `s = sigma_min (1 + delta)`, the error grows roughly like `1 / delta`: up to
  1.4e-10 at `delta = 1e-6` and 1.3e-6 at `1e-10`. At `cond(Sigma) >= 1e10` an `OPTIMAL` result is still not
  validated. See
  [docs/guide/validation.md](docs/guide/validation.md#the-kkt-residual-and-its-noise).

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
