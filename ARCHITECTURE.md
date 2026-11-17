# Architecture

A map of how the codebase fits together, not a manual. Start here, then follow the file names into the source:
the code carries its own meaning through naming and flow. What it deliberately does not carry is the evidence
behind a decision (measurements, rejected approaches); that lives in `docs/`, and
[STYLEGUIDE.md](STYLEGUIDE.md#comments) says why.

This document describes what exists today: nine objectives (`MIN_VARIANCE`, `MAX_RETURN`, `MAX_SHARPE`, `UTILITY`,
`IMPORTANCE`, `KELLY`, `TRACKING_ERROR`, `MAX_DIVERSIFICATION` and `RISK_PARITY`) with scalar or per-asset weight
bounds, an optional gross exposure limit (the leverage) and an optional risk-free asset (cash; `RISK_PARITY` takes
neither), eight of them reduced to a dense Goldfarb-Idnani quadratic program behind a C API (`MAX_RETURN` along a
parametric family of them, `MAX_SHARPE` and `MAX_DIVERSIFICATION` through a homogeneous one), and `RISK_PARITY`, a
non-convex least squares, solved by Newton on the Spinu form and, with a target or a bound that cuts it off, by a
sequence of such QPs and an augmented Lagrangian; portfolio return, volatility, utility, tracking error, the Sharpe and
diversification ratios and the risk contributions of given weights; and a Python package over the C API. Every
objective is solved: `PO_UNSUPPORTED` no longer occurs for a valid objective.

## Layers

```
  Python (python/src/portopt, ctypes + numpy)
                 |
        C API  include/portopt/portopt.h        <-- the only stable contract
                 |
        src/capi   handles, error boundary
          |                 \
        src/solve            src/metrics
        requirements,        portfolio return
        dispatch             and volatility
          |       \              |
        src/riskparity          |
        Spinu Newton, QP        |
        sequence, Newton-KKT    |
          |                     |
        src/problems            |
        spec and data           |
        to a QP and back        |
          |       \              |
          |     src/qp           |
          |     Goldfarb-Idnani  |
          |        |             |
        src/core   spec, data, result, failure
          |        |
        src/linalg   Cholesky, triangular solves, Householder QR
```

Dependencies point one way, in the order `core` and `qp` <- `problems` <- `riskparity` <- `solve` <- `capi`:
`core` (objective, spec, data, result, failure) <- `problems` (mean_variance,
tracking_error, max_return, parametric_path, homogeneous, constraints: translate a spec and data
to a QP and read the solution back) <- `riskparity` (spinu: the positive-orthant solution of
`RISK_PARITY` by Newton; objective, polyhedral_descent: the derivatives of f and the Gauss-Newton QP sequence over
the polyhedron; augmented_lagrangian: the target volatility; kkt, newton_polish: the KKT residual against its
noise and Newton on the KKT system; risk_parity: the policy over them, 7 modules, one below the limit of 8 per folder)
<- `solve` (requirements, solve: the checks of the request and the dispatch to an objective, including the one to
`riskparity`) <- `capi`. `linalg` (Cholesky, triangular solves and Householder QR) depends on nothing; `qp` (Goldfarb-Idnani in `goldfarb_idnani.*` and `active_set.hpp`) depends
only on `linalg`. `core/data.cpp` uses the Cholesky factor to reject a covariance that is not positive definite, and
`problems` uses `qp`. `metrics` depends on `core` and is used by `riskparity` and `capi`. Nothing in `src/` includes from `src/capi/`
except `capi` itself. Everything below the C API is internal and may change in any release; only the header is a
contract, and the Python package is a consumer of it like any other wrapper would be.

Build targets (`CMakeLists.txt`):

| Target | What it is |
|---|---|
| `portopt_core` | Static library: `linalg`, `qp`, `core`, `problems`, `riskparity`, `solve`, `metrics` |
| `portopt_capi` | OBJECT library: every `src/capi/*.cpp`, `error.cpp` included, compiled once |
| `portopt` | The shared library, built from the `portopt_capi` objects; links `portopt_core` |
| `portopt_tests` | GoogleTest executable; links `portopt_capi` and `portopt_core`, so a test sees the same thread-local error as the entry points it calls (static `libstdc++` on MinGW) |
| `portopt_c_header_check` | Compiles the public header as C11 |
| ctest `style` | `tests/style/check_style.py` |
| `portopt_benchmarks` | Google Benchmark suite (`BM_Solve/*`, `BM_Kernel/*`); opt-in (`PO_BUILD_BENCHMARKS`, preset `gcc-bench`); links `portopt_core` directly and counts allocations with a replaced `operator new` that no other target has |
| `portopt_profile` | One pass over every benchmark case: the per-case table and the dump of solutions in hex-float for `quick.py compare-solutions`; opt-in like `portopt_benchmarks` |

## Where each responsibility lives

| Path | Responsibility |
|---|---|
| `include/portopt/portopt.h` | The C API: `po_` functions, `PO_` constants, opaque handles, no public struct |
| `src/capi/handles.hpp` | The handle types (`po_spec`, `po_data`, `po_result`), each wrapping one core value, and the null-handle and failure helpers |
| `src/capi/error.hpp`, `error.cpp` | Thread-local `po_last_error()`, `set_error`, and `guard`, the one place exceptions are caught |
| `src/capi/spec_api.cpp`, `data_api.cpp`, `solve_api.cpp`, `result_api.cpp`, `metrics_api.cpp`, `version_api.cpp` | The `extern "C"` functions, one file per handle family |
| `src/core/objective.hpp` | `po::Objective`, matching the `PO_OBJECTIVE_*` values |
| `src/core/spec.hpp` | `po::Spec`: objective, scalar bounds, optional per-asset bounds, optional target return, risk aversion, importance, maximum volatility, leverage, target volatility and risk budget; a copyable value |
| `src/core/data.hpp`, `data.cpp` | `po::Data`: mu, covariance, its Cholesky factor, the optional risk-free rate and the optional return series with its benchmark; `Data::create` validates and symmetrizes the covariance and factors it once |
| `src/core/result.hpp` | `po::Result` and `SolveStatus`: weights, status, multipliers, iterations; the multiplier convention is stated there |
| `src/core/failure.hpp` | `Failure` (kind and message), the error type of `std::expected` below the C API |
| `src/linalg/cholesky.*` | Cholesky factor of a symmetric matrix (`nullopt` if not positive definite to working precision) and the two triangular solves; no explicit inverse |
| `src/linalg/qr.*` | Householder QR of a tall matrix without forming `a'a`: the rank decision on `diag(R)`, least squares, and `R'` as the lower factor of `a'a` |
| `src/qp/problem.hpp` | `qp::Problem` (G, an optional precomputed factor of G, g, equality and inequality rows, and the optional `GrossExposure` family) and `qp::Solution` (with the sign vectors and multipliers of the active generated rows), with the stationarity convention |
| `src/qp/goldfarb_idnani.*` | `qp::solve`: the dual active-set loop, the tolerances, the generation of the rows of the gross exposure family and the final refinement |
| `src/qp/active_set.hpp` | The active-set factorization (J, R), Givens rotations, the add and drop updates and the refinement step |
| `src/solve/solve.*` | `po::solve`: `check_request`, then dispatch to the objective |
| `src/riskparity/spinu.*` | `spinu`: Newton on `F(y) = 1/2 y' Sigma y - b' log y`, `y > 0`, `w = y / 1'y`; Cholesky and two triangular solves per step, Armijo far from the minimizer and the full step near it, the stop on the relative size of the step |
| `src/riskparity/objective.*` | `linearize`: `f = sum_i (c_i - b_i)^2`, its gradient `2 J'r`, the Jacobian `J` of the risk contributions and the exact Hessian; `value_change`, `f(w + s) - f(w)` computed from `s` |
| `src/riskparity/polyhedral_descent.*` | `descend`: proximal Gauss-Newton steps over the polyhedron (budget, target return, bounds), each a QP of `qp::solve` posed in the step, with a Levenberg-Marquardt proximal weight in the units of `2 J'J` (so the iterates do not depend on the scale of Sigma), Armijo on the Lagrangian of the step's multipliers, an optional augmented-Lagrangian volatility term and a stop reason; the QP multipliers at the fixed point are the user's |
| `src/riskparity/kkt.*` | The KKT system of `RISK_PARITY` at a point (`evaluate_kkt`: the gradients of the equalities, the Hessian of the Lagrangian, the stationarity and the gaps), each residual with its noise (how far one rounding of every weight and of every term moves it), `residual_of` in units of that noise and `is_kkt_point`, the one decision of what counts as a KKT point |
| `src/riskparity/newton_polish.*` | `newton_polish`: Newton on the KKT system with the exact Hessian of the Lagrangian over the free weights (Householder QR of the system), the multipliers estimated again by least squares at every step, a bound released when its multiplier has the wrong sign, and a step taken only when it decreases the KKT residual in noise units |
| `src/riskparity/augmented_lagrangian.*` | `augmented_lagrangian`: the multiplier method of Andreani, Birgin, Martinez and Schuverdt (2007, Algorithm 3.1) for the row `vol(w) = s` alone, each subproblem a `descend` over the polyhedron; penalty and multiplier in the units of the relative gap `vol(w) / s - 1` |
| `src/riskparity/risk_parity.*` | `RISK_PARITY`: `risk_budget_of` (the setter or `1 / n`) and `solve_risk_parity`, the policy: the Spinu point (the answer when it is inside the box and no target is set), the minimum-variance QP for feasibility and `sigma_min` (`INFEASIBLE`, or the case `s = sigma_min` with `lambda_vol = -inf`), then the descent or the augmented Lagrangian from the point of the polyhedron nearest to the Spinu point, the polish, and `is_kkt_point` deciding `OPTIMAL` against `MAX_ITER` |
| `src/solve/requirements.*` | `check_request`: the objective is implemented, every setter applies to it, the data it requires is present and per-asset bounds match the weight count (the objective x setter matrix of [docs/guide/math.md](docs/guide/math.md#which-setters-and-data-each-objective-takes)) |
| `src/problems/mean_variance.*` | One builder for `MIN_VARIANCE`, `UTILITY`, `IMPORTANCE` and `KELLY`: `G = a Sigma` with factor `sqrt(a) L` and `g = -b mu~`, `(a, b)` chosen by the objective; `solve_mean_variance` |
| `src/problems/tracking_error.*` | `TRACKING_ERROR`: `G = 2 R~'R~` with the factor `sqrt(2) R_qr'` from one Householder QR of `R~`, which decides the column rank; `g = -2 R~'(b - rf 1)` |
| `src/problems/parametric_path.*` | The path `Q(t)` of `MAX_RETURN`: the shift of the means along a face (the budget, a cash bound or a row of the leverage) and its sign certificate, solving at `t`, the segment (active set) and constant-weight decisions with their rounding bounds, the multiplier slope by least squares, the limit certificate of an inactive cap, the vertex exit, the crossing on a segment |
| `src/problems/max_return.*` | `MAX_RETURN`: the policy over the path (climb, bisection, the case `s = sigma_min`) and the volatility multiplier `-s / t*` |
| `src/problems/homogeneous.*` | `MAX_SHARPE` and `MAX_DIVERSIFICATION`: the homogeneous QP `min y' Sigma y` with `v'y = 1` and the rows `1'y >= 0` and, with a leverage, the family with offset `L 1`; the status from the active set (`INFEASIBLE`, `UNBOUNDED`, `OPTIMAL`) and the multipliers mapped back to `w = y / kappa` |
| `src/problems/constraints.*` | The constraint rows (budget, target return, bounds), their bookkeeping (`ConstraintRows`) and `to_result`, the remapping to the user's variables, with an optional scale for the multipliers; `lower_of` and `upper_of`, the bound of a weight, and `kMinimumVolatilityTolerance`, shared by `MAX_RETURN` and `RISK_PARITY` |
| `src/metrics/portfolio.*` | `portfolio_return`, `portfolio_volatility`, `utility`, `sharpe_ratio`, `diversification_ratio`, `tracking_error`, `risk_contributions` and `risk_budget_residual` for a weight vector (cash last) |
| `python/src/portopt/_native.py` | Library discovery (`PORTOPT_LIB`), ABI check, prototypes, handle ownership, status to exception |
| `python/src/portopt/_spec.py`, `_data.py`, `_result.py` | `Objective` and `Spec`, a frozen dataclass that builds its `po_spec` handle in `__post_init__` (bounds, `leverage`, `target_volatility` and `risk_budget` included; the bounds and the budget canonicalized to tuples so that equality and hash work); `Data` and its metrics (`risk_contributions` too); `SolveStatus`, `Result` and `solve` |
| `python/hatch_build.py` | Tags the wheel with the platform of the bundled library |
| `tests/unit/{capi,core,linalg,metrics,problems,qp,riskparity,solve}` | GoogleTest suites mirroring `src/`; `capi/c_header_check.c` compiles the header as C11 |
| `tests/unit/fuzz` | The fuzz battery: a generator of problem shapes (`random_problems.hpp`: condition number, scale, mean pattern, box, cash, leverage) and references computed by another route (`frontier_reference.hpp`, `sign_pattern_reference.hpp`, `ratio_reference.hpp`), used by the `MaxReturnFuzz`, `LeverageFuzz` and `RatioFuzz` suites; `tests/unit/fuzz/risk_parity/` holds the risk parity references (the scan of the feasible curve, the KKT check with its rounding floor, and the refinement of a result in `long double` that measures its forward error), the `RiskParityDescentFuzz` sweep of the polyhedral descent, the `RiskParitySolveSweep` of `po::solve` and the `RiskParityFuzz` battery |
| `tests/unit/kkt.hpp`, `user_kkt.hpp`, `random_covariance.hpp`, `capi_scope.hpp` | KKT residuals of a QP (the generated rows of the family included), the KKT conditions in the user's variables for every objective, random covariances of a chosen condition number, RAII for C handles |
| `tests/benchmarks` | The benchmark and profile tools: `instances.hpp` and `cases.hpp` (deterministic instances, cached, and the objective x variant x size table), `bm_solve.cpp`, `bm_kernel.cpp`, `bench_main.cpp`, `profile.cpp`, `alloc_counter.*` and `quick.py` (`run`, `compare`, `compare-solutions`, `selftest`); how to use them is in [docs/performance/measuring.md](docs/performance/measuring.md) |
| `tests/style/check_style.py` | Style rules that can be checked without a compiler; runs as the `style` ctest |
| `python/tests` | Python tests; they skip cleanly (pytest exit 5) when the native library is missing |
| `docs/guide/math.md` | Notation, the QP translation of each objective, the gross exposure family, the parametric path of MAX_RETURN, the homogenization of MAX_SHARPE and MAX_DIVERSIFICATION, RISK_PARITY (Spinu, the QP sequence, the augmented Lagrangian, the multipliers) and the multiplier convention |
| `docs/performance/measuring.md` | How to measure: the `gcc-bench` and `msvc-profile` presets, the quick suite, A/B by medians against noise, the comparison of solutions, the uProf recipe and the guards an optimization must keep green |
| `docs/performance/benchmarks.md`, `baseline.md`, `profile-1.md` | Results of the suite and the machine they come from; the plan 6 baseline (every case by objective, variant and size, and the kernels); where the time goes before plan 6 (uProf shares, calls per solve, assembly) and the hot-spot hypotheses it confirms or refutes |
| `docs/guide/validation.md` | Validation policy, tolerances, the numerical decisions of the solver and of each objective, and the fuzz battery |

## Stability rules of the C API

Each rule answers a concrete break that happened in the author's earlier libraries (ExcelReader, FixedWidthParser,
VolEngine).

| What went wrong before | Rule here |
|---|---|
| Loose positional parameters and `_ex` twins of functions | No public struct. Configuration goes through setters on opaque handles. New configuration is a new function, never a new parameter or an `_ex` variant. |
| An enum value removed | Enums are `typedef int32_t` plus `#define`; values are only added, never removed or renumbered. |
| Errors reported in many ways | Every fallible function returns `po_status` and puts the message in `po_last_error()`, which is thread-local and `""` after a success. |
| An exception crossing the boundary | None does: `guard` (`src/capi/error.cpp`) turns `std::bad_alloc` into `PO_OUT_OF_MEMORY` and anything else into `PO_INTERNAL_ERROR`. |
| Mutable global state | None. Handles are independent: freeing one never invalidates another. A `po_result` is immutable. |
| A convention assumed silently | One multiplier convention, written in the header and in [docs/guide/math.md](docs/guide/math.md#multiplier-convention). |
| ABI drift | Semantic versioning, with `PO_ABI_VERSION` in the header and `po_abi_version()` at run time; the contract freezes at 1.0 and 0.x releases may break it. |

Types on the boundary are `int32_t`, `size_t`, `double` and pointers. A null handle is `PO_INVALID_HANDLE`, a null
output pointer or a null required input array `PO_INVALID_ARGUMENT` (the one exception: `mu` and `cov` of
`po_data_create` may be `NULL`), and every `*_free` accepts `NULL`.

## Anatomy of a `po_solve`

The Python `solve` builds the handles below and then makes one call; any C client does the same.

1. **Handles.** `po_spec_create` and the spec setters validate bounds, the target return, the risk aversion, the
   importance, the maximum volatility, the leverage, the target volatility and the risk budget as they are set (a rejected call leaves the spec unchanged). `po_data_create` copies mu and the covariance, rejects a non-finite
   entry, a covariance asymmetric beyond 1e-10 of the larger diagonal entry or one that is not positive definite
   (`Data::create`, by Cholesky, with a pivot floor relative to the diagonal; see
   [docs/guide/validation.md](docs/guide/validation.md#positive-definiteness)), and stores the covariance
   symmetrized together with its Cholesky factor L: this is the only place positive definiteness is decided.
   `po_data_set_risk_free` adds cash as the last weight and `po_data_set_returns` stores the return series and the
   benchmark of `TRACKING_ERROR`.
2. **Validation.** `po_solve` (`src/capi/solve_api.cpp`) sets `*out` to NULL, checks the handles inside `guard` and
   calls `po::solve` (`src/solve/solve.cpp`). Its `check_request` (`src/solve/requirements.cpp`) rejects, in
   this order, a setter outside the objective's row of the matrix, data the objective requires (`RISK_PARITY` rejects
   cash too) and per-asset bounds whose size is not the weight count, cash included, or a risk budget whose size is not
   the asset count (`PO_INVALID_ARGUMENT`; a core `Unsupported` failure remains for an objective without a row, which no
   valid objective is). Each is an error status, not a result. `solve` then dispatches on the objective.
3. **Translation to a QP.** `mean_variance_problem` (`src/problems/mean_variance.cpp`) takes the risky weights as
   the QP variables, with G = a Sigma (passed with its factor sqrt(a) L) and g = -b mu~; `(a, b)` is (2, 0) for
   `MIN_VARIANCE`, (delta, 1) for `UTILITY`, (2 (1 - lambda), lambda) for `IMPORTANCE` and (1 / (1 + rf), 1) for
   `KELLY`. `tracking_error_problem` builds G = 2 R~'R~ (passed with the factor sqrt(2) R_qr' of a Householder QR
   of R~) and g = -2 R~'(b - rf 1) from the return series, and `solve_max_return` does not build one QP but `Q(t)`,
   with G = Sigma and g = -t (mu~ - c r), shifted along a face r of the rows (the budget, a cash bound or a row of the
   leverage; c the midrange of the r_i mu~_i), solved at several `t` (`src/problems/parametric_path.cpp`).
   `solve_risk_parity` (`src/riskparity/risk_parity.cpp`) builds no single QP: it solves Spinu's convex problem by Newton,
   and otherwise a sequence of QPs over the polyhedron of the linear rows (see
   [docs/guide/math.md](docs/guide/math.md#risk_parity)).
   `solve_max_sharpe` and `solve_max_diversification` (`src/problems/homogeneous.cpp`) build the homogeneous QP in
   y = kappa w: G = 2 Sigma, v'y = 1 as the equality, one row per finite bound written over y and the row 1'y >= 0, and
   with a leverage the family with offset L 1. The leverage of any other objective is the same `GrossExposure` on the
   weights (a = 0, b = -L). The builders other than the homogeneous one use `add_constraint_rows`
   (`src/problems/constraints.cpp`), which appends the budget row (no cash), the target-return row and one row per
   finite bound, recording in `ConstraintRows` which row is which. With cash, the cash weight 1 - 1'w is eliminated,
   so its bounds become rows over the risky weights. The rows are in [docs/guide/math.md](docs/guide/math.md#the-qp-rows).
4. **Goldfarb-Idnani.** `qp::solve` (`src/qp/goldfarb_idnani.cpp`) uses the factor of G it is given (it factors
   G by Cholesky only when none is given), starts at the unconstrained minimum and adds constraints one at a time:
   the equalities first, then the most violated inequality, each with a primal step bounded by the dual step that drops an active inequality whose multiplier
   reaches zero. The active set keeps the factorization J'N = [R; 0] updated by Givens rotations and never forms an
   inverse. Five numerical decisions are fixed by tests (details in
   [docs/guide/validation.md](docs/guide/validation.md#numerical-decisions-of-the-qp-solver)): feasibility checked
   again at the final x, no generated row equal to an active row, a scale-free dependence test, tolerances relative to
   the row and the point, and one step of iterative refinement over the active set at the optimum. With a
   `GrossExposure` family the most violated inequality may be a generated row: the row `(a - s)'x >= b` at
   `s = sign(x)`, which competes with the explicit rows by the same rule and, once active, is numbered after them (at
   most n are active at once; [the family](docs/guide/validation.md#the-gross-exposure-family)). A violated constraint
   with no feasible step is `Infeasible`; more than 10 (n + m) + 10 iterations is `MaxIterations`, the family
   counting as n rows.
5. **Remapping to the user's variables.** `to_result` (`src/problems/constraints.cpp`) turns the QP solution into
   a `Result`: with cash, the last weight is 1 - 1'w; the budget and target-return multipliers are read from the
   equality rows (NaN when the constraint does not exist); each bound multiplier is the lower-row multiplier minus the
   upper-row multiplier, 0 for a weight with no finite bound. A status other than `Optimal` (for a QP, `Infeasible`
   or `MaxIterations`, whose dual iterate is not primal feasible) gives NaN weights, NaN multipliers and a NaN
   residual, never a partial point; the same rule holds for the `MaxIterations` of `RISK_PARITY`. The QP objective equals the user's minimization form up to a constant, so
   the multipliers need no rescaling, with one exception: `MAX_RETURN` divides the multipliers of `Q(t*)` by `t*`
   (`multiplier_scale`) and sets the volatility multiplier to `-s / t*`; for an inactive cap the multipliers are
   the slope `b` of `N b = -mu~`, with scale 1. The leverage multiplier is `-sum_j zeta_j` over the multipliers of the
   active generated rows (0 when none is active, NaN without a leverage). `MAX_SHARPE` and `MAX_DIVERSIFICATION` map
   the multipliers of the homogeneous QP back with the factor `kappa / (2 sigma_y^3)`
   ([docs/guide/math.md](docs/guide/math.md#max_sharpe-and-max_diversification)); a QP that is `Infeasible` is
   `INFEASIBLE`, and an active row `1'y >= 0` is `UNBOUNDED`, with NaN weights and multipliers like `INFEASIBLE`.
6. **Publishing.** The `Result` is moved into a freshly allocated `po_result` and its pointer written to `*out` only
   once everything has succeeded, so a failure leaves `*out` NULL. The Python layer reads the result into a frozen
   dataclass, with read-only weight and multiplier arrays, and frees the handle.

Every `extern "C"` function runs inside `guard`, the single place that turns exceptions into statuses.
