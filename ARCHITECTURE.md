# Architecture

A map of how the codebase fits together, not a manual. Start here, then follow the file names into the source:
the code carries its own meaning through naming and flow. What it deliberately does not carry is the evidence
behind a decision (measurements, rejected approaches); that lives in `docs/`, and
[STYLEGUIDE.md](STYLEGUIDE.md#comments) says why.

This document describes what exists today: six objectives (`MIN_VARIANCE`, `MAX_RETURN`, `UTILITY`, `IMPORTANCE`,
`KELLY` and `TRACKING_ERROR`) with scalar or per-asset weight bounds and an optional risk-free asset (cash), each
reduced to a dense Goldfarb-Idnani quadratic program behind a C API (`MAX_RETURN` along a parametric family of
them); portfolio return, volatility, utility and tracking error of given weights; and a Python package over the C
API. The other objectives of the enum return `PO_UNSUPPORTED`.

## Layers

```
  Python (python/src/portopt, ctypes + numpy)
                 |
        C API  include/portopt/portopt.h        <-- the only stable contract
                 |
        src/capi   handles, error boundary
          |                 \
        src/problems         src/metrics
        spec and data        portfolio return
        to a QP and back     and volatility
          |       \              |
          |     src/qp           |
          |     Goldfarb-Idnani  |
          |        |             |
        src/core   spec, data, result, failure
          |        |
        src/linalg   Cholesky, triangular solves, Householder QR
```

Dependencies point one way: `core` (objective, spec, data, result, failure) <- `problems` (requirements,
mean_variance, tracking_error, max_return, parametric_path, constraints, solve: translate a spec and data to a QP
and read the solution back) <- `capi`. `linalg` (Cholesky, triangular solves and Householder QR) depends on nothing; `qp` (Goldfarb-Idnani in `goldfarb_idnani.*` and `active_set.hpp`) depends
only on `linalg`. `core/data.cpp` uses the Cholesky factor to reject a covariance that is not positive definite, and
`problems` uses `qp`. `metrics` depends on `core` and is used by `capi`. Nothing in `src/` includes from `src/capi/`
except `capi` itself. Everything below the C API is internal and may change in any release; only the header is a
contract, and the Python package is a consumer of it like any other wrapper would be.

Build targets (`CMakeLists.txt`):

| Target | What it is |
|---|---|
| `portopt_core` | Static library: `linalg`, `qp`, `core`, `problems`, `metrics` |
| `portopt_capi` | OBJECT library: every `src/capi/*.cpp`, `error.cpp` included, compiled once |
| `portopt` | The shared library, built from the `portopt_capi` objects; links `portopt_core` |
| `portopt_tests` | GoogleTest executable; links `portopt_capi` and `portopt_core`, so a test sees the same thread-local error as the entry points it calls (static `libstdc++` on MinGW) |
| `portopt_c_header_check` | Compiles the public header as C11 |
| ctest `style` | `tests/style/check_style.py` |

## Where each responsibility lives

| Path | Responsibility |
|---|---|
| `include/portopt/portopt.h` | The C API: `po_` functions, `PO_` constants, opaque handles, no public struct |
| `src/capi/handles.hpp` | The handle types (`po_spec`, `po_data`, `po_result`), each wrapping one core value, and the null-handle and failure helpers |
| `src/capi/error.hpp`, `error.cpp` | Thread-local `po_last_error()`, `set_error`, and `guard`, the one place exceptions are caught |
| `src/capi/spec_api.cpp`, `data_api.cpp`, `solve_api.cpp`, `result_api.cpp`, `metrics_api.cpp`, `version_api.cpp` | The `extern "C"` functions, one file per handle family |
| `src/core/objective.hpp` | `po::Objective`, matching the `PO_OBJECTIVE_*` values |
| `src/core/spec.hpp` | `po::Spec`: objective, scalar bounds, optional per-asset bounds, optional target return, risk aversion, importance and maximum volatility; a copyable value |
| `src/core/data.hpp`, `data.cpp` | `po::Data`: mu, covariance, its Cholesky factor, the optional risk-free rate and the optional return series with its benchmark; `Data::create` validates and symmetrizes the covariance and factors it once |
| `src/core/result.hpp` | `po::Result` and `SolveStatus`: weights, status, multipliers, iterations; the multiplier convention is stated there |
| `src/core/failure.hpp` | `Failure` (kind and message), the error type of `std::expected` below the C API |
| `src/linalg/cholesky.*` | Cholesky factor of a symmetric matrix (`nullopt` if not positive definite to working precision) and the two triangular solves; no explicit inverse |
| `src/linalg/qr.*` | Householder QR of a tall matrix without forming `a'a`: the rank decision on `diag(R)`, least squares, and `R'` as the lower factor of `a'a` |
| `src/qp/problem.hpp` | `qp::Problem` (G, an optional precomputed factor of G, g, equality and inequality rows) and `qp::Solution`, with the stationarity convention |
| `src/qp/goldfarb_idnani.*` | `qp::solve`: the dual active-set loop, the tolerances and the final refinement |
| `src/qp/active_set.hpp` | The active-set factorization (J, R), Givens rotations, the add and drop updates and the refinement step |
| `src/problems/solve.*` | `po::solve`: `check_request`, then dispatch to the objective |
| `src/problems/requirements.*` | `check_request`: the objective is implemented, every setter applies to it, the data it requires is present and per-asset bounds match the weight count (the objective x setter matrix of [docs/guide/math.md](docs/guide/math.md#which-setters-and-data-each-objective-takes)) |
| `src/problems/mean_variance.*` | One builder for `MIN_VARIANCE`, `UTILITY`, `IMPORTANCE` and `KELLY`: `G = a Sigma` with factor `sqrt(a) L` and `g = -b mu~`, `(a, b)` chosen by the objective; `solve_mean_variance` |
| `src/problems/tracking_error.*` | `TRACKING_ERROR`: `G = 2 R~'R~` with the factor `sqrt(2) R_qr'` from one Householder QR of `R~`, which decides the column rank; `g = -2 R~'(b - rf 1)` |
| `src/problems/parametric_path.*` | The path `Q(t)` of `MAX_RETURN`: the shift of the means along the budget, solving at `t`, the segment (active set) and constant-weight decisions with their rounding bounds, the multiplier slope by least squares, the limit certificate of an inactive cap, the vertex exit, the crossing on a segment |
| `src/problems/max_return.*` | `MAX_RETURN`: the policy over the path (climb, bisection, the case `s = sigma_min`) and the volatility multiplier `-s / t*` |
| `src/problems/constraints.*` | The constraint rows (budget, target return, bounds), their bookkeeping (`ConstraintRows`) and `to_result`, the remapping to the user's variables, with an optional scale for the multipliers |
| `src/metrics/portfolio.*` | `portfolio_return`, `portfolio_volatility`, `utility` and `tracking_error` for a weight vector (cash last) |
| `python/src/portopt/_native.py` | Library discovery (`PORTOPT_LIB`), ABI check, prototypes, handle ownership, status to exception |
| `python/src/portopt/_spec.py`, `_data.py`, `_result.py` | `Objective` and `Spec`, a frozen dataclass that builds its `po_spec` handle in `__post_init__`; `Data` and its metrics; `SolveStatus`, `Result` and `solve` |
| `python/hatch_build.py` | Tags the wheel with the platform of the bundled library |
| `tests/unit/{capi,core,linalg,metrics,problems,qp}` | GoogleTest suites mirroring `src/`; `capi/c_header_check.c` compiles the header as C11 |
| `tests/unit/kkt.hpp`, `random_covariance.hpp`, `capi_scope.hpp` | KKT residuals of a QP, random covariances of a chosen condition number, RAII for C handles |
| `tests/style/check_style.py` | Style rules that can be checked without a compiler; runs as the `style` ctest |
| `python/tests` | Python tests; they skip cleanly (pytest exit 5) when the native library is missing |
| `docs/guide/math.md` | Notation, the QP translation of each objective, the parametric path of MAX_RETURN and the multiplier convention |
| `docs/guide/validation.md` | Validation policy, tolerances and numerical decisions of the solver |

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
   importance and the maximum volatility as they are set (a rejected call leaves the spec unchanged). `po_data_create` copies mu and the covariance, rejects a non-finite
   entry, a covariance asymmetric beyond 1e-10 of the larger diagonal entry or one that is not positive definite
   (`Data::create`, by Cholesky, with a pivot floor relative to the diagonal; see
   [docs/guide/validation.md](docs/guide/validation.md#positive-definiteness)), and stores the covariance
   symmetrized together with its Cholesky factor L: this is the only place positive definiteness is decided.
   `po_data_set_risk_free` adds cash as the last weight and `po_data_set_returns` stores the return series and the
   benchmark of `TRACKING_ERROR`.
2. **Validation.** `po_solve` (`src/capi/solve_api.cpp`) sets `*out` to NULL, checks the handles inside `guard` and
   calls `po::solve` (`src/problems/solve.cpp`). Its `check_request` (`src/problems/requirements.cpp`) rejects, in
   this order, an objective this version does not solve (`PO_UNSUPPORTED`), a setter outside the objective's row of
   the matrix, data the objective requires and per-asset bounds whose size is not the weight count, cash included
   (`PO_INVALID_ARGUMENT`). Each is an error status, not a result. `solve` then dispatches on the objective.
3. **Translation to a QP.** `mean_variance_problem` (`src/problems/mean_variance.cpp`) takes the risky weights as
   the QP variables, with G = a Sigma (passed with its factor sqrt(a) L) and g = -b mu~; `(a, b)` is (2, 0) for
   `MIN_VARIANCE`, (delta, 1) for `UTILITY`, (2 (1 - lambda), lambda) for `IMPORTANCE` and (1 / (1 + rf), 1) for
   `KELLY`. `tracking_error_problem` builds G = 2 R~'R~ (passed with the factor sqrt(2) R_qr' of a Householder QR
   of R~) and g = -2 R~'(b - rf 1) from the return series, and `solve_max_return` does not build one QP but `Q(t)`,
   with G = Sigma and g = -t (mu~ - c 1), c the midrange of mu~ where the budget row holds 1'w = 1 (0 with cash),
   solved at several `t` (`src/problems/parametric_path.cpp`). Each uses
   `add_constraint_rows`
   (`src/problems/constraints.cpp`) appends the budget row (no cash), the target-return row and one row per finite
   bound, recording in `ConstraintRows` which row is which. With cash, the cash weight 1 - 1'w is eliminated, so its bounds become rows over the risky weights. The
   rows are in [docs/guide/math.md](docs/guide/math.md#the-qp-rows).
4. **Goldfarb-Idnani.** `qp::solve` (`src/qp/goldfarb_idnani.cpp`) uses the factor of G it is given (it factors
   G by Cholesky only when none is given), starts at the unconstrained minimum and adds constraints one at a time:
   the equalities first, then the most violated inequality, each with a primal step bounded by the dual step that drops an active inequality whose multiplier
   reaches zero. The active set keeps the factorization J'N = [R; 0] updated by Givens rotations and never forms an
   inverse. Three numerical decisions are fixed by tests (details in
   [docs/guide/validation.md](docs/guide/validation.md#numerical-decisions-of-the-qp-solver)): a scale-free
   dependence test, tolerances relative to the row and the point, and one step of iterative refinement over the
   active set at the optimum. A violated constraint with no feasible step is `Infeasible`; more than
   10 (n + m) + 10 iterations is `MaxIterations`.
5. **Remapping to the user's variables.** `to_result` (`src/problems/constraints.cpp`) turns the QP solution into
   a `Result`: with cash, the last weight is 1 - 1'w; the budget and target-return multipliers are read from the
   equality rows (NaN when the constraint does not exist); each bound multiplier is the lower-row multiplier minus the
   upper-row multiplier, 0 for a weight with no finite bound. A status other than `Optimal` (for a QP, `Infeasible`
   or `MaxIterations`, whose dual iterate is not primal feasible) gives NaN weights, NaN multipliers and a NaN
   residual, never a partial point. The QP objective equals the user's minimization form up to a constant, so
   the multipliers need no rescaling, with one exception: `MAX_RETURN` divides the multipliers of `Q(t*)` by `t*`
   (`multiplier_scale`) and sets the volatility multiplier to `-s / t*`; for an inactive cap the multipliers are
   the slope `b` of `N b = -mu~`, with scale 1.
6. **Publishing.** The `Result` is moved into a freshly allocated `po_result` and its pointer written to `*out` only
   once everything has succeeded, so a failure leaves `*out` NULL. The Python layer reads the result into a frozen
   dataclass, with read-only weight and multiplier arrays, and frees the handle.

Every `extern "C"` function runs inside `guard`, the single place that turns exceptions into statuses.
