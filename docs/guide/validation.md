# Validation

How PortOpt's results are checked, which tolerances the tests use and the numerical decisions behind them.

## Policy

There is no external oracle. The equations are the author's, and a library that implemented the same equations
would only repeat them, so agreement with it would prove nothing about whether they are the right ones. The tests
instead check each result against conditions that hold if and only if it is correct:

1. **KKT conditions.** For the QP the solver returns, `tests/unit/kkt.hpp` computes the four residuals:
   stationarity `G x + g - E'y - C'z`, primal feasibility, dual feasibility (`z >= 0`) and complementarity
   `z_i (C_i x - c_i)`. For a solved MIN_VARIANCE the same conditions are then checked in the user's variables
   (`expect_user_kkt` in `tests/unit/problems/min_variance_test.cpp`): the equalities and bounds hold, the sign of
   each bound multiplier matches the active side, and the stationarity equation of the
   [multiplier convention](math.md#multiplier-convention) holds. A solution that satisfies the KKT conditions of a
   convex QP is its global optimum.
2. **Closed forms.** Where a closed form exists, the result must match it: the unconstrained minimum and its
   multiplier `2 / A`, the parabola of the efficient frontier, the capital market line with cash, all cash without a
   target return, and a target return equal to the common mean of equal-mean assets.
3. **Consistency.** The same problem is checked at two levels: the QP solution returned by `qp::solve` must satisfy
   KKT in its own variables, and the user-level `Result` built from it must satisfy them in the user's. The C API
   tests solve the same problems through `po_solve` and read the result through the accessors, and the Python tests
   compare the package with closed forms and `numpy`.
4. **Invariants.** Weights sum to one (cash included); the target return is met; a bound multiplier is 0 for a
   weight without a finite bound and its sign matches the active side; an infeasible problem has NaN weights and
   multipliers; a result outlives the spec and data it was solved from; a rejected setter argument leaves the spec
   unchanged; parallel solves over the same read-only handles agree.

## Tolerances

| Check | Tolerance |
|---|---|
| KKT residuals (stationarity, primal, dual, complementarity), relative | 1e-10 |
| KKT residuals of a covariance with condition number 1e10 | 1e-8 |
| Closed forms (unconstrained weights and multiplier, equal-mean target) | 1e-12 |
| Frontier variance and capital market line volatility, relative | 1e-10 |
| Bound feasibility (`w >= lower`, a capped weight at its cap) | 1e-14 |

The KKT residuals are relative: stationarity is divided by `max(1, |G|max |x|max, |g|max)`, primal and dual
feasibility by `max(1, |rhs|max)` and complementarity by `max(1, |z|max) max(1, |rhs|max)`. NaN always fails, so a
solution full of NaN cannot pass a residual check. The random covariances come from `tests/unit/random_covariance.hpp`,
built from a random orthogonal basis and eigenvalues spread over a chosen condition number (100 by default).

## Positive definiteness

A covariance is accepted only if `linalg::cholesky` factors it, and that decision is taken once, in
`Data::create`. The factor rejects a pivot `p_j <= c n eps a_jj`, with `c = 16` (`kPivotFactor` in
`src/linalg/cholesky.hpp`), `eps = DBL_EPSILON` and `a_jj` the original diagonal entry, whatever the sign of
`p_j`. The rounding error of a Cholesky pivot is bounded by a small multiple of `n eps a_jj`, so a pivot below the
floor carries no information: in a singular matrix (a duplicated asset, a sample covariance of fewer periods than
assets) its sign is decided by rounding, and a test against 0 accepted such matrices at random. On sample
covariances of rank T - 1 the largest noise pivot measured was 1.3 `n eps a_jj`; on the test covariances of
condition number 1e10 the smallest pivot is above 4e5 `n eps a_jj`, so `c = 16` separates them with wide margins.
`Data` keeps the factor L, and MIN_VARIANCE hands `sqrt(2) L` to the QP as the factor of `G = 2 Sigma`
(`qp::Problem::hessian_factor`); the solver does not factor G again and does not decide definiteness a second
time. Later objectives with `G = k Sigma` pass `sqrt(k) L`.

## Numerical decisions of the QP solver

Three decisions were taken while validating the Goldfarb-Idnani implementation (`src/qp/`). Each is pinned by a test.

- **A scale-free dependence test.** A constraint can be added only if its normal is independent of the active
  ones. With `d = J'n+` (`J` the basis of the active-set factorization, `n+` the normal being added), the part of `d`
  outside the active columns, `d2`, must satisfy `||d2|| > 1e-12 ||d||`. The test is relative to `||d||`, so it does
  not change when G and g are scaled. A constraint that fails it is redundant (a consistent one is accepted as
  already satisfied, an inconsistent one is `Infeasible`). A target-return row equal to a multiple of the budget row,
  as when every asset has the same mean, goes through this path.
- **Tolerances relative to the point.** A constraint counts as violated when its slack is below
  `-1e-12 max(1, |rhs|, ||row|| |x|max)`; it is never compared with an absolute number.
- **One refinement step at the optimum.** The active-set solution carries a rounding error in the active residuals
  `N'x - b`. One step of iterative refinement over the active set (`refine_toward_active` in
  `src/qp/active_set.hpp`) corrects `x` and the multipliers with the factorization already in hand, so the equality
  constraints and the active bounds hold to the tolerances above without a second solve.

`tests/unit/qp/goldfarb_idnani_test.cpp` solves 50 random budget-and-box problems (2 to 30 variables) and requires
the KKT conditions at 1e-10, at least 50 active bounds in total so the bounds are exercised, a redundant consistent
equality at the objective scale 1e-4, and G and g scaled by 1e-4 and 1e8.

## Where the tests are

| Suite | What it checks |
|---|---|
| `tests/unit/linalg/cholesky_test.cpp` | The factor and the triangular solves; rejection of an indefinite matrix and of a pivot at rounding level |
| `tests/unit/qp/goldfarb_idnani_test.cpp` | Unconstrained, equality-only and bounded problems; infeasible and inconsistent ones; the scale tests above; a provided factor of G gives the same solution |
| `tests/unit/core/data_test.cpp` | Rejection of non-finite, asymmetric, indefinite, singular (a duplicated asset, a sample covariance of fewer periods than assets) and mismatched input; symmetrization within tolerance; the stored factor |
| `tests/unit/problems/min_variance_test.cpp` | KKT in the user's variables, the closed forms, the frontier and the capital market line, per-asset bounds, cash bounds, an ill-conditioned covariance |
| `tests/unit/problems/solve_test.cpp` | Argument checks of `po::solve` and their messages |
| `tests/unit/metrics/portfolio_test.cpp` | Portfolio return and volatility, with and without cash |
| `tests/unit/capi/` | The C boundary: status codes, `po_last_error()`, null handles, ABI constants, the metrics; `c_header_check.c` compiles the header as C11 |
| `python/tests/test_min_variance.py` | The Python package against closed forms, `numpy` and the status contract |
