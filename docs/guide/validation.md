# Validation

How PortOpt's results are checked, which tolerances the tests use and the numerical decisions behind them.

## Policy

There is no external oracle. The equations are the author's, and a library that implemented the same equations
would only repeat them, so agreement with it would prove nothing about whether they are the right ones. The tests
instead check each result against conditions that hold if and only if it is correct:

1. **KKT conditions.** For the QP the solver returns, `tests/unit/kkt.hpp` computes the four residuals:
   stationarity `G x + g - E'y - C'z`, primal feasibility, dual feasibility (`z >= 0`) and complementarity
   `z_i (C_i x - c_i)`. For every solved objective the same conditions are then checked in the user's variables
   (`expect_user_kkt` in `tests/unit/user_kkt.hpp`, used by the suite of each objective): the equalities and bounds
   hold, the sign of each bound multiplier matches the active side, and the stationarity equation of the
   [multiplier convention](math.md#multiplier-convention) holds, with the gradient of that objective's `f` and, for
   MAX_RETURN, the volatility term. A solution that satisfies the KKT conditions of a
   convex QP is its global optimum.
2. **Closed forms.** Where a closed form exists, the result must match it: the unconstrained minimum and its
   multiplier `2 / A`, the parabola of the efficient frontier, the capital market line with cash, all cash without a
   target return, and a target return equal to the common mean of equal-mean assets. The same for the other
   objectives: the unconstrained `UTILITY` weights and its cash form, the exact replication of a benchmark that is
   a portfolio of the assets (`TRACKING_ERROR`, with zero error), `MAX_RETURN` on the upper branch of the frontier
   and, with cash, on the capital market line with `lambda_vol` equal to minus its slope.
3. **Consistency.** The same problem is checked at two levels: the QP solution returned by `qp::solve` must satisfy
   KKT in its own variables, and the user-level `Result` built from it must satisfy them in the user's. The C API
   tests solve the same problems through `po_solve` and read the result through the accessors, and the Python tests
   compare the package with closed forms and `numpy`.
4. **Consistency between objectives.** The objectives built on one builder must agree with each other, which
   checks the scales `(a, b)` of the builder and not only each objective against its own KKT conditions:
   `IMPORTANCE` with `lambda = 0` is `MIN_VARIANCE` and with `lambda` in (0, 1) is `UTILITY` with
   `delta = 2 (1 - lambda) / lambda`; `KELLY` is `UTILITY` with `delta = 1 / (1 + rf)` (delta = 1 without cash);
   `UTILITY` with a large `delta` approaches `MIN_VARIANCE`; `MAX_RETURN` with the cap at the volatility of the
   unconstrained minimum is `MIN_VARIANCE`, and with equal means too. A `MAX_RETURN` cap below the minimum
   volatility is infeasible.
5. **Derivatives.** A multiplier is the derivative of the optimal `f` with respect to its right-hand side, so the
   test solves again with the right-hand side moved and compares: the volatility multiplier of `MAX_RETURN` against
   a central difference of the optimal `-mu~'w` in the cap `s` (relative tolerance 1e-5, step `1e-6 s`), and, at
   `s = sigma_min` on a vertex of the box, against a one-sided difference of order two, since the optimum is not
   differentiable there from the left.
6. **Invariants.** Weights sum to one (cash included); the target return is met; a bound multiplier is 0 for a
   weight without a finite bound and its sign matches the active side; an infeasible problem has NaN weights and
   multipliers; a result outlives the spec and data it was solved from; a rejected setter argument leaves the spec
   unchanged; parallel solves over the same read-only handles agree.

## Tolerances

| Check | Tolerance |
|---|---|
| KKT residuals (stationarity, primal, dual, complementarity), relative | 1e-10 |
| KKT residuals of a covariance with condition number 1e10 | 1e-8 |
| KKT residuals of `MAX_RETURN` with `cond(Sigma) = 1e7` (`MaxReturn.ConditioningSweep`) | 1e-8 |
| `TRACKING_ERROR` objective with `cond(R) = 1e8` against an independent least squares, relative | 1e-10 |
| Closed forms (unconstrained weights and multiplier, equal-mean target) | 1e-12 |
| Frontier variance and capital market line volatility, relative | 1e-10 |
| Volatility multiplier of `MAX_RETURN` against a finite difference, relative | 1e-5 |
| `vol(w) = s` of a solved `MAX_RETURN`, relative | 1e-12 |
| Bound feasibility (`w >= lower`, a capped weight at its cap) | 1e-14 |

The KKT residuals (`tests/unit/kkt.hpp`) are relative, and the stationarity and dual scales have no floor of 1,
so that a problem whose objective is scaled by 1e-4 is held to the same relative accuracy as one at scale 1:

- **Stationarity** `|G x + g - E'y - C'z|max` is divided by the largest of `|Gx|max`, `|g|max`, `|E'y|max` and
  `|C'z|max`, the four terms that cancel. Only an all-zero set of terms (scale 0, residual 0) is left undivided.
- **Primal feasibility** (the equality gaps and the inequality violations) is divided by `max(1, |rhs|max)`.
- **Dual feasibility** `max(-z_i)` is divided by the multiplier magnitude `|z|max`.
- **Complementarity** is measured per row, as `|z_i s_i| / (max(|z_i|, tiny) max(1, |c_i|, sum_j |C_ij x_j|))`, with
  `s_i = C_i x - c_i` and `tiny` the smallest normal double. For a row with a nonzero multiplier this is the slack
  relative to the magnitude of the row at the point, and for a row whose multiplier is exactly 0 it is 0. The
  solver gives an inactive row a multiplier of exactly 0, so a nonzero multiplier claims an active row, and the
  claim to check is that the row holds with equality to the accuracy rounding allows for that row. A product
  `|z_i s_i|` divided by global scales would let the slack of a row with a small multiplier, or of a row far smaller
  than the largest right-hand side, grow unnoticed.

`expect_user_kkt` holds the user's variables to the same standard, with no floor of 1 either:

- **Stationarity** is divided by the largest magnitude of a single term: each term of the gradient of `f` (passed as
  its terms, such as `delta Sigma w` and `-mu~`, so that a gradient that cancels at an unconstrained optimum is not
  compared with its own rounding) and each constraint term (`lambda_budget 1`, `lambda_return mu~`,
  `lambda_vol Sigma w / vol`, `beta`, `beta_0`).
- **Bound complementarity** is the per-row measure above with the row scale `max(|bound|, |w|max)`: the rounding of a
  weight at its bound comes from the whole solution, so the scale is that of all the weights, not of the one weight
  (which is 0 at a bound of 0).
- **The volatility cap**, when its multiplier is finite: `vol(w) <= s (1 + tol)`, `lambda_vol <= 0`, and
  `|lambda_vol (vol - s)|` relative to `|lambda_vol| s`.

NaN always fails, so a solution full of NaN cannot pass a residual check. With these measures the worst residuals in
the suite are, for the QP, 1.7e-14 (stationarity) at the 1e-10 tolerance and 3.9e-12 for the covariance of condition
number 1e10 at 1e-8; in the user's variables, 4.0e-14 at 1e-10 outside the conditioning sweep of `MAX_RETURN`, and in
that sweep 3.8e-11 up to `cond(Sigma) = 1e6` (at 1e-10) and 2.8e-10 at 1e7 (at 1e-8). The last one is the solver's
own: on the same problems `qp::solve` alone reaches a stationarity of 2.2e-10 along `Q(t)` at 1e7 (5.4e-12 at
`t = 0`), because it starts from the unconstrained minimizer `t Sigma^-1 mu~`. The random covariances come from `tests/unit/random_covariance.hpp`, built from a
random orthogonal basis and eigenvalues spread over a chosen condition number (100 by default).

## Positive definiteness

A covariance is accepted only if `linalg::cholesky` factors it, and that decision is taken once, in
`Data::create`. The factor rejects a pivot `p_j <= c n eps a_jj`, with `c = 16` (`kPivotFactor` in
`src/linalg/cholesky.hpp`), `eps = DBL_EPSILON` and `a_jj` the original diagonal entry, whatever the sign of
`p_j`. The rounding error of a Cholesky pivot is bounded by a small multiple of `n eps a_jj`, so a pivot below the
floor carries no information: in a singular matrix (a duplicated asset, a sample covariance of fewer periods than
assets) its sign is decided by rounding, and a test against 0 accepted such matrices at random. On sample
covariances of rank T - 1 the largest noise pivot measured was 1.3 `n eps a_jj`; on the test covariances of
condition number 1e10 the smallest pivot is above 4e5 `n eps a_jj`, so `c = 16` separates them with wide margins.
`Data` keeps the factor L, and every objective with `G = a Sigma` (the four built by
`src/problems/mean_variance.cpp`, and the path `Q(t)` of `MAX_RETURN`) hands `sqrt(a) L` to the QP as the factor of
`G` (`qp::Problem::hessian_factor`); the solver does not factor G again and does not decide definiteness a second
time. `TRACKING_ERROR` has no covariance: it factors `R~` once by Householder QR (`src/linalg/qr.cpp`), never
`R~'R~`, whose condition number is the square of that of `R~`. The Hessian factor is `sqrt(2) R_qr'`, and the column
rank is decided on the diagonal of `R_qr`: full when every `|R_jj| > c max(T, N) eps max_j |R~_j|`, `c = 16`
(`kRankFactor` in `src/linalg/qr.hpp`), a bound on the rounding of Householder QR relative to the largest column.
With the former Cholesky of `R~'R~`, `cond(R)` above about 1e7 was rejected and accepted problems near that limit lost
up to 1e-3 of the objective; with the QR, 300 random series with `cond(R)` up to 1e8 are all solved, with a largest
relative objective excess of 2.0e-10 against a least-squares reference.

## Numerical decisions of MAX_RETURN

The path `Q(t)` (`src/problems/parametric_path.cpp`, the policy in `src/problems/max_return.cpp`) takes four
decisions, each pinned by a test in `tests/unit/problems/max_return_test.cpp`.

- **The rounding of a weight along the path.** The solver starts from the unconstrained minimizer `t Sigma^-1 mu~`, so
  a weight of `Q(t)` is known to `16 n eps (max(1, |w|max) + t |Sigma^-1 mu~|max)` (`kNoiseFactor = 16`). On
  equal-mean paths, where `w(t)` is constant, the measured change of the weights between doublings stayed below 0.3
  of `n eps (max(1, |w|max) + t |Sigma^-1 mu~|max)` for `n` from 6 to 50 and condition numbers from 1e2 to 1e8, so
  16 leaves a margin of about 50. Two points on the same active set have constant weights when the equality and
  active rows fix `w` (a vertex, decided by their count), or else when the weights differ by no more than that
  rounding. The former fixed tolerance, 1e-14 of the weight scale, is 0 for weights at 0 and is exceeded by the
  rounding as `t` grows, which made a constant path climb until `Q(t)` broke down.
- **The multiplier slope.** Where `w` is constant the slope `b` of the QP multipliers solves `N b = -mu~` over the
  active normals, by Householder least squares; a slope below `-16 n eps |mu~|max` is a decrease. Differencing two
  QP solutions instead gave slopes whose rounding reached 1e6 eps of the multiplier scale at a condition number of
  1e8, and a stationarity of 2.5e-10 for an inactive cap at 1e7.
- **The crossing.** The rounding of the volatility of one QP solution `Q(t*)` was measured at about
  `1e-15 cond(Sigma)` relative (up to 1e-8 at 1e7 and 1e-7 at 1e8), so the former test
  `|vol(w(t*)) - s| <= 1e-12 s` failed on ill-conditioned problems until the bracket shrank to adjacent doubles and
  the search stopped at `MAX_ITER` (4% of long-only problems at 1e4, 50% at 1e7). `Q(t*)` now only confirms that
  `t*` lies on the bracket's segment: the same active set and a volatility within `1e-6` of the bracket's volatility
  range from `s` (`kCrossingSanity`; an error of the affine model is a sizable fraction of the range). The reported
  point is the exact root of the segment, the interpolation of its two ends, whose volatility is `s` to rounding
  (within the 1e-12 of the table, also at 1e7).
- **Breakdown.** Feasibility is decided at `t = 0`; a `Q(t)` that fails at `t > 0` is `MAX_ITER`, never
  `INFEASIBLE`.

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
| `tests/unit/linalg/qr_test.cpp` | Householder QR: `Q'a = R`, the lower factor of `a'a`, least squares, rank deficiency (duplicated, combined and zero columns), and a rank threshold relative to the columns (`cond = 1e10` full rank, `1e15` not) |
| `tests/unit/qp/goldfarb_idnani_test.cpp` | Unconstrained, equality-only and bounded problems; infeasible and inconsistent ones; the scale tests above; a provided factor of G gives the same solution |
| `tests/unit/core/data_test.cpp` | Rejection of non-finite, asymmetric, indefinite, singular (a duplicated asset, a sample covariance of fewer periods than assets) and mismatched input; symmetrization within tolerance; the stored factor |
| `tests/unit/problems/min_variance_test.cpp` | KKT in the user's variables, the closed forms, the frontier and the capital market line, per-asset bounds, cash bounds, an ill-conditioned covariance |
| `tests/unit/problems/utility_test.cpp` | `UTILITY`, `IMPORTANCE` and `KELLY`: closed forms, defaults, KKT, and the identities between them and `MIN_VARIANCE` |
| `tests/unit/problems/tracking_error_test.cpp` | Exact replication, constrained least squares, cash, target return and long-only KKT, a rank-deficient series rejected, an ill-conditioned one (`cond(R) = 1e8`) accepted with the objective of an independent least squares, no covariance needed |
| `tests/unit/problems/max_return_test.cpp` | Frontier and capital market line, an infeasible cap, the cap at the minimum (including a vertex of the box), an inactive cap, equal means, the multiplier against finite differences, daily scale, random long-only problems, a conditioning sweep (`cond(Sigma)` 1e2 to 1e7), a risk-free rate above every mean (all cash), equal means with tight upper bounds or with cash |
| `tests/unit/problems/requirements_test.cpp` | The objective x setter matrix, its messages and the order of the checks |
| `tests/unit/problems/solve_test.cpp` | Argument checks of `po::solve` and their messages |
| `tests/unit/metrics/portfolio_test.cpp` | Portfolio return and volatility, utility and tracking error, with and without cash |
| `tests/unit/capi/` | The C boundary: status codes, `po_last_error()`, null handles, the new setters and `po_data_set_returns`, ABI constants, the metrics; `c_header_check.c` compiles the header as C11 |
| `python/tests/test_min_variance.py`, `test_qp_objectives.py` | The Python package against closed forms, `numpy` and the status contract; the frozen `Spec`, the objective identities, `MAX_RETURN` on the capital market line, the input checks of `Data` |
