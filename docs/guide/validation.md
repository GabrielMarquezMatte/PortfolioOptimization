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
| Fuzz battery: KKT residuals and agreement with the reference, `cond(Sigma)` <= 1e4 / 1e6 / 1e8 (and above) | 1e-10 / 1e-9 / 1e-7 |

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

The path `Q(t)` (`src/problems/parametric_path.cpp`, the policy in `src/problems/max_return.cpp`) takes seven
decisions, each pinned by a test in `tests/unit/problems/max_return_test.cpp` or
`tests/unit/fuzz/max_return_fuzz_test.cpp`.

- **The shift of the means.** Where the rows hold `1'w = 1` (no cash), `Q(t)` is solved for `v = mu~ - c 1`, `c` the
  midrange of `mu~`: the same `w(t)`, and the budget multiplier is mapped back by `-c`
  (`MaxReturn.BudgetMultiplierFollowsCommonMean` checks it through `expect_user_kkt` and against means lifted by 1).
  The solver starts from `t Sigma^-1 v`; with the common mean in `g` that start is `t Sigma^-1 mu~`, whose rounding
  swamped `Q(t)` at the `t` (1e12 to 1e14) where nearly equal means at `cond(Sigma) = 1e8` reach the cap: the
  volatility of `Q(t)` jumped by 1 to 5% between adjacent doubles. On the fuzz battery the shift turned 7 of the 8
  such cases from `MAX_ITER` into `OPTIMAL` within their band.
- **The inactive cap: the limit certificate.** The cap is inactive at a point of the climb when the point satisfies the
  KKT conditions of `max mu~'w` over the rows (every row within the QP's feasibility tolerance, no active slope below
  `-16 n eps |v|max`, and `|v + N b|max <= 16 n eps` times its largest term). Those conditions prove the point optimal
  for every larger `t`, so no comparison of weights between two points is needed. On the fuzz battery and the
  `MaxReturn` unit tests (8976 evaluations of the stationarity) the residual of the accepted certificates was at most 0.64 `n eps` of the
  largest term (1033 accepted, 321 of them exactly 0, with equal means) and that of the rejected ones at least
  8.6e4 `n eps`: the factor 16 leaves 25 above the first and 5400 below the second. The former decision, constant
  weights between two doublings within `16 n eps (max(1, |w|max) + t |Sigma^-1 mu~|max)`, was up to `cond(Sigma)`
  too loose: with nearly equal means a slowly moving `w` was declared constant and returned as `OPTIMAL` with a bound
  violated by up to 0.043 or the return off by 1e-5 (15 failing cases of `AgainstFrontierReference` and 29 of
  `NeverSilentlyWrong` before the change). Comparing the weights with the noise of the stationarity equation instead,
  `16 n eps (max(1, |w|max) + t |mu~|max)`, missed truly constant paths: 24 cases with equal means (mostly at
  `cond(Sigma) = 1e8`) climbed until `Q(t)` broke down.
- **The rounding of a weight at `s = sigma_min`.** Only `at_minimum` still compares weights, with the noise of the
  stationarity equation, `16 n eps (max(1, |w|max) + t |v|max)`. On the 320 shapes of `fuzz_shapes(400, 1)` with a
  positive `sigma_min`, solved at `s = sigma_min`, the statuses and the kind of `lambda_vol` (`-inf`, finite, NaN) are
  those of the former noise model case by case (303 `OPTIMAL`: 210 `-inf` and 93 finite; 17 `INFEASIBLE`, where the
  minimum volatility of `Q(0)` differs by more than 1e-12 from that of the `MIN_VARIANCE` solve that set the cap).
- **The multiplier slope.** Where `w` is constant the slope `b` of the QP multipliers solves `N b = -v` over the
  active normals, by Householder least squares; a slope below `-16 n eps |v|max` is a decrease. Differencing two
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
- **Adjacent ends across a breakpoint.** When the bracket shrinks to adjacent doubles on two active sets, the lower end
  is the answer only when `|vol - s| <= 1e-5 s` (`kAdjacentEnd`); otherwise the path has broken down (`MAX_ITER`).
  `vol(t)` is continuous, so at an end that is the answer the deviation is the rounding of `vol(Q(t))`. Measured on
  the fuzz battery and the unit suite: with every bisection forced to adjacent doubles, the 692 ends whose result
  passed every check of `AgainstFrontierReference` deviated by at most 7.2e-8 (at `cond(Sigma) = 1e8`, in line with
  `1e-15 cond(Sigma)`); in the normal search, the one legitimate end across a breakpoint
  (`MaxReturn.EqualMeansWithCash`) deviated by 1.6e-15, and the 9 ends reached where `Q(t)` is lost in its rounding
  (`t` from 1.8e8 to 1.2e12) by 1.6e-3 to 0.23. Those 9 were returned as `OPTIMAL` with `lambda_vol < 0` and a
  volatility up to 23% below `s`. 1e-5 is 140 times the largest rounding and 160 times below the smallest wrong
  deviation.
- **Breakdown.** Feasibility is decided at `t = 0`; a `Q(t)` that fails at `t > 0` is `MAX_ITER`, never
  `INFEASIBLE`.

## Numerical decisions of the QP solver

Four decisions were taken while validating the Goldfarb-Idnani implementation (`src/qp/`). Each is pinned by a test.

- **Feasibility at the final x.** When no row is violated, the refinement step below moves `x` onto the active rows
  by their residuals, which carry the rounding of every step from the unconstrained start. Far from it (a start of
  3e17 in `GoldfarbIdnani.FinalFeasibilityAtLargeScale`: `Q(t)` of `MAX_RETURN` at `t = 1e12` with a nearly riskless
  asset of variance 1e-10) the active rows were off by O(10) before the refinement, and the move left inactive rows
  violated by up to 28. Every row is therefore checked again at the final `x`, with the tolerance of that `x`, and a
  violated row is added like any other. On 48 such problems (12 seeds, `t` from 1e10 to 1e13), the 8 `Optimal`
  answers that violated a bound (by 0.007 to 28) are all feasible with the check (the largest violation 2.5e-4 of the
  tolerance); the others are unchanged.
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

## Fuzz battery

`tests/unit/fuzz/` holds a generator of problem shapes and references computed by another route, used by the
`MaxReturnFuzz` tests.

- **Shapes** (`random_problems.hpp`). `fuzz_shapes(count, seed)` cycles through every combination of `cond(Sigma)` in
  {1e2, 1e4, 1e6, 1e8}, the annual and daily scales (`Sigma` times 1 or 1e-4, `mu` times its square root), four mean
  patterns (spread in [0.03, 0.15]; equal; nearly equal, with a relative spread log-uniform in [1e-8, 1e-5]; below
  the risk-free rate: `rf` above every mean with cash, every mean negative without), four boxes (none, long-only,
  tight with `upper = 2 / n`, and long-only with a third of the weights fixed by `lower == upper`) and cash, so any
  256 consecutive shapes cover them all; `n` (3 to 24) and the covariance seed are random.
- **Reference** (`frontier_reference.hpp`). `max_return_reference` is the largest target return `r` with
  `MIN_VARIANCE(target_return = r)` feasible and of volatility at most `s`, by doubling and bisection on `r` to
  1e-13 relative; it returns the return of that portfolio, not `r`, because the solver accepts a target-return row to
  its feasibility tolerance, with a floor of 1e-12 that is 1e-9 of a daily mean. NaN when no feasible portfolio has
  volatility at most `s`.
- **Bands.** KKT residuals (`expect_user_kkt`) and the agreement of the return (relative to the larger of the
  reference and `|mu|max`) are held to 1e-10 for `cond(Sigma) <= 1e4`, 1e-9 up to 1e6 and 1e-7 up to 1e8, by the
  condition number of the covariance actually solved; above 1e8 the 1e-7 band applies.
- **`MaxReturnFuzz.AgainstFrontierReference`**: `fuzz_shapes(400, 1)` without the fixed box (300 shapes), caps at 0.1,
  0.3, 0.5, 0.7 and 0.9 of the way from `sigma_min` to the largest asset volatility: `OPTIMAL` (or `INFEASIBLE` exactly
  when the reference is NaN), the return of the reference, `expect_user_kkt` and `vol <= s (1 + tol)`, all within the
  band. Of the 1500 cases, 1499 pass, 794 of them with an inactive cap; the worst agreement with the reference is
  1.1e-12 at `cond <= 1e4`, 8.2e-12 at 1e6 and 8.1e-10 at 1e8.
- **`MaxReturnFuzz.NeverSilentlyWrong`**: the nearly-equal-mean shapes of `fuzz_shapes(400, 2)` (100) with a nearly
  riskless asset appended (variance log-uniform in [1e-10, 1e-6], mean 0.4 times the smallest), caps at 0.1, 0.5 and
  0.9 of the way to the largest asset volatility and at 1.01, 2 and 10 times `sigma_min`. The covariance solved has a
  condition number from about 4e2 to 4e16, so most cases take the 1e-7 band. The result must be `OPTIMAL` with the
  return of the reference, every weight within its bounds and `vol <= s (1 + tol)`, or `MAX_ITER`. Of the 528 cases,
  413 are `OPTIMAL` (worst agreement 7.5e-9, at a band of 1e-7) and 115 `MAX_ITER`.
- **The known residual.** One case of `AgainstFrontierReference` stops at `MAX_ITER` and is pinned by name: shape seed
  11394049047870518593 (n 10, `cond(Sigma) = 1e8`, daily scale, nearly equal means, long-only, cash), cap fraction
  0.9. Its cap is reached near `t = 1.2e12`, where `Q(t)` is lost in its rounding (the volatility jumps between
  0.90 s and 1.04 s at adjacent `t`). Before the adjacent-end check it came out `OPTIMAL` at `vol = 0.973 s` with
  `lambda_vol = -8.0e-12`: volatility complementarity 2.7e-2, return 1.6e-9 below the reference. With cash there is no
  budget row and no shift of the means, so nothing removes the common mean from the start `t Sigma^-1 mu~`. Any other
  `MAX_ITER` fails the test.
- **Runtime.** 1.3 s and 0.7 s in `gcc-release`, 15.6 s and 8.8 s in `gcc-debug`; both presets run the full counts.

## Where the tests are

| Suite | What it checks |
|---|---|
| `tests/unit/linalg/cholesky_test.cpp` | The factor and the triangular solves; rejection of an indefinite matrix and of a pivot at rounding level |
| `tests/unit/linalg/qr_test.cpp` | Householder QR: `Q'a = R`, the lower factor of `a'a`, least squares, rank deficiency (duplicated, combined and zero columns), and a rank threshold relative to the columns (`cond = 1e10` full rank, `1e15` not) |
| `tests/unit/qp/goldfarb_idnani_test.cpp` | Unconstrained, equality-only and bounded problems; infeasible and inconsistent ones; the scale tests above; a provided factor of G gives the same solution; feasibility at the final x from a start of 3e17 |
| `tests/unit/core/data_test.cpp` | Rejection of non-finite, asymmetric, indefinite, singular (a duplicated asset, a sample covariance of fewer periods than assets) and mismatched input; symmetrization within tolerance; the stored factor |
| `tests/unit/problems/min_variance_test.cpp` | KKT in the user's variables, the closed forms, the frontier and the capital market line, per-asset bounds, cash bounds, an ill-conditioned covariance |
| `tests/unit/problems/utility_test.cpp` | `UTILITY`, `IMPORTANCE` and `KELLY`: closed forms, defaults, KKT, and the identities between them and `MIN_VARIANCE` |
| `tests/unit/problems/tracking_error_test.cpp` | Exact replication, constrained least squares, cash, target return and long-only KKT, a rank-deficient series rejected, an ill-conditioned one (`cond(R) = 1e8`) accepted with the objective of an independent least squares, no covariance needed |
| `tests/unit/problems/max_return_test.cpp` | Frontier and capital market line, an infeasible cap, the cap at the minimum (including a vertex of the box), an inactive cap, equal means, the multiplier against finite differences, daily scale, random long-only problems, a conditioning sweep (`cond(Sigma)` 1e2 to 1e7), a risk-free rate above every mean (all cash), equal means with tight upper bounds or with cash, the budget multiplier under a common shift of the means |
| `tests/unit/fuzz/max_return_fuzz_test.cpp` | `MAX_RETURN` over the fuzz battery against the frontier reference, and nearly equal means with a nearly riskless asset: correct or `MAX_ITER`, never silently wrong |
| `tests/unit/problems/requirements_test.cpp` | The objective x setter matrix, its messages and the order of the checks |
| `tests/unit/problems/solve_test.cpp` | Argument checks of `po::solve` and their messages |
| `tests/unit/metrics/portfolio_test.cpp` | Portfolio return and volatility, utility and tracking error, with and without cash |
| `tests/unit/capi/` | The C boundary: status codes, `po_last_error()`, null handles, the new setters and `po_data_set_returns`, ABI constants, the metrics; `c_header_check.c` compiles the header as C11 |
| `python/tests/test_min_variance.py`, `test_qp_objectives.py` | The Python package against closed forms, `numpy` and the status contract; the frozen `Spec`, the objective identities, `MAX_RETURN` on the capital market line, the input checks of `Data` |
