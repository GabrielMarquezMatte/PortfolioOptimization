# Validation

How PortOpt's results are checked, which tolerances the tests use and the numerical decisions behind them.

## Policy

There is no external oracle. The equations are the author's, and a library that implemented the same equations
would only repeat them, so agreement with it would prove nothing about whether they are the right ones. The tests
instead check each result against conditions that hold if and only if it is correct:

1. **KKT conditions.** For the QP the solver returns, `tests/unit/kkt.hpp` computes the four residuals:
   stationarity `G x + g - E'y - C'z - sum_j zeta_j (a - s_j)`, primal feasibility, dual feasibility (`z >= 0`,
   `zeta >= 0`) and complementarity `z_i (C_i x - c_i)`, with the active generated rows `(a - s_j)'x >= b` of the
   gross exposure family entering like explicit rows and the family's own violation `b - a'x + ||x||_1` in the
   primal residual. For every solved objective the same conditions are then checked in the user's variables
   (`expect_user_kkt` in `tests/unit/user_kkt.hpp`, used by the suite of each objective): the equalities and bounds
   hold, the sign of each bound multiplier matches the active side, and the stationarity equation of the
   [multiplier convention](math.md#multiplier-convention) holds, with the gradient of that objective's `f` and, for
   MAX_RETURN, the volatility term, and with a leverage the term `lambda_lev g`, `g` in the subdifferential of
   `||w||_1`. A solution that satisfies the KKT conditions of a convex QP is its global optimum.
2. **Closed forms.** Where a closed form exists, the result must match it: the unconstrained minimum and its
   multiplier `2 / A`, the parabola of the efficient frontier, the capital market line with cash, all cash without a
   target return, and a target return equal to the common mean of equal-mean assets. The same for the other
   objectives: the unconstrained `UTILITY` weights and its cash form, the exact replication of a benchmark that is
   a portfolio of the assets (`TRACKING_ERROR`, with zero error), `MAX_RETURN` on the upper branch of the frontier
   and, with cash, on the capital market line with `lambda_vol` equal to minus its slope, and the tangency and the most
   diversified portfolios (`MAX_SHARPE` and `MAX_DIVERSIFICATION` without bounds).
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
   differentiable there from the left. The leverage multiplier against a central difference of the optimal `f` in
   `L` (relative tolerance 1e-6, step `1e-4 L`): on a fixed active set `w(L)` is affine and `f*` quadratic in `L`, so
   the difference is exact up to rounding.
6. **Invariants.** Weights sum to one (cash included); the target return is met; a bound multiplier is 0 for a
   weight without a finite bound and its sign matches the active side; an infeasible problem has NaN weights and
   multipliers; a result outlives the spec and data it was solved from; a rejected setter argument leaves the spec
   unchanged; parallel solves over the same read-only handles agree.

**Where the numbers of this document come from.** The random covariances, means and budgets of the tests come from
`tests/unit/random_covariance.hpp` and `tests/unit/fuzz/random_problems.hpp`. Since plan 8 they use no libm
function and no `<random>` distribution (`po::portable::exp` and `log`, Marsaglia's polar method for the normals,
rejection for the integers), so a test instance has the same bits everywhere: verified locally on Windows (MinGW and MSVC) and on Linux (GCC 14,
glibc 2.28), and the same tests run in CI on all three platforms ([Other platforms](#fuzz-battery)). Every count, worst agreement, residual per `cond(Sigma)` and run time below
that comes from such an instance and that a test or a `DISABLED_Measure*` test of the repository reproduces was
measured again on those instances (`gcc-release`, Windows 10, GCC 16.2 from MSYS2 UCRT64, 2026-11-19, one run of
each, on a machine that was not isolated, so the times are indicative). Figures that came from an instrumented
copy of the code or from a probe that the repository does not keep were measured on the instances of the
earlier generators (`std::log`, `std::exp`, `std::cos` and the `<random>` distributions); they are marked
"earlier instances", they justify a constant and argue a margin, and a change of that constant should measure
them again. The counts that a test asserts are those of this document.

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
| QP with the gross exposure family: value against the sign-pattern brute force, relative; the 1-norm of x at most `L (1 + tol)` | 1e-10; 1e-12 |
| Leverage multiplier against a finite difference, relative | 1e-6 |
| `MAX_SHARPE` against 10 000 feasible portfolios: none with a larger ratio, relative | 1e-12 |
| `MAX_SHARPE` and `MAX_DIVERSIFICATION` weights against the closed form (C++ and Python), absolute | 1e-12 |
| `sharpe_ratio` and `diversification_ratio` of the Python package against `numpy`, absolute | 1e-14 |
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
- **The leverage**: `sum_i |w_i| <= L (1 + tol)` over the risky weights whatever the multiplier (also where it is
  NaN, at a vertical frontier of `MAX_RETURN`), and when the multiplier is finite `lambda_lev <= 0` and
  `|lambda_lev (sum_i |w_i| - L)|` relative to `|lambda_lev| L`. In the stationarity, the
  term is `lambda_lev sign(w_i)` where `w_i` is not 0; where `w_i` is 0, `g_i` may be any value in [-1, 1], so the
  residual there is the part of its magnitude beyond `|lambda_lev|`. A weight within `tol |w|max` of 0 counts as 0:
  the family holds such a weight at 0 only to rounding, and its sign is the rounding's.

The bands assume weights of order 1 (`||w||_1 = O(1)`, as in every test problem). The solver is backward stable:
its rounding is that of a nearby problem, so the forward error of the weights grows with `||w||_1`, and a problem
whose exact weights are large (near the boundary of `MAX_SHARPE` below, with weights of 1e5 and more) is held to its
band relative to those weights, not to 1.

NaN always fails, so a solution full of NaN cannot pass a residual check. With these measures the worst residuals
of the unit tests (the fuzz batteries have their own, below) are, for the QP, 1.5e-14 (stationarity) at the 1e-10
tolerance and 3.3e-11 for the covariance of condition number 1e10 at 1e-8; in the user's variables, 2.4e-13 at
1e-10 outside the conditioning sweep of `MAX_RETURN`, and in that sweep 3.4e-15 up to `cond(Sigma) = 1e6` (at 1e-10)
and 3.7e-15 at 1e7 (at 1e-8). On the earlier instances the sweep reached 3.8e-11 and 2.8e-10, the last one the
solver's own: on the same problems `qp::solve` alone reached a stationarity of 2.2e-10 along `Q(t)` at 1e7
(5.4e-12 at `t = 0`), because it starts from the unconstrained minimizer `t Sigma^-1 mu~`; the instances of the
sweep no longer show it, and the cause stays. The random covariances come from `tests/unit/random_covariance.hpp`,
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
`Data` keeps the factor L, and every objective with `G = a Sigma` (the four built by
`src/problems/mean_variance.cpp`, and the path `Q(t)` of `MAX_RETURN`) hands `sqrt(a) L` to the QP as the factor of
`G` (`qp::Problem::hessian_factor`); the solver does not factor G again and does not decide definiteness a second
time. `TRACKING_ERROR` has no covariance: it factors `R~` once by Householder QR (`src/linalg/qr.cpp`), never
`R~'R~`, whose condition number is the square of that of `R~`. The Hessian factor is `sqrt(2) R_qr'`, and the column
rank is decided on the diagonal of `R_qr`: full when every `|R_jj| > c max(T, N) eps max_j |R~_j|`, `c = 16`
(`kRankFactor` in `src/linalg/qr.hpp`), a bound on the rounding of Householder QR relative to the largest column.
With the former Cholesky of `R~'R~`, `cond(R)` above about 1e7 was rejected and accepted problems near that limit lost
up to 1e-3 of the objective; with the QR, 300 random series with `cond(R)` up to 1e8 are all solved, with a largest
relative objective excess of 2.0e-10 against a least-squares reference (the pivots and the excess of this
paragraph are from the earlier instances).

## Numerical decisions of MAX_RETURN

The path `Q(t)` (`src/problems/parametric_path.cpp`, the policy in `src/problems/max_return.cpp`) takes eight
decisions, each pinned by a test in `tests/unit/problems/max_return_test.cpp` or
`tests/unit/fuzz/max_return_fuzz_test.cpp`.

- **The shift along a face.** A face is a row `r'w >= rho` (or `= rho`) whose normal has every entry -1 or 1, so that it
  holds the whole portfolio: the budget `1'w = 1`, a cash bound, or a generated row `s'w <= L` of the leverage.
  `Q(t)` is solved for `v = mu~ - c r` with the face held as an equality, `c` the midrange of the `r_i mu~_i` (which
  minimizes `|v|max`): on the face `-t v'w = -t mu~'w + t c rho`, so the optimum restricted to the face is the same,
  and its face multiplier in `Q(t)` is that of `v` less `t c`. The solver starts from `t Sigma^-1 v`; with the common
  mean in `g` that start is `t Sigma^-1 mu~`, whose rounding swamps `Q(t)` where nearly equal means reach the cap at a
  large `t`.
  - Where the rows hold `1'w = 1` (no cash) the budget is the face of every point, with no certificate: every
    feasible `w` lies on it. The volatility of `Q(t)` jumped by 1 to 5% between adjacent doubles at `t` from 1e12
    to 1e14 without it; on the fuzz battery it turned 7 of the 8 such cases from `MAX_ITER` into `OPTIMAL` within
    their band (earlier instances). `MaxReturn.BudgetMultiplierFollowsCommonMean` checks the mapping through `expect_user_kkt` and
    against means lifted by 1.
  - Without a budget (cash) the common mean falls on the inequality face that holds the portfolio once it binds: the
    cash bound `w_cash >= 0`, or the leverage, which with cash acts as a budget. The path first tries the face it
    certified last, then the active face of largest multiplier in `Q(t)` itself, and keeps `Q(t)` unshifted when the
    certificate fails (next decision). Without the shift, in `MaxReturnFuzz.LeverageFuzz` the leverage binds at `t`
    near 1 and the spread of the means (1e-8 to 1e-5 relative) moves `w` along it until the cap at `t` from 1e4 to
    beyond 1e15: on the earlier instances 41 cases of 15 shapes, all cash with nearly equal means, failed (25 `MAX_ITER`, 16 accepted at an
    adjacent end with `|vol - s| / s` from 3.4e-10 to 6e-6, above their band; at `t = 2.6e10` on one shape
    `vol(Q(t))` moved by 1.9e-5 between `t` 2e-9 apart). With it, one case remained (the known residual below),
    solved since by the tangent crossing.
  - The multipliers of every point are those of `Q(t)` itself whatever its shift, so that points solved with
    different shifts compare and interpolate; the slope `b_v` of a point is kept in its own frame (`N b_v = -v`) and
    mapped only where the user's multipliers are formed, since `b = b_v - c` at the face adds a common mean of order
    `|mu~|` to entries as small as the spread, whose round trip lost the stationarity of the limit certificate.
    `MaxReturn.FaceShiftLeverageMultiplier` (`lambda_lev = -mu~_best` at the vertex `L e_best`) and
    `MaxReturn.FaceShiftCashMultiplier` (`beta_0 = mu~_best`, `beta_i = mu~_best - mu~_i`) check the mapping in closed
    form and through `expect_user_kkt`; both fail when the face multiplier or the slope is not mapped.
  - The rule of the midrange cannot shift along a face kept as an inequality: `c` above the face multiplier makes
    the shifted optimum leave the face. On the first measured shape the face multiplier was 0.1158017126 and the
    midrange 0.1158017167, and the shifted solution left the face (volatility 1.02 against 5.35) at every `t` tried.
    Holding the face as an equality makes any `c` exact and moves the decision to the sign of the face multiplier.
- **The face certificate.** The optimum restricted to an inequality face solves `Q(t)` exactly when its face
  multiplier in `Q(t)` is nonnegative: then the KKT conditions of `Q(t)` hold with the face among the active rows. The
  decision is the sign itself, with no tolerance. Over the `MaxReturn` unit tests and the three `MaxReturnFuzz`
  tests (23 494 shifted solves of an inequality face, 434 with `c = 0`): 21 726 accepted and 1452 rejected, with 316
  shifted solves not `Optimal` (fallback as well). Relative to `|c|` (`c` is negative for some faces, the leverage
  face among them, so the scale is its magnitude), the face multiplier was within rounding of 0 (`|u / t| <= 1.4e-14
  |c|`, a face whose multiplier is 0 in exact arithmetic) in 3 accepted and 1 rejected solves; both outcomes are
  exact to rounding there, the rejected one falling back to `Q(t)` unshifted. Above that, up to 1e-5 `|c|`, there are
  28 accepted solves, from 4.2e-14 to 7.6e-6 `|c|` (from 4.4e-13 on each is twice the one before: the bisection of `t`
  closing on the exit of the face of one case, which halves the multiplier at every step), and one rejected, at
  -8.8e-6 `|c|`, eight orders above the rounding. On the earlier
  instances the smallest accepted multiplier outside the rounding was 9.0e-4 `|c|` and the rejected one closest to 0
  -5.0e-5 `|c|`.
- **The inactive cap: the limit certificate.** The cap is inactive at a point of the climb when the point satisfies the
  KKT conditions of `max mu~'w` over the rows (every row within the QP's feasibility tolerance, no active slope below
  `-16 n eps |v|max`, and `|v + N b|max <= 16 n eps` times its largest term). Those conditions prove the point optimal
  for every larger `t`, so no comparison of weights between two points is needed. The generated rows of the leverage
  enter as ordinary active rows (their normals `a - s_j` in `N`, and the family checked by the rule of its generated
  row at `sign(w)`), and the certificate is judged in the frame of the point's shift. Over the `MaxReturn` unit tests and the three `MaxReturnFuzz` tests (23 242 evaluations
  of the stationarity, 1234 of them with every term 0, equal means): without generated rows, the residual of the
  accepted certificates was at most 0.99 `n eps` of the largest term (2504 accepted) and that of the rejected ones at
  least 7.4e4 `n eps` (18 589); with generated rows active, at most 2.0 `n eps` (434) and at least 3.4e9 `n eps`
  (1715). The factor 16 leaves 16 above the largest accepted and 4.6e3 below the smallest rejected without them, 8.0
  and 2.1e8 with them. The former decision (earlier instances), constant
  weights between two doublings within `16 n eps (max(1, |w|max) + t |Sigma^-1 mu~|max)`, was up to `cond(Sigma)`
  too loose: with nearly equal means a slowly moving `w` was declared constant and returned as `OPTIMAL` with a bound
  violated by up to 0.043 or the return off by 1e-5 (15 failing cases of `AgainstFrontierReference` and 29 of
  `NeverSilentlyWrong` before the change). Comparing the weights with the noise of the stationarity equation instead,
  `16 n eps (max(1, |w|max) + t |mu~|max)`, missed truly constant paths: 24 cases with equal means (mostly at
  `cond(Sigma) = 1e8`) climbed until `Q(t)` broke down.
- **The rounding of a weight at `s = sigma_min`.** Only `at_minimum` still compares weights, with the noise of the
  stationarity equation, `16 n eps (max(1, |w|max) + t |v|max)`. On the earlier instances, the 320 shapes of `fuzz_shapes(400, 1)` with a
  positive `sigma_min`, solved at `s = sigma_min`, the statuses and the kind of `lambda_vol` (`-inf`, finite, NaN) are
  those of the former noise model case by case (303 `OPTIMAL`: 210 `-inf` and 93 finite; 17 `INFEASIBLE`, where the
  minimum volatility of `Q(0)` differs by more than 1e-12 from that of the `MIN_VARIANCE` solve that set the cap).
- **The multiplier slope.** Where `w` is constant the slope `b` of the QP multipliers solves `N b = -v` over the
  active normals, by Householder least squares; a slope below `-16 n eps |v|max` is a decrease. For the generated
  rows (52 672 slopes over the same runs, 13 403 of them decreasing), the decreasing ones were at most -1.9e7 `n eps
  |v|max` and the others at least -0.30: 1.2e6 and 53 from the threshold. Differencing two
  QP solutions instead gave slopes whose rounding reached 1e6 eps of the multiplier scale at a condition number of
  1e8, and a stationarity of 2.5e-10 for an inactive cap at 1e7.
- **The crossing.** The rounding of the volatility of one QP solution `Q(t*)` was measured (earlier instances) at about
  `1e-15 cond(Sigma)` relative (up to 1e-8 at 1e7 and 1e-7 at 1e8), so a test `|vol(w(t*)) - s| <= 1e-12 s` on a QP
  solution fails on ill-conditioned problems (4% of long-only problems at 1e4, 50% at 1e7, ending `MAX_ITER`). The
  Goldfarb-Idnani solution of `Q(t)` starts from `t Sigma^-1 v`, and on the known residual below that is 6.9e8, whose
  `eps` is the 1.5e-7 of noise in `vol(Q(t))`; the error lies along the path (the solution of a `t` wrong by 1.8e-7),
  not off the face. The crossing therefore comes from the exact tangent of the face, anchored at a solve of the face
  and certified by the KKT conditions of `Q(t*)` at one point (docs/guide/math.md, and
  `docs/research/2026-11-18-max-return-path-and-portable-generators.md`). Its volatility is `s` to the rounding of
  its evaluation: on the known shape, under both covariances, `|vol / s - 1| <= 5.6e-16` at every cap.
- **History: interpolation and adjacent ends.** Until the tangent crossing, the reported point interpolated the two
  QP solutions at the ends of a bracket on one active set, confirmed by a solve of `Q(t*)` within `1e-6` of the
  bracket's volatility range (`kCrossingSanity`); at adjacent doubles across a breakpoint the lower end was accepted
  within `1e-5 s` (`kAdjacentEnd`), a constant that no noise-relative form could replace (the measurements are in this
  file at commit 7287d84). Both are gone: the tangent does not difference two QP solutions, whose rounding made the
  direction wrong by 1e6 to 1e7 relative near the crossing, and an end is accepted only through a certified crossing.
  Identifying a segment by the sign vectors of the generated rows also went: on the known residual the solver
  alternates between two sets of sign vectors of the same face, so no bracket ever held one "active set".
- **Breakdown.** Feasibility is decided at `t = 0`; a `Q(t)` that fails at `t > 0` is `MAX_ITER`, never
  `INFEASIBLE`.

## Numerical decisions of the QP solver

Five decisions were taken while validating the Goldfarb-Idnani implementation (`src/qp/`). Each is pinned by a test;
the counts of problems that failed before each decision were measured on the earlier instances.

- **Feasibility at the final x.** When no row is violated, the refinement step below moves `x` onto the active rows
  by their residuals, which carry the rounding of every step from the unconstrained start. Far from it (a start of
  3e17 in `GoldfarbIdnani.FinalFeasibilityAtLargeScale`: `Q(t)` of `MAX_RETURN` at `t = 1e12` with a nearly riskless
  asset of variance 1e-10) the active rows were off by O(10) before the refinement, and the move left inactive rows
  violated by up to 28. Every row is therefore checked again at the final `x`, with the tolerance of that `x`, and a
  violated row is added like any other. On 48 such problems (12 seeds, `t` from 1e10 to 1e13), the 8 `Optimal`
  answers that violated a bound (by 0.007 to 28) are all feasible with the check (the largest violation 2.5e-4 of the
  tolerance); the others are unchanged.
- **A generated row equal to an active row.** The family never generates a row whose normal and right-hand side equal
  those of an active explicit row, as it never generates one of an active sign vector. The MAX_RETURN path holds a
  face of the family as the equality `-s'x = -L`; the steps from a start far away (`t = 1e10`) leave on that equality
  a residual from 1e-11 to 1e-4, far above the tolerance of the final `x` (about 1.5e-12), and the family,
  evaluated at `sign(x) = s`, generated the equality's own row: dependent on the active set and violated, it was
  declared infeasible (9 of the 25 restricted solves of one shape probed, `t` from 1e9 to 3e11). With the rule, the
  failures left by the face shift went from 16 on 7 shapes to 3 on 2 (`GoldfarbIdnani.FamilyRowHeldAsEquality`). An
  explicit inequality face is moved to the equality rows instead of copied, so no explicit row is ever duplicated;
  that took the 3 to 1.
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

## The gross exposure family

The leverage `sum_i |w_i| <= L` is the family `a'x - ||x||_1 >= b` of the rows `(a - s)'x >= b` over every `s` in
`{-1, 1}^n` (`qp::GrossExposure`; the user's leverage is `a = 0`, `b = -L` over the risky weights). The solver never
lists the `2^n` rows; it generates the one it needs. The decisions, each measured with an instrumented copy of the
solver over `GoldfarbIdnani.GrossExposure*` without `GrossExposureFamilyWithManyExchanges`, the `Leverage` suite and
`LeverageFuzz.AgainstSignPatterns` (1981 solves with the family), and over 54 larger solves outside the battery (earlier
instances; n = 50, 150 and 300; `cond(Sigma)` = 1e2, 1e6 and 1e8; `MIN_VARIANCE` and `UTILITY`; L = 1.1, 1.5 and 3;
all `OPTIMAL` and within their band by `expect_user_kkt`):

- **The generated row and its violation.** `s'x` is largest at `s = sign(x)`, so `(a - s)'x >= b` with that `s` is the
  most violated row of the family, and its slack `a'x - ||x||_1 - b` is the family's. It competes with the most
  violated explicit inequality by the rule of an explicit row, applied to the generated row: violated when the slack
  is below `-1e-12 max(1, |b|, ||a - s|| |x|max)`; the explicit row wins a tie. No new constant. Over the 1992 solves,
  every generated row came from a violation of at least 2.1e6 times that tolerance (5708 generations), while the
  slack that the rule tolerated as rounding was at most 0.12 of it in the middle of a solve and 3.4e-4 of it at the
  final `x`: the rule leaves 8 above the largest rounding and 2.1e6 below the smallest violation.
- **The tie at `x_i = 0`.** `s_i = +1` where `x_i = 0` exactly; the sign of a weight within the tolerance of 0 is that
  of its rounding. Both choices give the same slack, and a sign that is wrong for a weight held at 0 moves the slack
  of the generated row by `2 |x_i|`, which the same tolerance absorbs: no separate tolerance decides what is zero.
  On the earlier instances 208 generations had an exact zero in `x` and 698 a weight within the tolerance of 0. At the final `x` of
  `LeverageFuzz`, `sign(x)` was an active generated row in 81 solves (largest negative family slack 3.0e-4 of the
  tolerance) and differed from every active row, by such a weight, in 50 (largest 2.2e-4). Outside the battery the
  rounding of a weight held at 0 grows with n: at n = 300 and `cond(Sigma) = 1e6` the largest tolerated final slack,
  from such a weight, was 0.84 of the tolerance, and the smallest violation that generated a row 1.1 times it, while
  the slack of `sign(x)` where it is an active row stayed at most 1.0e-3 (at `cond(Sigma) = 1e8`). A weight whose
  rounding passes the tolerance only makes the solver generate the row of the other sign, which can hold the weight
  at 0 together with the first one (a subgradient strictly inside [-1, 1], as in `Leverage.ZeroWeightsSubgradient`);
  the 54 solves were all `OPTIMAL` within their band.
- **Discarded and generated again.** A generated row that leaves the active set, or that was redundant and never
  entered it, is discarded; the family can generate it again later. The active generated rows are numbered after the
  explicit rows, and the final check after the refinement step (above) checks the family too.
- **The iteration limit.** The family counts as 2n rows in the limit `10 (n + rows) + 10`: at most n generated rows
  are active at once, but a generated row can leave and be generated again, and a face of the family held by
  several sign vectors makes the solver alternate between them. With n rows the limit was 490 on the shape of
  `GoldfarbIdnani.GrossExposureFamilyWithManyExchanges` (a solve of `Q(t)` at `t = 9.2e12` of the known shape under
  the Linux covariance) and the solve took 556 iterations; with 2n rows the limit is 730 and, with the portable
  `hypot`, the solve takes 516 (0.71 of the limit, a final stationarity of 4.3e-10 at the 1e-9 the test asks). The
  largest count measured on the rest of the battery was 46 iterations (0.11 of the limit) and, on the earlier
  instances, 1627 (0.27 of the n-row limit) at n = 300, `cond(Sigma) = 1e8`, `L = 1.1`.
- **What the solution reports.** `qp::Solution` has the sign vectors of the active generated rows, in lexicographic
  order, and their multipliers `zeta_j`, with `G x + g = E'y + C'z + sum_j zeta_j (a - s_j)`. The active set is
  identified by the active explicit rows and these sign vectors. In the user's variables
  `lambda_lev = -sum_j zeta_j / scale`: the term `-sum_j zeta_j s_j` is `lambda_lev g` with
  `g = sum_j zeta_j s_j / sum_j zeta_j` in the subdifferential of `||w||_1`. It is NaN without a leverage and 0 when
  no generated row is active. A family whose bound is never reached leaves the solution bit for bit that of the
  problem without it (`GoldfarbIdnani.GrossExposureInactive`).

## Numerical decisions of MAX_SHARPE and MAX_DIVERSIFICATION

The homogeneous QP (`src/problems/homogeneous.cpp`, the equations in
[math.md](math.md#max_sharpe-and-max_diversification)) adds one decision of its own, the noise of `kappa`; the others
are the solver's. Each is measured over `RatioFuzz.AgainstFrontierReference` (1200 cases, 1067 of them with an
`Optimal` QP) with an instrumented copy of the builder, and the noise of `kappa` also at the boundary:

- **UNBOUNDED is the active set.** The row `1'y >= 0` is in the active set of the optimum, as the solver reports it
  (`qp::Solution::active_inequalities`), and the solver adds that row only when it is violated beyond its feasibility
  tolerance. The row was active in 28 of the 1067 solves, every one of them `UNBOUNDED` by the reference, with
  `kappa = 1'y` at most 2.4e-16 of `sqrt(n) |y|max` (the scale of the row at the point); in the other 1039, `kappa` was
  at least 0.145 of it. The row alone does not decide the boundary, below.
- **The noise of kappa: the boundary of its sign.** Without bounds, `kappa = 1' Sigma^-1 v / v' Sigma^-1 v` is fixed by
  the data, and `1' Sigma^-1 v = 0` is the border between `OPTIMAL` (it is positive) and `UNBOUNDED` (it is not). When
  `1' Sigma^-1 v` is zero to rounding, for example `rf` equal to the return of the global minimum-variance portfolio
  (`B / A` in the notation of [math.md](math.md#min_variance)) with unbounded short selling, the QP returns
  `kappa = 1'y` of either sign within its rounding with the row `1'y >= 0` inactive. Before this decision a positive one
  was `OPTIMAL` with weights `y / kappa` of order `1 / eps`: against an exact rational reference, `|w|` from 1e14 to
  3e15 with weights wrong by 0.37 to 2.8 relative at `rf = r_gmv`, `|w| = 5e11` and 83% wrong at `cond(Sigma) = 1e6`,
  `|w|` about 1e9 and 44 to 134% wrong at 1e8, several of them with an exact `kappa <= 0`, which is `UNBOUNDED`.
  - **The model.** Where only `v'y = 1` is active, `Sigma y = sigma_y^2 v` and `kappa = sigma_y^2 z'v` with
    `z = Sigma^-1 1`. A backward error `|dSigma| <= n eps |Sigma|` of the solve and the rounding `|dv| <= eps m` of `v`
    (`m = |mu| + |rf|` for `MAX_SHARPE`, `sigma` for `MAX_DIVERSIFICATION`) move it to first order by at most
    `kappa_noise = n eps (|z|'|Sigma||y| + sigma_y^2 |z|'m)`, with `z` from two triangular solves on the factor of
    `Data` (`kappa_noise` in `src/problems/homogeneous.cpp`). On an unbounded cone, `kappa <= 16 kappa_noise`
    (`kKappaNoiseFactor`) with the row inactive is the closure point `kappa = 0`: `UNBOUNDED`. It replaces the former
    rule `kappa <= 0`, which it contains. The model is a first-order bound in the spirit of Higham, with `|Sigma|` where the
    rigorous backward error of a Cholesky solve is stated with `|R'||R|`, `R` the factor (stated from memory, not verified in
    the book); the margin that the measurements below show (0.034 `kappa_noise` at the boundary) is what the guard rests on.
  - **Only on unbounded cones.** `kappa_noise` grows like `n^2 eps cond(Sigma) kappa`. For box-feasible points
    (`n = 150`, earlier instances), `16 kappa_noise / kappa` has a median of 2e-4 at `cond(Sigma) = 1e8`, 1.5e-2 at 1e10 and 2.0 at 1e12.
    Applied everywhere, the threshold returned `UNBOUNDED` on bounded feasible sets, where it is impossible:
    - at 1e12 and `n = 150`, 5 of 6 boxes [-1, 1];
    - at 1e13, some long-only and leveraged problems;
    - at 1e14, every bounded type.

    At `kappa = 0` the rows read `y_i >= 0` (finite lower), `y_i <= 0` (finite upper), `1'y = 0` and, with a leverage,
    `||y||_1 <= 0`. So the cone is bounded, and `kappa` stays away from 0, exactly when one of two things holds:
    - a leverage is set;
    - no pair `i != j` has `upper_i = +inf` and `lower_j = -inf` (`cone_bounded`).

    There only the row and `kappa <= 0` decide. On the reviewer's construction (earlier instances), with means of either sign, `n` 50 and
    150, `cond(Sigma)` 1e10 to 1e14 and both scales, 240 bounded solves are all `OPTIMAL`:
    - a box, long-only with and without a cap, a leverage alone and with a box;
    - the worst bound violation is 1.1e-15 and the worst `||w||_1 / L - 1` is 8.9e-16.

    The battery and the boundary tests have no bounded cone near `kappa = 0`, and their statuses do not change.
  - **A partially bounded boundary.** Some weights are bounded and the cone is still unbounded
    (`MaxSharpe.PartiallyBoundedAtMinimumVarianceReturn`):
    - free assets at `rf = r_gmv` of their block, with three assets in [0, 1] whose lower rows are active as
      `kappa -> 0`;
    - `cond(Sigma)` 1e2 to 1e8 and both scales;
    - `rf` 0 to 16 ulps away and moved down by 3 to 3e4 `cond(Sigma) eps`;
    - 192 solves against a long double reference.

    The model uses `z = Sigma^-1 1` of the whole covariance, not of the free block, and still holds:
    - every exact `kappa <= 0` is `UNBOUNDED`;
    - every `OPTIMAL` result has an exact `kappa` of at least 36 `kappa_noise`, the lower rows of the bounded assets
      active, and weights whose relative error bounds the error of `kappa` by 1.0e-2 `kappa_noise`;
    - the 56 `UNBOUNDED` results with an exact `kappa > 0` have at most 11.0 `kappa_noise`;
    - none above the threshold is `UNBOUNDED`.
  - **The boundary side.** The error of the computed `kappa` against its exact value, in units of `kappa_noise`, with
    the row `1'y >= 0` inactive, is at most 0.034: on `MaxSharpe.RiskFreeAtMinimumVarianceReturn` at `rf = r_gmv`,
    where the data make `kappa` exactly 0 (16 problems, `cond(Sigma)` 1e2 to 1e8, both scales; at most 0.024), on the
    same test with `rf` 1, 3 and 16 ulps away, where `kappa` is that of the construction (96 solves; at most 0.030), on
    `MaxDiversification.InverseCovarianceTimesVolatilitiesSumsToZero`, also exactly 0 (24 problems; at most 0.016), and
    on the partially bounded boundary above (192 solves; at most 0.034). The reviewer's probes of the earlier
    instances (`n = 6`, `cond(Sigma)` 1e2 to 1e8, both scales, 3 draws each, `rf = r_gmv` computed in floating point,
    at 0 and 1 ulp from it and moved by 3 to 3e4 times `cond(Sigma) eps` relative, 312 solves against an exact
    rational `kappa`; 187 with the row inactive, at most 0.029) are not repeated. The factor 16 is 470 times the
    largest error.
  - **The battery side.** The 1039 `OPTIMAL` solves of `RatioFuzz` have `kappa` at least 3.6e6 `kappa_noise`, 2.3e5
    times the threshold, and the factor changes no status of the battery. Outside the boundary tests, the unbounded
    cones of the unit suites of `MAX_SHARPE` and `MAX_DIVERSIFICATION` were at least 1.1e11 on the earlier instances.
  - **What it gives** (earlier instances). On the probes no result is `OPTIMAL` with an exact `kappa <= 0` (several were before). The 99
    `OPTIMAL` results have an exact `kappa` of at least 23 `kappa_noise`, weights from 11 to 1.3e12 and a relative error
    of the weights of at most 6.5e-4, consistent with an error of `kappa` of at most 0.029 `kappa_noise` over a `kappa`
    above 16: backward-stable results whose band is relative to weights far from order 1 (see "Tolerances"). The other
    side of the same honesty: 54 problems whose exact `kappa` is positive but at most 16 `kappa_noise` are
    `UNBOUNDED`, with exact weights of at least 3.2e6, since their sign is within the rounding of the data.
- **Stationarity above `cond(Sigma)` = 1e8: a conditioning limit** (earlier instances). Where many weights are interior, the relative
  stationarity residual of `expect_user_stationarity` grows about as 0.1 `cond(Sigma) eps`:
  - the box [-1, 1] at `n = 150`: 1.5e-9 to 3e-9 at 1e8, 1.2e-7 to 2.5e-7 at 1e10, 9e-6 to 2e-5 at 1e12 and 7e-4 to
    1.7e-3 at 1e14;
  - the tangency without bounds, the same values, 1.1e-5 to 1.9e-5 at 1e12;
  - a leverage alone: 1.4e-6 to 4.5e-6 at 1e12.

  It is not along `v`, so it does not come from the mapping `nu = 2 sigma_y^2`. It does not depend on the noise
  threshold. It is the limit of a backward-stable solve in double: `(Sigma + dSigma) w = c v` with
  `|dSigma| <= n eps |Sigma|` leaves a residual of up to `n eps |Sigma||w| / |Sigma w|` relative to `|Sigma w|`, and
  that quotient reaches `cond(Sigma)` along the eigenvectors of small eigenvalues. The 1e-7 band holds such problems at 1e8
  and not at 1e10. Long-only problems, whose active bounds dominate the scale, stay at most 5.4e-9 at
  1e12 and 8.3e-8 at 1e14. So above 1e8 the bounded-cone tests (`BoundedConeAtLargeCondition`) check the status, the
  budget, the bounds and the leverage at the 1e-7 band, and the whole KKT only for long-only at 1e12 (margin 16: the worst stationarity of the test is 6.2e-9 at the band of 1e-7; on the earlier instances it was 1.4e-9 at
n = 150 and, at 1e13, 3.9e-8, which would leave 2.6). No
  band is added for the stationarity there.
- **Cash bounds.** With cash the cash weight is held at exactly 0, and the request is `INFEASIBLE` when the cash bounds
  exclude 0: `lower_0 > 0` or `upper_0 < 0`, compared with 0 exactly.
- **The multipliers.** The mapping uses `nu = 2 sigma_y^2`, which complementarity gives in exact arithmetic, rather
  than the QP's `nu`, and computes `sigma_y = |L'y|` from the factor of `Data`. The difference between the two is the
  rounding of the active rows, and enters the user's stationarity as a term along `v`: the worst stationarity of the
  battery is at most 0.075 of its band (below).

## Numerical decisions of RISK_PARITY

Without `target_return`, `target_volatility` or an active bound the solution is the positive-orthant point of Spinu:
the minimizer `y` of `F(y) = 1/2 y' Sigma y - b' log y` over `y > 0`, and `w = y / 1'y`. `F` is strictly convex and
coercive, its gradient is `Sigma y - b / y` and its Hessian `Sigma + diag(b / y^2)`. At the minimizer
`y_i (Sigma y)_i = b_i`, so the risk contributions of `w` are `b` and `sum_i b_i = 1` gives `y' Sigma y = 1`.

**The start.** `y_i = sqrt(b_i / Sigma_ii)`, scaled to `y' Sigma y = 1`. The first factor is the exact minimizer of a
diagonal `Sigma` (`w_i` proportional to `sqrt(b_i) / sigma_i`), the second is the value of `y' Sigma y` at every
minimizer, so the start has the scale of the answer at any scale of `Sigma` (annual or daily).

**The step.** One Cholesky factorization of the Hessian per iteration and two triangular solves (`L z = -g`,
`L' d = z`), never an inverse. Far from the minimizer (`delta = max_i |d_i / y_i| > 1/4`) the step is damped: Armijo
backtracking (`1e-4`, halving) that must also decrease `F` strictly, starting at 0.99 of the way to the boundary of
the orthant, so `y > 0` always. The damped step of the self-concordant analysis is not used: the barrier with
`b_i < 1` is not standard self-concordant. Close to it (`delta <= 1/4`) the full Newton step is taken without testing
`F`: near the minimizer `F` is flat to its rounding while the weights still move, so a test on `F` stops long before
the weights are right.

**The stop.** On the size of the step just taken, `delta = max_i |d_i / y_i|`, which is dimensionless and bounds the
relative change of every weight at once: after a step of size `delta` the relative error of each weight is of order
`delta^2`. It stops when `delta <= 16 eps` (`eps = DBL_EPSILON`), or when `delta` stops halving (`delta >= delta_prev /
2`) at `delta <= 1e-6`, the rounding floor of the step, which depends on `cond(Sigma)` and sits above `16 eps` in about
two thirds of the runs. A test on the Newton decrement `lambda^2 <= f n eps (1 + |F|)` was measured first and
rejected: it bounds `F - F*`, not the weights. The error of weight `i` is at most `lambda / sqrt(b_i)`, and with
`b_i` down to about `1e-7` (uneven budgets at `n = 100`) the relative error of a weight reached `1.6e-7` at `cond(Sigma) = 1e4`
with `f = 16` and the full step after the test, against `1.1e-10` when run to the floor.

**Measurement.** `n` in {1, 2, 5, 10, 20, 50, 100}, `cond(Sigma)` in {1e2, 1e4, 1e6, 1e8}, both scales (`Sigma` and
`1e-4 Sigma`), equal budgets and uneven ones (`b_i` log-uniform in `[1e-6, 1]`, normalized), seeds 1 to 10,
`random_covariance`: 1120 runs, none failed to stop. The reference is Newton on the same `Sigma` in `long double`
run to a step below `1e-17`; the error is `max_i |w_i - w_i*| / w_i*`. Windows, GCC 16.2 (MSYS2 UCRT64).

Newton iterations, mean (max) over `cond(Sigma)` up to 1e8:

| n | equal budget | uneven budget |
|---|---|---|
| 1 | 1.0 (1) | 1.0 (1) |
| 2 | 1.4 (4) | 10.5 (19) |
| 5 | 9.0 (18) | 13.2 (21) |
| 10 | 10.4 (17) | 16.8 (28) |
| 20 | 10.7 (16) | 19.4 (32) |
| 50 | 11.4 (16) | 23.6 (35) |
| 100 | 12.5 (17) | 29.6 (53) |

The ceiling is 200 iterations, 3.8 times the largest count (53).

The constant `16` in `delta <= 16 eps` (runs ending by that test rather than by the floor, mean iterations, largest
relative error of a weight at `cond(Sigma) = 1e4` and `1e8`; the 960 runs with `n >= 2`):

| factor | by the test | mean iterations | error at 1e4 | error at 1e8 |
|---|---|---|---|---|
| 1 | 152 | 14.29 | 1.3e-12 | 1.1e-9 |
| 16 | 347 | 14.04 | 1.3e-12 | 1.1e-9 |
| 256 | 559 | 13.74 | 1.3e-12 | 1.1e-9 |
| 4096 | 751 | 13.42 | 1.4e-12 | 1.1e-9 |

The error does not depend on the constant up to 4096, because it is the floor of the step that decides; the constant
only trades the iterations spent at the floor (0.9 per run between 1 and 4096). 16 is kept as the smallest multiple
of `eps` that is above the rounding of a step in the well-conditioned runs. At the stall the largest step size
measured was `3.2e-9`, so the `1e-8` ceiling on the stall exit is above every measured stall (by a factor of 3.2)
without accepting a step that is still converging. The ceiling was `1e-6` before; it was lowered because the Spinu
point is the one `OPTIMAL` answer that does not pass the final KKT gate, and with `1e-8` every test outcome is
unchanged.

Largest relative error of a weight against the long double reference (at `n` up to 100), against the bands of the
tests:

| cond(Sigma) | weight, measured | band |
|---|---|---|
| 1e2 | 1.9e-13 | 1e-10 |
| 1e4 | 1.3e-12 | 1e-10 |
| 1e6 | 3.7e-11 | 1e-9 |
| 1e8 | 1.1e-9 | 1e-7 |

The error grows with `cond(Sigma)` as the conditioning of the weights does, not with the number of iterations. For a
diagonal `Sigma` the start is the exact answer and one step confirms it (relative error `1e-15`).

**The scale of the stationarity check.** The check of `grad f = budget 1 + ... + beta` is relative to the largest
term of the gradient and of the multipliers. At a point with `f = 0` (the Spinu solution) every term of the gradient
is proportional to the gap `c - b`, so they vanish together and that scale is rounding. There, and only there (when
`max_i |c_i - b_i|` is at most the tolerance of the check), the scale is the largest `|grad f_j|` a gap of 1 can
produce, so the check says the gap is below the tolerance times that unit. At any larger gap the terms are the scale
and the band is not relaxed: a stationarity error of `1e-8` of the gradient at a gap of `1e-3` is rejected at
tolerance `1e-10` (`RiskParityReference.KktRejectsAStationarityErrorAwayFromTheFloor`).

### Derivatives of f

`linearize` (`src/riskparity/objective.cpp`) returns `f`, `grad f = 2 J'r`, the Jacobian
`J = (diag(p) + diag(w) Sigma - 2 c p') / v` (`p = Sigma w`, `v = w'p`) and, on request, the exact Hessian
`2 J'J + 2 sum_i r_i hess(c_i)`, with
`sum_i r_i hess(c_i) = (diag(r) Sigma + Sigma diag(r)) / v - 2 (a p' + p a') / v^2 - 2 k Sigma / v + 8 k p p' / v^2`,
`a = r o p + Sigma (w o r)` and `k = c'r`. Measured on 25 points (`n` in {2, 3, 7, 10, 30}, `cond(Sigma) = 1e4`,
weights with a short position, budgets in `[0.2, 1]` normalized):

| quantity | reference | largest error | band of the test |
|---|---|---|---|
| gradient | three-term formula in `long double` (`risk_parity_gradient`) | 1.7e-15 of max abs gradient | 1e-12 |
| Jacobian | central differences of `c`, step `1e-6` | 7.0e-11 of max abs J | 1e-6 |
| Hessian | central differences of the reference gradient, step `1e-6` | 6.6e-10 of max abs H | 1e-6 |
| `J w = 0` | homogeneity of degree 0 | 5.1e-16 of `sum_j abs(J_ij w_j)` | 1e-12 |
| `value_change` | `g's + 1/2 s'Hs`, abs(s) from `1e-8` to `1e-15` (n = 7) | 9.3e-16 of `abs(g's)` | 1e-12 |

`value_change` computes `f(w + s) - f(w)` from `s` (`dc_i = (dN_i v - N_i dv) / (v (v + dv))` with
`dN_i = s_i p_i + w_i q_i + s_i q_i`, `q = Sigma s`), so it keeps its relative accuracy at any size of `s`. A
difference of two values of `f` loses every digit below about `eps f`. The steps of the test start at `1e-8`
because at `abs(s) = 1e-6` the cubic remainder of the expansion was already `9.1e-10` of `abs(g's)`: the band
would have measured truncation, not rounding.

**`Sigma w` by a compensated product.** Each entry of `p = Sigma w` (and of `Sigma s`) is the Dot2 of Ogita, Rump
and Oishi (2005, "Accurate sum and dot product", SIAM J. Sci. Comput. 26(6)): every product and partial sum is
split exactly into its rounded value and its error (TwoProduct, TwoSum), and the errors are added
back once, so the result is as accurate as a sum in twice the working precision, rounded once. Plain `double`,
portable to every compiler; the accuracy does not come from an extended type (`long double` is `double` on MSVC).
TwoProduct is Dekker's product (Veltkamp's split and four exact products) where it is exact, both factors of
magnitude in `[2^-485, 2^500]` or zero, and `std::fma(a, b, -fl(a b))` elsewhere: the two give the same bits,
and the first needs no call into `libm` in a build without FMA instructions.
The entries of `Sigma w` cancel by up to `cond(Sigma)`, and an error of `eps cond(Sigma)` in `p` reaches
`grad f` at the same relative size: summed plainly in `double`, converged descents (earlier instances) at `cond(Sigma) = 1e6` and
`1e8` (LongOnly, uneven budgets, `n = 10`) failed the stationarity of `expect_risk_parity_kkt` by `1.4e-8` and
`1.1e-5` against bands of `1e-9` and `1e-7`, with the descent's own residual below them. With the compensated
product none of the converged runs of the sweep below fails, at any `cond(Sigma)`. `f` differs from
`metrics::risk_budget_residual`, which sums plainly, by the rounding of the metric (`2.4e-15` relative measured, earlier instances).

The error-free transformations hold only if every operation is rounded as written. `portopt_core` is compiled with
`-ffp-contract=off` on GCC and Clang: with GCC 16.2 at `-O3 -std=c++23`, `a * b + c` still becomes one `vfmadd`
when the target has FMA (`-mfma`, `-march=native`), which would fold a product into the following TwoSum; with
`-ffp-contract=off` it stays `vmulsd` and `vaddsd`. No build uses `-ffast-math` or `-fassociative-math`
(disabled at `-O3`). With Visual Studio 2022 and later, MSVC's default `/fp:precise` neither contracts nor
reassociates (Microsoft's documentation; the MSVC build runs the same tests, and `Determinism.*` and
`Generators.HaveTheCommittedBits` check its bits against the committed ones).

### The polyhedral descent

`descend` (`src/riskparity/polyhedral_descent.cpp`) minimizes `L = f + L3` over the polyhedron of the budget, the
target-return row and the bounds by proximal Gauss-Newton steps, each a QP of `qp::solve` with
`Q = 2 J'J + tau I (+ rho grad h3 grad h3')`, then Armijo (`1e-4`, halving, at most 60) along the segment to the
QP solution. Three numerical decisions, each forced by a measurement (earlier instances):

- **The QP is posed in the step `d`**, with the rows shifted to the iterate (`rhs - N w`), not in `x = w + d` with
  the linear term `grad L - Q w`. In `x` the solver returns `d` as the difference of two points of the size of `w`,
  rounded at `eps cond(Q) abs(w)`: on a problem whose solution has `f = 0` (`n = 3`, `cond(Sigma) = 1e4`, no
  bounds), the descent in `x` stalled at `f = 6.5e-23` with a relative stationarity of 1.19; posed in `d` it
  reached `f = 2.5e-26` and stopped converged in 6 iterations. The multipliers are the same in both.
- **Armijo is tested on the Lagrangian `l = L - u'(N w - rhs)` of the step's multipliers**, with slope
  `-d'Qd`, not on `L` with slope `grad L'd`. On the polyhedron the change of `L` is at most that of `l` (the
  equality part of `u'N s` is 0 and the active-bound part is at most 0), so every accepted step still decreases
  `L`. Near a point with large multipliers the rounding of the step across the rows moves `L` by `u'N s`, more
  than the decrease: with the test on `L` the descent stopped at a relative stationarity of `3.3e-9` (`n = 5`,
  three weights at a cap of 0.24, uneven budget); with the test on `l` it reached `3.2e-13` in 9 iterations.
- **The stop** is the relative KKT residual at the iterate, the largest of the stationarity
  `max abs(grad L - N'u)` (relative to the largest term of `grad f`, of the volatility term and of `N'u`, and,
  when `sqrt(f)` is at most the tolerance, to the unit-gap scale `2 max_j sum_i abs(J_ij)` of the section above)
  and the slack of every row with a nonzero multiplier relative to `max(abs(rhs), max abs(w))`. Not the change of
  `w` or of `f`, which stop at points that are not stationary.

**Monotone `L`** (earlier instances). In exact arithmetic `L` decreases at every step. The computed `L` can still grow once the
decrease is below the rounding of the iterate across the rows, which moves `L` by up to about
`eps sum_k abs(u_k) sum_j abs(w_j)`. Measured over 60 runs (`n = 8`, box `[0, 0.15]`, uneven budget, `cond(Sigma)
= 1e2`, seeds 1 to 30, with and without a volatility term): 1325 steps, 129 computed increases, the largest 0.20
of that unit. `PolyhedralDescent.MonotoneDescent` allows exactly one unit.

**The proximal weight.** `J` does not change when `Sigma` is scaled (`c` is homogeneous of degree 0 in
`Sigma`), so the weight of the proximal term is set in the units of `2 J'J`, with the Levenberg-Marquardt update:

- `tau_0 = proximal tr(2 J'J) / n` at the start, `proximal = 0.05`;
- after an iteration whose full step (`gamma = 1`) Armijo accepts, `tau /= 3`; after one it shortens, `tau *= 4`;
- at every iteration `tau >= 100 n eps max_i (2 J'J)_ii`, above the rounding of `2 J'J` (`linalg::cholesky`
  rejects a pivot at or below `16 n eps a_jj`);
- `Q` is factored by `linalg::cholesky` before the QP, which receives the factor; while it is rejected,
  `tau *= 10`, at most 8 times, after which the descent stops with `NotPositiveDefinite` rather than throwing.

Without a volatility term every decision depends on `Sigma` only through `J`, so daily and annual units give the
same iterates (`PolyhedralDescent.ScalingSigmaLeavesIteratesUnchanged`, `Sigma` and `1e-4 Sigma`: same stop,
same iteration count, weights within `1e-10`).

The first rule measured (earlier instances), `tau = 0.05 tr(Sigma) / (2n)` (the `riskParityPortfolio` package,
`riskParityPortfolioSCA`), is not in those units. On a sweep like the one below (2400 runs, Box None with
`target_return` included, tolerance `1e-12`): at scale 1 and `cond(Sigma) = 1e8`, `tau` dominated `2 J'J`, every
step was a short gradient step and 22 of 300 runs converged within 2000 iterations; at scale `1e-4` and
`cond(Sigma) = 1e8`, `tau` fell below the rounding of `2 J'J` and `qp::solve` threw on 16 of 300 runs. In total
1640 of 2400 converged. That is the rule of the plan, `tau = 0.05 tr(Sigma) / (2n)`, and the Levenberg-Marquardt `tau`
above replaces it. A fixed scale-free `tau = 0.05 tr(2 J'J) / n` converged on 1470 (1568 with 0.005): the
adaptive update, not only the units, is what makes the descent converge.

**Sweep** (`RiskParityDescentFuzz.DISABLED_MeasureSweep`, run with `--gtest_also_run_disabled_tests`; 139 s for the two values of `proximal` below).
Every combination of `cond(Sigma)` in {1e2, 1e4, 1e6, 1e8}, scale {1, 1e-4}, equal and uneven budgets, box None,
LongOnly and Capped, with and without `target_return`, except Box None with `target_return`: that is an unbounded
line on which `f` can decrease towards its value at infinity, the contributions of the direction of the line (on
one `n = 3` case the iterates moved past `abs(w_i) = 1` while a scan over `abs(w_i) <= 5` found `f = 0.095`). `n`
in {3, 5, 10, 20, 50}, seeds 1 to 5: 2000 runs, started at the minimum-variance point of the polyhedron, at most
2000 iterations, tolerance a tenth of the band of `cond(Sigma)`, every converged run checked by
`expect_risk_parity_kkt` at the band. The tenth: the descent's KKT residual and that of `expect_risk_parity_kkt`
agree to a few percent, and a stop at the band itself landed at `1.03e-10` against `1e-10` on 4 runs.
`proximal = 0.05`, with the converged count at `0.005` for comparison:

| cond(Sigma) | scale | box | converged | cap reached | no decrease | KKT failed | mean (max) iterations | converged, 0.005 |
|---|---|---|---|---|---|---|---|---|
| 1e2 | 1 | None | 50 / 50 | 0 | 0 | 0 | 10.4 (16) | 50 |
| 1e2 | 1 | LongOnly | 99 / 100 | 1 | 0 | 0 | 13.4 (84) | 98 |
| 1e2 | 1 | Capped | 100 / 100 | 0 | 0 | 0 | 20.1 (119) | 100 |
| 1e2 | 1e-4 | None | 50 / 50 | 0 | 0 | 0 | 10.4 (16) | 50 |
| 1e2 | 1e-4 | LongOnly | 99 / 100 | 1 | 0 | 0 | 13.3 (84) | 98 |
| 1e2 | 1e-4 | Capped | 100 / 100 | 0 | 0 | 0 | 20.1 (119) | 100 |
| 1e4 | 1 | None | 50 / 50 | 0 | 0 | 0 | 16.8 (23) | 50 |
| 1e4 | 1 | LongOnly | 92 / 100 | 4 | 4 | 0 | 15.6 (89) | 93 |
| 1e4 | 1 | Capped | 100 / 100 | 0 | 0 | 0 | 21.0 (216) | 100 |
| 1e4 | 1e-4 | None | 50 / 50 | 0 | 0 | 0 | 16.8 (23) | 50 |
| 1e4 | 1e-4 | LongOnly | 93 / 100 | 3 | 4 | 0 | 16.8 (120) | 93 |
| 1e4 | 1e-4 | Capped | 100 / 100 | 0 | 0 | 0 | 21.0 (216) | 100 |
| 1e6 | 1 | None | 50 / 50 | 0 | 0 | 0 | 23.7 (39) | 50 |
| 1e6 | 1 | LongOnly | 89 / 100 | 4 | 7 | 0 | 17.8 (112) | 88 |
| 1e6 | 1 | Capped | 98 / 100 | 0 | 2 | 0 | 18.2 (97) | 99 |
| 1e6 | 1e-4 | None | 50 / 50 | 0 | 0 | 0 | 23.7 (39) | 50 |
| 1e6 | 1e-4 | LongOnly | 91 / 100 | 2 | 7 | 0 | 24.4 (543) | 91 |
| 1e6 | 1e-4 | Capped | 100 / 100 | 0 | 0 | 0 | 18.1 (97) | 100 |
| 1e8 | 1 | None | 43 / 50 | 7 | 0 | 0 | 306.2 (1800) | 43 |
| 1e8 | 1 | LongOnly | 94 / 100 | 3 | 3 | 0 | 29.9 (1140) | 93 |
| 1e8 | 1 | Capped | 98 / 100 | 0 | 2 | 0 | 15.9 (85) | 100 |
| 1e8 | 1e-4 | None | 43 / 50 | 7 | 0 | 0 | 306.2 (1800) | 43 |
| 1e8 | 1e-4 | LongOnly | 95 / 100 | 3 | 2 | 0 | 18.4 (68) | 95 |
| 1e8 | 1e-4 | Capped | 100 / 100 | 0 | 0 | 0 | 15.8 (85) | 99 |

| n | converged | cap reached | no decrease | mean (max) iterations | converged, 0.005 |
|---|---|---|---|---|---|
| 3 | 396 / 400 | 0 | 4 | 7.8 (27) | 398 |
| 5 | 398 / 400 | 2 | 0 | 12.8 (112) | 398 |
| 10 | 387 / 400 | 0 | 13 | 38.5 (994) | 390 |
| 20 | 379 / 400 | 16 | 5 | 45.5 (1319) | 379 |
| 50 | 374 / 400 | 17 | 9 | 53.6 (1800) | 368 |

In total 1934 of 2000 converge at `proximal = 0.05` and 1933 at `0.005`, so 0.05 is kept: on these instances the two
are tied. No run failed to factor `Q` or to solve the QP, and no converged run fails `expect_risk_parity_kkt` at
either value. The other 66 stop at the cap (35) or with no decrease (31), all LongOnly, Capped at `cond(Sigma)` 1e6 and
1e8 (4 runs) or without a box at 1e8; they are left to the Newton-KKT polish that follows the descent. The small
differences between the scales on LongOnly are rounding: the two scales start from minimum-variance points that
differ in their last digits. `RiskParityDescentFuzz.ConvergedRunsAreKkt` runs the same grid at `n` in
{3, 5, 10, 20} and seeds 1 and 2 (640 runs), checks every converged run with `expect_risk_parity_kkt` and
requires the 623 converged runs measured there.

**The identity `w'grad f = 0`** in `expect_risk_parity_kkt` is checked relative to `sum_i abs(w_i)` times the
largest term of `grad f`, not to `sum_i abs(w_i grad f_i)`. At a KKT point where the active bounds are `w_i = 0`,
`grad f_i = budget + beta_i` with `beta_i w_i = 0`, so that sum is `abs(budget) sum_i abs(w_i)`, and the identity
itself makes `budget` vanish to the KKT residual. With the old scale, on the earlier instances, 6 converged runs of the first sweep failed
on it (relative residual `4.8e-6` to `4.5e-4`, all LongOnly with uneven budgets); with the new one none does. A
wrong gradient still fails: with the second of the three terms of the reference gradient scaled by 1.001, 394 of
the 624 runs of `ConvergedRunsAreKkt` that converged when this was measured on the earlier instances (its floor is now
623) fail the identity (the rest sit at `f = 0`, where every term
vanishes and an error proportional to the gap cannot be seen). `RiskParityDescentFuzz.HomogeneityIdentityCatchesAWrongGradient`
keeps that check on one converged LongOnly point with weights at 0 (`n = 10`, `cond(Sigma) = 1e6`). Windows, GCC 16.2 (MSYS2 UCRT64).

### The KKT residual and its noise

`evaluate_kkt` and `residual_of` (`src/riskparity/kkt.cpp`) measure a candidate `(w, lambda, beta)` against
`grad f = sum_k lambda_k grad h_k + beta`, the equalities `h_k = 0` (budget, target return, `vol(w) = s`), the
bounds and the complementarity of `beta`. Every residual is divided by its **noise**: how far it moves when each
weight is rounded once (`abs(dw_j) <= eps abs(w_j)`) and each of its terms is rounded once,

- stationarity of weight `j`: `eps ((abs(H) abs(w))_j + 2 sum_i abs(J_ij) (abs(c_i) + b_i) + sum_k abs(lambda_k
  grad h_k)_j + abs(beta_j))`, `H` the Hessian of the Lagrangian;
- budget: `eps (sum_i abs(w_i) + 1)`; target return: `eps (sum_i abs(mu_i w_i) + abs(r))`; volatility:
  `eps (sum_i abs(w_i p_i) / vol + s)`;
- a bound: `eps max(abs(bound), max abs(w))`.

`is_kkt_point` is the one decision of `OPTIMAL` against `MAX_ITER`: the largest ratio is at most
`kKktNoiseFactor = 16`, and the budget row is resolved, `eps (sum_i abs(w_i) + 1) <= 1e-12` (the feasibility
tolerance of the QP solver, `qp::kRelativeTolerance`). Measurements (Windows, GCC 16.2, MSYS2 UCRT64):

- **The noise is tight** (earlier instances). Moving every weight of an accepted point by one ulp in a random direction (20 trials
  per point; `n = 20`, uneven budgets, `cond(Sigma)` 1e2 to 1e8, the three boxes, with and without a target
  volatility: 24 points) moved the stationarity by 0.51 to 0.84 of its noise. A residual of 1 is what one ulp of
  the weights does.
- **The residual is computed far below its noise** (earlier instances). Against the same residual in `long double` the computed
  stationarity differed by at most `1.4e-4` of its noise (`n` in {10, 20, 50}, `cond(Sigma) = 1e8`, no box, with
  a target volatility: 12 points). `Sigma w` is the compensated product of the section above.
- **The factor 16.** With the augmented Lagrangian stopping at its own tighter factor of 2 (below) and the
  multipliers settled at the final point, the largest residual of the 2770 `OPTIMAL` results of the sweep below
  is 9.30 (a target return alone, no box, `n = 5`, `cond(Sigma) = 1e6`; 8.47 with Capped and a target return at
  `n = 20`, 8.10 at `n = 3` and `cond(Sigma) = 1e8`; 1.95 on the hardest class, no box with a target volatility
  at 1e8, where it was 1.87 on the earlier instances). Before those two changes it was 15.04 on the earlier
  instances, and 15.37 on these, with the stop at 16: the iteration stopped at the first iterate within 16 and
  the polish did not improve it, so 16 was then the factor, not a margin. It is a margin of 1.7 over the largest
  accepted residual of the plan-8 sweep (8.6 over the 1.87 of the earlier instances), which comes from the shapes
  with a target return and no volatility, where the stop of the augmented Lagrangian does not act. The rejected
  points are far away: a target volatility out of reach leaves a gap of at least `4.4e-3 s`, far above its noise.
- **The budget row.** On a line without bounds (target return, no box) `f` is homogeneous of degree 0 and tends
  to a limit along the line, with a gradient that vanishes like `1 / abs(w)`. The descent follows it and stops
  "converged" at a relative stationarity of `1e-11` with weights of size `1e10` to `1e13`: every residual is
  within its noise there, because the noise grows with the weights, and `1'w = 1` holds to `1e-3` at best. In
  the 28 such runs of the sweep (of 240 with a target return and no box) the descent did not end at an accepted
  point; the 7 that the rule alone turns away have `sum_i abs(w_i)` between 3.0e11 and 3.7e13, and the largest of
  an accepted result is 30.8. The two groups are separated by more
  than ten decades, and the line `eps (sum_i abs(w_i) + 1) <= 1e-12`, at a gross exposure of `1e-12 / eps = 4.5e3`,
  falls inside the gap: 146 times above the largest accepted and `6.7e7` times below the smallest rejected. `1e-12`
  is also the tolerance to which the QP solver holds `1'w = 1`, so a point beyond it no longer resolves the row the
  polyhedron was built on. The rule is in the contract of `solve_risk_parity` (`src/riskparity/risk_parity.hpp`): a
  legitimate portfolio with a gross exposure above about 4.5e3 is `MAX_ITER`. Without the rule 7 of those 28 are
  `OPTIMAL`, with `sum_i abs(w_i)` between 3.0e11 and 3.7e13 (uneven budgets at `n` 10 and 20, an equal one at
  `n = 5`); the other 21 are `MAX_ITER` either way. With it all 28 are `MAX_ITER`
  (`RiskParity.TargetReturnWithoutMinimizerIsMaxIter`, `Kkt.WeightsThatDoNotResolveTheBudgetAreNotAKktPoint`).

**The stationarity band is not reached with a target volatility and free weights.** The noise above is in the
units of the residual, not of the gradient. Where the weights are held by bounds, or the point is the Spinu
solution, the noise is far below the band of `cond(Sigma)` and `expect_risk_parity_kkt` passes. With a target
volatility and few or no active bounds the solution moves towards the minimum-variance portfolio, long and
short, and there `abs(H) abs(w)` is larger than the terms of `grad f` by up to `4e13` at `cond(Sigma) = 1e8`
(`2 J'J` grows like `cond(Sigma)^2` while the terms of `grad f = 2 J'r` cancel): one ulp of the weights moves
the stationarity by more than the band. The measured `eps (abs(H) abs(w))_j` over the largest term (earlier instances), `n = 20`,
uneven budget:

| cond(Sigma) | no box, target return | no box, target volatility | no box, both | LongOnly, target volatility | Capped, target volatility | band |
|---|---|---|---|---|---|---|
| 1e6 | 6.8e-12 | 1.6e-7 | 1.5e-6 | 1.2e-10 | 7.7e-16 | 1e-9 |
| 1e8 | 4.2e-12 | 8.4e-3 | 2.4e-3 | 1.0e-8 | 9.7e-16 | 1e-7 |

One ulp of the weights there moves the stationarity by more than the band, so the band alone cannot be the
test of these points. `expect_risk_parity_kkt` passes the stationarity residual `r` in one of two ways, both
computed in the test in `long double` and without `src/riskparity/kkt.cpp`:

- relative to the largest term of `grad f` or of the constraints, at the band of `cond(Sigma)`; or
- entry by entry, `abs(r_j) <= 2 floor_j`, with the **rounding floor**
  `floor_j = eps ((abs(hess f - lambda_vol hess vol) abs(w))_j + sum_k abs(lambda_k grad h_k)_j + abs(beta_j))`:
  how far `r_j` moves when every weight and every multiplier of the result, all doubles, is rounded once. The
  Hessian is the central difference of the reference gradient (step `1e-12` of the largest weight), which
  agrees with the analytic Hessian of the library within `1e-6` at `cond(Sigma)` up to `1e6`
  (`RiskParityReference.HessianByDifferencesMatchesTheObjective`) and within `3.2e-7` of its largest entry at
  `1e8` (three seeds, `n = 6`; the test asks `1e-3`:
  `RiskParityReference.HessianByDifferencesMatchesTheObjectiveAtCondition1e8`). The gap is not the analytic
  Hessian being rounded by `eps cond(Sigma)`, as this section once said: against 80 digits (earlier instances) the analytic one is
  right to `1.3e-15` of its largest entry at every `cond(Sigma)` up to `1e10` at the points of the measurement
  of the Newton multipliers below (eta there), so the gap is that of the central differences. At the random
  weights of this test the two Hessians were not measured apart.

The multipliers in `r` are not always those of the result. Where the free weights (those whose `beta` is 0)
determine the equality multipliers (at least as many free weights as equalities, their gradients over the free
weights independent by the `1e-12` rule of the QP solver), the test refits them at the weights of the result,
by least squares over the free weights in `long double` (modified Gram-Schmidt), with the `beta` of a held
weight what its row then leaves: the multipliers the KKT gate of the library checks
(`least_squares_fit`). Elsewhere (every weight held, dependent rows) `r` takes the multipliers of the result.
The multipliers of the result are then checked against the KKT point refined from it in `long double` (the six
Newton steps of the forward error below): `max_k abs(lambda_k - lambda_ref_k) / max_k abs(lambda_ref_k)` may
exceed that of the fit by at most the band (`expect_multipliers_toward_the_reference`; not checked where the
refinement meets a singular system or `long double` is `double`). The reason: the polish returns multipliers
corrected by a Newton solve (the section on the polish), and with the exact multipliers the stationarity at
weights off the KKT point is `H dw`, which the least-squares fit absorbs in part. On every result measured here
the corrected multipliers stay within the floor anyway (at most 0.27 of it in the measurement of the Newton
multipliers below), so the refit decides nothing there; it is a guard. The guard cuts both ways
(`RiskParityReference.KktUsesLeastSquaresMultipliersForTheFloor`): at `s = sigma_min (1 + 1e-10)` (`n = 5`,
`cond(Sigma) = 1e4`) the fit is `5.1e-8` from the refined multipliers, and multipliers half as far from them,
all of it in the budget multiplier, leave a stationarity at `4.4e6` times its floor, yet pass; at
`s = 1.2 sigma_min` of another shape, multipliers `1e-6` from the refined ones fail, by the comparison with the
reference alone.

**Portability of the reference.** Every sum of the reference that cancels (`Sigma w`, `w' Sigma w`, `c'(c - b)`,
the residual of the Spinu reference in `spinu_test.cpp`) is a compensated dot product (Dot2: `fma` and TwoSum),
templated on the type, so it is as accurate as a sum in twice the working precision whether `long double` is
the 64-bit-mantissa type of GCC or `double` (MSVC). In `double` the entries of `Sigma w` at `cond(Sigma) = 1e8`
agree with those in `long double` within one ulp
(`RiskParityReference.CompensatedDotInDoubleIsAccurateAtCondition1e8`). Two checks need more than `double` and
are skipped, with a message, where `long double` is `double` (`kLongDoubleIsExtended`): the rounding-floor
branch of the stationarity (the floor is `eps` of the magnitudes of the terms, about what summing them in
`double` adds), and the refinement of a point beyond `double` for the forward error. The tests run under GCC (MSYS2 UCRT64) and under MSVC, where the two
checks above are skipped.

The equalities, the bounds, the sign and complementarity of the multipliers and the identity `w'grad f = 0`
are checked at the band as before. The factor 2 is not the solver's 16: the largest `abs(r_j) / floor_j` of the
218 results of the sweep below that are above the band is 0.29 (1.36 over the 251 of the earlier instances; there one
ulp of the weights moved the residual by up to 0.84 of the solver's noise, which is a little larger than this floor). The floor is not a wider band:
the weights of such a result moved by 64 ulps each are more than 8 floors away and are rejected, while their
equalities still hold at the band (`RiskParityReference.KktRejectsAStationarityAboveItsRoundingFloor`).

**Forward error of the results above the band.** A stationarity at its floor bounds the backward error, not
the distance to the solution, which is larger by the conditioning of the KKT system. So every result that
passes by the floor alone is also compared with the KKT point of its own active set refined in `long double`
(`tests/unit/fuzz/risk_parity/risk_parity_refinement.hpp`: six Newton steps on the KKT system from the result,
with the reference gradient and its Hessian by differences; the refined point is at most `1.7e-4` of the floor
from stationarity and satisfies the equalities at the band). Of the 2770 `OPTIMAL` results of the sweep below,
218 are above the band, all with a target volatility, none at `cond(Sigma) = 1e2`:

| cond(Sigma) | OPTIMAL | above the band | largest `abs(r_j) / floor_j` | largest weight error | band | largest value error over `band (1 + f)` |
|---|---|---|---|---|---|---|
| 1e2 | 686 | 0 | | | 1e-10 | |
| 1e4 | 694 | 14 | 0.19 | 1.1e-13 | 1e-10 | 1.3e-7 |
| 1e6 | 698 | 98 | 0.29 | 1.1e-12 | 1e-9 | 1.0e-6 |
| 1e8 | 692 | 106 | 0.27 | 7.0e-6 | 1e-7 | 7.4e-6 |

The weight error is `max_i abs(w_i - w_ref_i) / max_i abs(w_i)` and the value error `abs(f(w) - f(w_ref))`.
By `cond(Sigma)`, `n` and box, the largest weight error (number of results above the band):

| cond(Sigma) | box | n = 3 | n = 5 | n = 10 | n = 20 | n = 50 |
|---|---|---|---|---|---|---|
| 1e4 | LongOnly | | 9.2e-17 (2) | | | |
| 1e4 | None | | | 4.7e-16 (4) | 1.5e-14 (6) | 1.1e-13 (2) |
| 1e6 | LongOnly | 1.2e-16 (2) | 5.3e-15 (6) | 1.3e-15 (4) | | |
| 1e6 | None | 9.0e-17 (4) | 1.4e-13 (12) | 3.7e-13 (24) | 1.1e-12 (22) | 6.1e-13 (24) |
| 1e8 | LongOnly | 8.7e-17 (3) | 2.1e-14 (4) | 5.9e-14 (5) | | 2.7e-13 (2) |
| 1e8 | None | 9.9e-17 (8) | 2.1e-14 (12) | 5.9e-12 (24) | 3.0e-12 (24) | 7.0e-6 (24) |

On the seeds of this sweep the weights are right at the band, by a factor of at least `1e3` except at
`cond(Sigma) = 1e8`, `n = 50`, no box, where the largest weight error, 7.0e-6, is 70 times above the band of 1e-7
(`n <= 20` stays under 6e-12 there); the 24 results of that class are the conditioning limit of the next paragraph.
That is not true of the class in general.
`RiskParitySolveSweep.OptimalRunsAreKkt` requires, on its grid of 576 shapes, that every `OPTIMAL` passes
`expect_risk_parity_kkt`, that at least 540 of the 576 are `OPTIMAL` (the 36 others are `MAX_ITER` without a point: 8 at
`n = 3` with a target return and no box, 12 at `n = 3` long only with both targets and 16 at `n = 5` capped with both
targets), and that every one above the band has a weight error at most the band and a value
error at most `band (1 + f)` against a refined point within `1e-2` of the floor. The refinement needs a `long
double` wider than `double`; where it is not (MSVC) that part is skipped.

**The weights at `cond(Sigma) = 1e8` with a target volatility are a conditioning limit.** With a target
volatility and few or no active bounds at `cond(Sigma)` around `1e8`, the weights are determined only to the
rounding of the KKT residual amplified by the conditioning of the reduced problem. A result at the rounding
floor can differ from the exact KKT point of its active set by up to about `1e-4` relative (the measured maximum
is `2.8e-5`, below), while its objective value is right to its band. Measured over the class (no box, target volatility, uneven budget, `cond(Sigma) =
1e8`; `RiskParityFuzz.DISABLED_MeasureWeightErrorDistribution`, one run):

| n | OPTIMAL | weight error: median | 90% | 99% | max | above `1e-7` | value error over `band (1 + f)`, max |
|---|---|---|---|---|---|---|---|
| 20 | 199 / 200 | 7.2e-13 | 4.7e-12 | 6.8e-7 | 1.3e-5 | 3 | 2.2e-5 |
| 50 | 30 / 30 | 2.0e-12 | 6.6e-10 | 2.9e-6 | 2.8e-5 | 2 | 7.8e-5 |

(All 30 seeds at `n = 50` ran within the 240 s cap; seed 27 of `n = 20` is `MAX_ITER` with every number NaN.) The worst
of the class at `n = 20` is seed 135: weight error `1.3e-5`, a refined point at `2.7e-5` of the floor and
`abs(f(w) - f(w_ref)) = 2.2e-12`; at `n = 50` it is seed 18, with `2.8e-5` (a refined point at `6.0e-5` of the floor). On the earlier instances seed 31696
(`n = 20`) was the one first seen, with a weight error of `2.3e-5` (`5.2e-5` with the augmented Lagrangian stopping
at 16), a stationarity of 0.22 of its rounding floor, a refined point at `1.3e-5` of the floor (both are KKT points of
the same active set to their rounding) and `abs(f(w) - f(w_ref)) = 1.6e-12` at `f = 2.6e-5`; its instance is another
since plan 8 (weight error `6.2e-12`).

A posteriori bounds were tried and do not predict the error: with `K` the KKT matrix of the active set,
`max_j sum_l abs(K^-1)_jl floor_l` (any point at the floor) gives 0.61 at seed 135 and a median of 0.015 over
the 199 results of `n = 20` (largest 1.7), and `max_j sum_l abs(K^-1)_jl abs(r_l)` (the residual `r` of the result itself)
gives 0.0068 there; the error is at most `1.9e-3` of the second over the 229 results. Both hold, neither is sharp (the
componentwise absolute values lose the cancellation that makes the true error small). The contract of
`solve_risk_parity` states the limitation. `RiskParityFuzz.WeightsAtCondition1e8BelowTheirMeasuredCeiling`
(seeds 1 to 4 of `n = 20`, with weight errors from 9.8e-14 to 6.2e-12, 0.8 s) checks, against the refined point, the value error at
`band (1 + f)` and, as the guard of the weights, the measured ceiling of `1e-4` (the largest of the 229 is
`2.8e-5`); the band of the weights is not applied to a result at its floor. It also checks the weight error
within 4 times the residual bound, a sanity bound only, since that bound is at least 5e2 times loose.

**The multipliers of the target-volatility path (measured, not derived).** The weights and the value above are
what the floor bounds; the multipliers are not. With a target volatility the least-squares multipliers at the
rounded weights lose the correction `H dw` that the Newton step would have put in them, and that error was the
one reported until the polish kept the Newton correction of the multipliers at its last point (the section on
the polish; the measurement below). Against the multipliers of the KKT point refined in `long double` from the
result (`max_k abs(lambda_k - lambda_ref_k) / max_k abs(lambda_ref_k)`; no box, equal budget, `n` in
{5, 8, 12, 20}, 4 seeds, `s` between `sigma_min` and `1.5 sigma_min`;
`RiskParityFuzz.TargetVolatilityMultipliersByCondition`, 0.5 s; one run, Windows, GCC 16.2, MSYS2 UCRT64; the
largest are at `n = 8`, seeds 405, 404 and 404; the least-squares column is that of the 8 shapes with `n` in {5, 8} of
the measurement below):

| cond(Sigma) | least squares (before) | largest multiplier error | band | ceiling pinned by the test |
|---|---|---|---|---|
| 1e4 | 3.7e-10 | 8.1e-13 | 1e-10 | 8.2e-12 |
| 1e6 | 4.2e-4 | 2.4e-12 | 1e-9 | 2.4e-11 |
| 1e8 | 2.2e-2 | 9.3e-11 | 1e-7 | 9.4e-10 |

With a box (`RiskParityFuzz.TargetVolatilityBoundMultipliersByCondition`, 0.3 s: LongOnly and Capped, equal and
uneven budgets, `n` in {5, 8}, 4 seeds, 32 results per `cond(Sigma)` with 60 to 68 weights held), the bound
multipliers of the held weights are compared by value with those of the refined point, `beta_ref_j` the entry of
`grad f - A'lambda` there, over `max(abs(lambda_budget*), abs(lambda_vol*) s, max abs(beta_ref))`, and their
sign and complementarity are checked by `expect_risk_parity_kkt`:

| cond(Sigma) | equalities | ceiling | held `beta` | ceiling | held `beta` of the least-squares fit |
|---|---|---|---|---|---|
| 1e4 | 1.4e-14 | 1.9e-13 | 8.6e-15 | 1.8e-13 | 5.7e-12 |
| 1e6 | 2.1e-13 | 2.6e-11 | 2.6e-12 | 3.6e-12 | 1.1e-8 |
| 1e8 | 2.0e-12 | 1.4e-9 | 1.6e-10 | 1.5e-9 | 3.1e-7 |

The last column is the `beta` of the least-squares fit at the same weights (`least_squares_fit`, without the
`H dw` term of `beta_j = (grad f - A'lambda_new)_j + (H dw)_j`): the test requires it more than 10 times above
the ceiling, so a `beta` without that term fails. Capped results hold most weights and leave few free ones: there
the least-squares and the corrected multipliers agree to `1e-15`; the difference is in the LongOnly shapes, which
hold one to four weights at 0.

The ceilings of the second table are those of the earlier instances (the measured largest times 10: measured, not
derived), which the new measurements, 14 to 700 times below them for the equalities and 1.4 to 21 times for `beta`
(9.4 at 1e8, under the factor of 10), still pass; the ceilings of the first table are the new largest times 10.
**The targets of plan 5 are missed.** The plan (and the research note) asked for `1e-12` up to `cond(Sigma) = 1e6` and
`1e-9` at `1e8` with a slack of 10 over the measured. At `1e4` the measured `8.1e-13` meets `1e-12`, not with the slack
(the ceiling `8.2e-12`). At `1e6` the measured `2.4e-12` (`n = 8`, seed 404; the other seven shapes are at most
`1.2e-12`) misses `1e-12` by a factor of 2.4, before any slack. At `1e8` the measured `9.3e-11` (`2.6e-10` in either
multiplier on the probe cases of the measurement below) is under `1e-9`, and so is the ceiling `9.4e-10`; on the
earlier instances, `5.3e-10` and `5.3e-9`, it was not. The measurement, not that model, sets the ceilings. The weights
of the same results are right to `6e-13` or better.

**Why the rank test of the QR does not veto the correction** (earlier instances; not remeasured). At `cond(Sigma) = 1e8` with `n >= 12` and no box the
scaled KKT matrix `[H, -A'; A, 0]` fails `has_full_column_rank` (`16 N eps` times its largest column norm) on about
half the shapes, yet the solution of that same QR is the right correction: with the test as a veto, `n = 12` of
`TargetVolatilityMultipliersByCondition` (seed 404) kept least-squares multipliers wrong by `4.1e-2` and `n = 20`
(seeds 403 to 405) by up to `1.3e-3`; corrected they are within `2.9e-11`. A step of the weights still requires the
rank test (`newton_step`); the correction requires a finite solution, and contraction, stationarity and sides
decide. A wider grid (`n` in {5, 8, 12, 20, 50}, no box and long only, both budgets, seeds 601 to 604, two at
`n = 50`; `po::solve` with the veto and without, the KKT point refined in `long double`; one run) shows the change
at `1e8` and its effect at `1e10`:

| cond(Sigma) | results `OPTIMAL` | refused by the rank test before | largest error, before | largest error, now | largest error, least squares |
|---|---|---|---|---|---|
| 1e4 | 72 / 72 | 0 | 2.7e-13 | 2.7e-13 | 1.5e-8 |
| 1e6 | 72 / 72 | 0 | 4.2e-12 | 4.2e-12 | 1.6e-5 |
| 1e8 | 72 / 72 | 12 (3 of 8 at `n = 12`, 5 of 8 at `n = 20`, 4 of 4 at `n = 50`, all no box) | 1.3 | 3.8e-10 | 1.3 |
| 1e10 | 60 / 72 | 17 (no box at `n >= 8`, one long only at `n = 20`) | 31 | 0.25 | 48 |

On the earlier instances, at `1e8`, of the 16 thinnest-class shapes (`n` in {12, 20}, uneven budget, 8 seeds) two fail the rank test: seed 7
of `n = 20` goes from `8.8e-4` to `2.5e-12`, and seed 4 from `1.3e-2` to `8.7e-5`, its weights being `2.3e-5` from
the refined point (the multipliers can only be as right as the weights). Every one of the 35 corrections the rank
test refused before was accepted without it except two at `1e10` that the contraction test refuses
(ratio 0.57 and 0.69), and every one accepted is nearer the refined KKT point than least squares. The weights and
the statuses of all 550 results of these measurements (527 `OPTIMAL`; the near-`sigma_min` cases below included)
are bit for bit the same
with and without the veto. At `1e10` the error left is that of the weights (`1.5e-2` at the worst result, `n = 20`).

`RiskParity.HomogeneityIdentity` (`budget = -target_return r - volatility s - sum_i beta_i w_i`) holds with the
corrected multipliers at `cond(Sigma)` `1e4` and `1e6` (relative gap `4.8e-17` and `4.4e-17`; box `[0, 0.3]`, both
targets), within `4.8e-16` of the scale: the largest measured times 10. On the earlier instances the least-squares
multipliers at the same weights left gaps of the same order as the corrected ones (`1.1e-15` and `1.0e-16`): the
identity checks that the corrected `beta` are consistent to rounding, not that they are nearer the KKT point. Above
`cond(Sigma) = 1e8` an `OPTIMAL` result is not validated: at `1e10` the correction reaches `2.1e-7` in either
multiplier on the 8 shapes of the measurement below and, on the earlier instances, `0.25` on the wider grid above,
where the weights are themselves off. `Data` accepts `cond(Sigma) >= 1e10`. The path with a target return alone is
not made worse by the correction: `RiskParity.TargetReturnOnlyCorrectionNeverWorsensTheFit` (`cond(Sigma) = 1e8`, `n = 8`,
LongOnly and Capped, both budgets, 4 seeds) finds the corrected multipliers within `3.1e-14` of the refined ones where
the least-squares fit is off by up to `1.3e-11`.
`RiskParity.MultipliersAreDerivatives` runs on well-conditioned data only.

**The QP path, for comparison (measured once).** `RiskParityMultiplierMeasure.DISABLED_MeasureQpMultipliers`
(`tests/unit/fuzz/risk_parity/risk_parity_multiplier_measure_test.cpp`, run with
`--gtest_also_run_disabled_tests`) solves `MIN_VARIANCE` with a target return half way between the return of the
minimum-variance portfolio and the largest mean, no box and long only, `n` in {5, 8, 20}, seeds 1 to 4, and
compares the multipliers of Goldfarb-Idnani with the KKT point of the result's active set solved in `long
double` (`2 Sigma w_F = lambda_budget 1 + lambda_return mu_F`, `1'w = 1`, `mu'w = r`, four Newton steps from
the result, the last changing the multipliers by at most `5.9e-18`; `beta` of a held weight is what its row
leaves, over `max(abs(lambda_budget*), abs(lambda_return*) max abs(mu), max abs(beta*))`). Item 5 of the
research note asked for it. Largest errors, one run:

| cond(Sigma) | no box: weights | multipliers | long only: weights | multipliers | held `beta` |
|---|---|---|---|---|---|
| 1e4 | 4.2e-14 | 6.4e-13 | 2.5e-15 | 7.9e-14 | 2.8e-14 |
| 1e6 | 1.4e-12 | 1.4e-11 | 3.5e-14 | 8.9e-15 | 1.6e-14 |
| 1e8 | 1.9e-10 | 1.2e-9 | 3.5e-14 | 7.2e-13 | 6.0e-13 |

The QP multipliers come out of the same factorization as the weights and are of the same order as the corrected
ones of `RISK_PARITY` (`1.2e-9` and `9.3e-11` at `1e8`), below `eps cond(Sigma)`, as the research note expected.
Nothing in the QP path changes; no test pins these numbers.

**Future work: phase 2 (a plan, not implemented).** The correction above is one Newton solve in `double`: what
limits it is the KKT residual evaluated in `double`, amplified by the KKT system (eta, in the measurement below,
shows the Hessian is not the limit). A second phase would evaluate that residual in double-double, with the weights carried as a pair
`w_ext = (w_hi, w_lo)` inside the polish (`Sigma w` and the compensated dot products of `objective.cpp` and
`kkt.cpp` over the pair), take two or three steps reusing the factorization of the last one, stop by the rule of
Demmel et al. (a correction at most `eps` of the multipliers, or a ratio of successive corrections at least
1/2), and still report `w_hi`; no `long double`, so that MSVC keeps it. It would be activated if a use needs
multipliers validated at `cond(Sigma)` around `1e10` (now up to `2.1e-7` on the 8 shapes of the measurement below
and, on the earlier instances, `0.25` on the wider grid above, there with weights `1.5e-2` off), or if a measurement at `1e8` exceeds the ceilings above (a regression of the one solve, or a class of
problems not covered by these cases); until then the ceilings above are the documented accuracy and `1e10` stays
not validated.

**The Newton multipliers at the returned point (measured, not derived).**
`RiskParityMultiplierMeasure.DISABLED_MeasureNewtonMultipliers`
(`tests/unit/fuzz/risk_parity/risk_parity_multiplier_measure_test.cpp`, 0.9 s, run with
`--gtest_also_run_disabled_tests`) takes the point `po::solve` returns, builds on it the Newton system of
`newton_polish` from the current `evaluate_kkt` (rows of `A` scaled to the largest entry of `H`, one Householder
QR) and keeps `lambda_LS + dlambda` with `w` unchanged, `lambda_LS` being the multipliers of the result. Both are
compared with the KKT point refined in `long double` from the result: `lambda_vol` relative to `abs(lambda_vol*)`,
`lambda_budget` relative to `max(abs(lambda_budget*), abs(lambda_vol*) s)`. The cases are those of the 80-digit
probe of `docs/research/2026-11-16-multiplier-accuracy.md` (no box, equal budget, `s = 1.2 sigma_min`) on the
covariances of `random_covariance` (not the probe's), and the 8 shapes with `n` in {5, 8} of
`RiskParityFuzz.TargetVolatilityMultipliersByCondition` (seeds 402 to 405), each at
`cond(Sigma)` in {1e4, 1e6, 1e8, 1e10}. One run (Windows, GCC 16.2, MSYS2 UCRT64); the probe cases:

| cond(Sigma) | n, seed | weight error | `lambda_LS` vol / budget | `lambda_LS + dlambda` vol / budget | `abs(dlambda) / abs(lambda_LS)` | stationarity over floor, `lambda_LS` / Newton |
|---|---|---|---|---|---|---|
| 1e4 | 5, 7 | 1.6e-16 | 9.0e-11 / 8.0e-11 | 1.1e-14 / 2.4e-14 | 8.8e-11 | 0.16 / 0.17 |
| 1e4 | 5, 11 | 1.7e-16 | 2.1e-11 / 5.2e-11 | 4.1e-14 / 5.1e-15 | 2.2e-11 | 0.17 / 0.08 |
| 1e4 | 8, 3 | 1.3e-16 | 1.9e-13 / 5.6e-11 | 1.7e-15 / 6.6e-15 | 9.5e-12 | 0.67 / 0.28 |
| 1e6 | 5, 7 | 1.4e-16 | 2.0e-7 / 3.9e-7 | 4.9e-13 / 1.3e-13 | 2.6e-7 | 0.18 / 0.22 |
| 1e6 | 5, 11 | 7.1e-16 | 1.7e-8 / 8.2e-8 | 1.3e-12 / 9.6e-13 | 2.1e-8 | 0.45 / 0.06 |
| 1e6 | 8, 3 | 5.5e-16 | 8.2e-11 / 6.9e-8 | 1.3e-13 / 1.4e-13 | 1.2e-8 | 0.14 / 0.05 |
| 1e8 | 5, 7 | 8.4e-17 | 1.7e-4 / 7.5e-4 | 1.6e-11 / 8.6e-11 | 4.2e-4 | 0.07 / 0.10 |
| 1e8 | 5, 11 | 5.8e-16 | 5.4e-5 / 5.4e-4 | 7.6e-11 / 2.6e-10 | 9.1e-5 | 0.60 / 0.09 |
| 1e8 | 8, 3 | 2.4e-15 | 4.0e-6 / 6.3e-4 | 1.3e-11 / 4.0e-11 | 1.1e-4 | 0.34 / 0.08 |
| 1e10 | 5, 7 | 2.0e-16 | 7.0e-2 / 0.64 | 1.4e-10 / 1.8e-9 | 0.41 | 0.01 / 0.02 |
| 1e10 | 5, 11 | 2.7e-15 | 1.9e-2 / 0.50 | 1.1e-9 / 7.4e-9 | 7.0e-2 | 0.15 / 0.04 |
| 1e10 | 8, 3 | 5.3e-14 | 19 / 40 | 3.5e-8 / 2.8e-8 | 1.1 | 0.45 / 0.12 |

The 8 shapes of the fuzz test with `n` in {5, 8}, largest by `cond(Sigma)` (combined is the error of that test,
`max_k abs(lambda_k - lambda_ref_k) / max_k abs(lambda_ref_k)`):

| cond(Sigma) | `lambda_LS`, combined | `lambda_LS + dlambda`, combined | `lambda_LS + dlambda`, either multiplier | `abs(dlambda) / abs(lambda_LS)` | stationarity over floor, Newton |
|---|---|---|---|---|---|
| 1e4 | 3.7e-10 | 8.1e-13 | 8.2e-13 | 3.6e-10 | 0.15 |
| 1e6 | 4.2e-4 | 2.4e-12 | 3.2e-12 | 4.4e-4 | 0.24 |
| 1e8 | 2.2e-2 | 9.3e-11 | 1.7e-10 | 2.7e-2 | 0.20 |
| 1e10 | 1.5e2 | 5.8e-8 | 2.1e-7 | 1.3 | 0.18 |

The `lambda_LS` column reproduces the table above (3.7e-10, 4.2e-4, 2.2e-2). At `cond(Sigma) = 1e10` the correction,
which does not require the rank test of the QR (`kRankFactor`), is accepted on all eight shapes and on the three probe
cases (the largest, `2.1e-7` in either multiplier, at `n = 5`, seed 405, whose weights are `3.1e-16` from the refined
point). On the earlier instances the scaled KKT matrix of the four shapes with `n = 8` failed that test, as
`newton_step` would find it, and the largest correction was `2.7e-5`; before the veto went, they gave no correction
(`n = 5` alone: `6.3e-7` combined, `2.3e-6` in either multiplier).

- **The Newton correction recovers the multipliers.** One solve at the returned point, all in `double`, brings
  them to at most `8.2e-13` at 1e4, `3.2e-12` at 1e6, `2.6e-10` at 1e8 and `2.1e-7` at 1e10 (in either multiplier,
  over the 8 shapes and the 3 probe cases of each), from `3.7e-10`, `4.2e-4`, `2.2e-2` and `1.5e2` (combined; at 1e8
  by 3.7 to 8.9 orders of magnitude in every case), with the weights untouched. That is inside the research note's
  prediction at 1e8 (its model put it at `1e-10` to `5e-10` on its own cases).
- **The 80-digit check** (earlier instances; the script is not kept). An independent refinement of the dumped cases (`PORTOPT_MULTIPLIER_DUMP=path` writes
  each problem, the point and the analytic Hessian in `%.17g`; a script outside the repository refines
  `(w, lambda_LS)` by 12 Newton steps in 80-digit decimal arithmetic, Jacobian by forward differences of step
  `1e-35`, the method of section 0 of the research note) agrees with the `long double` reference within
  `1.2e-16` in either multiplier in all 44 cases, `1e10` included, so the errors above are also those against
  80 digits.
- **eta, the error of the analytic Hessian** (earlier instances; the `eta` column of their tables is not repeated). The Lagrangian Hessian `evaluate_kkt` returns at
  `(w, lambda_LS)` differs from the 80-digit one at the same point by at most `1.3e-15` of its largest entry
  (`1.6e-15` entry by entry) at every `cond(Sigma)` up to 1e10: it is right to its rounding at these points,
  not off by `eps cond(Sigma)` as the research note modelled it. What limits `lambda_LS + dlambda` is the
  residual, evaluated in `double`, amplified by the KKT system, not `H`. (The `3.2e-7` of the difference check
  above, at random weights, was not split between the two Hessians.)
- **The stationarity with the new multipliers.** At `(w, lambda_LS + dlambda)` the stationarity is at most 0.28
  of its rounding floor in every case (with `lambda_LS` it reaches 0.70 at 1e4, 0.91 at 1e6, 1.24 at 1e8 and 0.94 at
  1e10), and the same in units of the noise of `kkt.cpp`: both the gate of `is_kkt_point` and that of
  `expect_risk_parity_kkt` would pass with them.
- **For a safeguard.** `abs(dlambda) / abs(lambda_LS)` is the error of `lambda_LS`: up to `2.7e-2` at 1e8 (one
  of the 11 cases above `1e-2`) and `4.4e-5` to 1.3 at 1e10 (10 of the 11 above `1e-2`). A second solve from
  `(w, lambda_LS + dlambda)` changes the multipliers by at most `5.6e-4` of the first correction (at 1e4, where the
  correction is itself at the rounding), `6.8e-6` at 1e6, `7.7e-8` at 1e8 and `5.2e-9` at 1e10.

### The Newton polish

`newton_polish` (`src/riskparity/newton_polish.cpp`) takes the point and the active set of the last QP of the
descent (the signs of `beta`) and iterates Newton on the KKT system over the free weights,
`[H, -A'; A, 0] [dw; dlambda] = -[grad f - A'lambda; h]`, `H` the exact Hessian of the Lagrangian
(`hess f - lambda_vol hess vol`), solved by one Householder QR of the matrix (rows of `A` scaled to the size of
`H`), never an inverse. Decisions:

- **Acceptance.** A step is taken only when it decreases the largest KKT residual in noise units, feasibility
  included; the first step that does not ends the polish, and so does a residual at or below 1 (`improved`
  says whether the residual ended below that of the start). From a descent stopped at a relative residual of `1e-6` (`n = 5`, box `[0, 0.24]`,
  three weights at the cap, with and without a target return) one or two steps reach a residual below 1 and a
  second call takes no step (`NewtonPolish.ReachesRoundingFloor`). At most 20 steps; the most used is 6 at the
  first effort of the solver on the grid below and, on the earlier instances, 10 at `s = sigma_min (1 + 1e-10)`.
- **Multipliers by least squares at every point.** Only `dw` is taken from the Newton system; the multipliers
  of the equalities are then the least-squares solution of the stationarity over the free weights at the new
  weights (Householder QR of the gradients), and those of the held weights are what their rows leave. With the
  multipliers of the Newton system itself, a target volatility near the minimum was not solved (earlier instances): at
  `s = sigma_min (1 + 1e-9)` (`n = 5`, `cond(Sigma) = 1e4`), where `lambda_vol s = -2.5e5`, the polish refused
  its first step from every point the augmented Lagrangian reached and the result was `MAX_ITER`; with the
  estimate it takes 4 to 8 steps to a residual of 0.3.
- **The multipliers are corrected once more at the end** (`corrected_multipliers`): at the point the polish
  ends at, one more Newton solve at the same weights gives `dlambda`, added to the multipliers of the kept
  equalities, and `beta_j = (grad f - A'lambda_new)_j + (H dw)_j` for a held weight. The weights, the status,
  the steps and the KKT gate (`is_kkt_point`, with the least-squares multipliers) do not change; the solver
  makes the correction only at a point that passed the gate, and `po::solve` returns the corrected
  multipliers. The correction is kept only if both solves give a finite solution (not the rank test of the QR,
  which a step of the weights requires: see "Why the rank test of the QR does not veto the correction" above), a
  second solve from `(w, lambda_new)` changes the multipliers by at most half the first correction, every held
  `beta` keeps the sign of its bound, and the stationarity at `(w, lambda_new)` is within 16 times its noise
  (the factor of the gate). The last rule is measured, not derived: near `sigma_min` the contraction alone is
  not enough. On the earlier instances, at `s = sigma_min (1 + 1e-10)` the correction was `6.1e-7` of `abs(lambda)` and
  its second solve contracted by `6.3e-7`, yet it took the multipliers from `5.3e-8` to `5.5e-7` of the refined
  ones and left a stationarity of 958 times its noise, against 0.027 with least squares. A cap on
  `abs(dlambda) / abs(lambda)` cannot tell that correction from good ones: the good corrections at
  `cond(Sigma) = 1e8` reach `2.7e-2`. Away from `sigma_min` every correction measured left at most 0.38 of the
  noise on the earlier instances (280 results at `cond(Sigma) <= 1e8`: the 216 of the grid above, 16 of the thinnest class and the 48
  shapes of the fuzz test, `n` up to 20) and 0.28 on the 44 cases of the measurement above; near `sigma_min` the
  rule refuses the worst corrections but does not separate good from bad (see "The multipliers near the minimum
  volatility" below). On the 18 results of the table below (9 values of `delta` from `1e-2` to
  `1e-10`, no box and long only) the safeguard kept the least-squares multipliers 6 times: at `delta = 1e-10` in both
  boxes, at `1e-8` and `1e-9` without a box and at `1e-6` and `1e-7` long only
  (`RiskParity.NearMinimumVolatilityKeepsLeastSquaresMultipliers`, which requires exactly that; twice, at
  `delta = 1e-10`, on the earlier instances); on the 44 cases of the measurement above it refuses none (on the earlier
  instances it refused the four with `cond(Sigma) = 1e10` and `n = 8` while the rank test was one of its rules). A point far from its KKT
  point (`s = 1.1 sigma_min` solved, `s = sigma_min (1 + 1e-4)` asked) has a contraction of 0.72 and is
  refused (`NewtonPolish.RejectsACorrectionThatDoesNotContract`), and so is a correction that gives a held
  weight the sign of the other side of its bound (`NewtonPolish.RejectsACorrectionThatChangesTheSideOfABound`:
  a weight held at 0 named at an upper bound). With a target return alone (`cond(Sigma) =
  1e8`, `n = 8`, LongOnly and Capped, both budgets, 4 seeds) the corrected multipliers stay within `3.1e-14` of the
  refined ones, never farther than the least-squares fit, which is off by up to `1.3e-11`
  (`RiskParity.TargetReturnOnlyCorrectionNeverWorsensTheFit`).
- **The settling is done at the start too**, whether or not a step follows: the weights the start names at a
  bound are put exactly there and the multipliers are the settled ones, so a point that no step improves does
  not keep the multipliers of the descent (`NewtonPolish.SettlesTheMultipliersWithoutAStep`: a budget
  multiplier off by half the noise of the stationarity left 0.48 of that noise along the budget row before
  this, and at most 0.0502 after, what the least-squares fit at those weights leaves).
- **A multiplier of the wrong sign releases its bound**, and a free weight that a step takes out of the box is
  held at the bound it crossed (`NewtonPolish.ReleasesABoundWhoseMultiplierHasTheWrongSign`).
- **Dependent equalities** (equal means with a target return) are dropped from the system by the rule of the
  QP solver (the part of the gradient outside the span of the earlier ones at most `1e-12` of it); their
  multipliers stay as the descent left them.

Approach to the minimum volatility, `s = sigma_min (1 + delta)` (`sample_data(5, 13)`: `n = 5`, `cond(Sigma) = 1e4`, equal
budget; every result is `OPTIMAL` and passes `expect_risk_parity_kkt` at `1e-10`):

| delta | no box: iterations | lambda_vol s | long only: iterations | lambda_vol s |
|---|---|---|---|---|
| 1e-2 | 353 | -66.8 | 53 | -1.66 |
| 1e-4 | 1196 | -5.2e3 | 190 | -47.2 |
| 1e-6 | 2264 | -6.4e4 | 278 | -490 |
| 1e-8 | 2563 | -6.5e5 | 305 | -4.9e3 |
| 1e-10 | 2692 | -6.5e6 | 320 | -4.9e4 |
| 3e-12 | 2577 | -3.8e7 | | |

`lambda_vol` grows like `1 / sqrt(delta)`, as the square-root shape of the optimal `f` near `sigma_min`
predicts. At `abs(s - sigma_min) <= 1e-12 s` (`problems::kMinimumVolatilityTolerance`, the feasibility tolerance
of the QP that computes `sigma_min`) no iteration is run: the result is the minimum-variance portfolio with
`lambda_vol = -inf` and the other multipliers NaN.

**The multipliers near the minimum volatility (measured, not derived).** The accuracy of the corrected multipliers
stated above holds for `s` at least `1e-3` relative above `sigma_min`. Nearer, it degrades roughly like `1 / delta`,
corrected or not, because the gradients of the budget and of the volatility turn parallel (the LICQ nearly fails).
Measured with `po::solve` against the KKT point refined in `long double` (no box, equal budget, `n` in {5, 8},
seeds 7, 11 and 13, `s = sigma_min (1 + delta)` with `sigma_min` from `frontier_point`; largest over the `OPTIMAL`
results, one run; of the 180 solves, 168 are `OPTIMAL`, 10 `MAX_ITER` (all at `delta <= 1e-7`)
and 2 `INFEASIBLE`, at `1e8` and `1e-10`, where the reference `sigma_min` is above the solver's, as before the
corrections of plan 5):

| delta | cond 1e4: returned | least squares | cond 1e6: returned | least squares | cond 1e8: returned | least squares |
|---|---|---|---|---|---|---|
| 1e-1 | 1.3e-14 | 6.5e-11 | 8.8e-13 | 2.0e-7 | 3.2e-10 | 3.6e-3 |
| 1e-2 | 2.3e-14 | 3.1e-11 | 2.1e-12 | 4.6e-7 | 8.6e-11 | 1.3e-3 |
| 1e-3 | 8.3e-14 | 6.0e-12 | 4.8e-13 | 1.6e-8 | 8.3e-12 | 2.7e-3 |
| 1e-4 | 3.2e-12 | 3.1e-12 | 1.1e-12 | 5.4e-9 | 9.7e-12 | 7.7e-6 |
| 1e-5 | 1.2e-11 | 6.1e-12 | 1.1e-11 | 2.2e-9 | 6.1e-12 | 4.9e-6 |
| 1e-6 | 2.3e-10 | 2.3e-10 | 1.6e-10 | 1.9e-10 | 1.3e-10 | 2.1e-6 |
| 1e-7 | 7.6e-10 | 7.6e-10 | 8.8e-10 | 1.1e-9 | 1.9e-9 | 5.4e-7 |
| 1e-8 | 4.2e-9 | 4.2e-9 | 7.8e-9 | 9.0e-9 | 7.1e-9 | 2.1e-8 |
| 1e-9 | 9.6e-8 | 9.6e-8 | 1.6e-7 | 1.6e-7 | 1.7e-7 | 1.7e-7 |
| 1e-10 | 4.3e-7 | 4.3e-7 | 6.2e-7 | 6.2e-7 | 7.9e-8 | 4.7e-7 |

The correction stops improving on least squares at `delta = 1e-4` at `cond(Sigma) = 1e4` and at `1e-6` at 1e6; at 1e8
it still improves it by a factor of 290 at `1e-7`, by 3 at `1e-8`, not at `1e-9` and by 6 at `1e-10`. Where it does not,
the two are within a factor of 1.3 of each other: the error is that of the point, not of the multiplier
estimate. `RiskParity.NearMinimumVolatilityKeepsLeastSquaresMultipliers` pins the statuses, the KKT band, the
multipliers within a rounding band of the least-squares fit and the six refusals of `sample_data(5, 13)`; no test pins
this table.

**The safeguard accepts corrections worse than least squares near the minimum.** Of the 168 results of the table, 24
corrections were accepted that leave the multipliers further from the refined point than least squares by more than
5% (19 of them by more than 50%), all at `delta` from `1e-4` to `1e-10`: by factors up to 698 (`cond(Sigma) = 1e4`,
`n = 5`, seed 13, `delta = 1e-5`: `1.2e-11` against `1.7e-14`), 645 (`1e8`, `n = 8`, seed 11, `delta = 1e-9`: `7.5e-8`
against `1.2e-10`), 161 (`1e8`, `n = 8`, seed 13, `delta = 1e-9`: `7.3e-8` against `4.5e-10`), 17.9 (`1e4`, `n = 5`,
seed 11, `delta = 1e-4`: `3.2e-12` against `1.8e-13`), 9.0 (`1e4`, `n = 5`, seed 7, `delta = 1e-5`) and 7.1 (`1e6`,
`n = 8`, seed 11, `delta = 1e-7`). The absolute errors are small (the largest, `7.5e-8`, sits where least squares
is at `1e-10`), and the safeguard cannot refuse them by the stationarity rule: on the earlier instances their
stationarity was 0.018 to 11.3 times its noise, while the 486 accepted corrections of all the measurements of this
section that leave the multipliers no further than least squares reached 7.9 (near `sigma_min`; at most 0.38 away
from it), and 456 of those were above 0.018; no measured threshold separated the two, so the factor stays that of the
gate (16). The stationarity of the 24 and the number of results the factor refuses were not measured again.

### The augmented Lagrangian

`augmented_lagrangian` (`src/riskparity/augmented_lagrangian.cpp`) is Algorithm 3.1 of Andreani, Birgin,
Martinez and Schuverdt (2007, "On augmented Lagrangian methods with general lower-level constraints", SIAM J.
Optim. 18(4)) for the single row `h3 = vol(w) - s`; the budget, the target return and the bounds are the
lower-level set of every subproblem, minimized by `descend`. The quantities that are bounded and grown are those
of the relative gap `vol(w) / s - 1`, `lambda s` and `rho s^2`, so the iterates do not change when `Sigma` and
`s` change units (`RiskParity.ScaleInvariance`: the same weights within `1e-10` at `1e-4 Sigma`).

| constant | value | source or measurement |
|---|---|---|
| first penalty `rho s^2` | 10 | measured below |
| growth of the penalty | 10, when the gap did not fall to 0.5 of the previous one | the paper (`gamma = 10`, `tau = 0.5`) |
| bound of `abs(lambda s)` | 1e8 | the largest of a converged run is `1.7e3`, and `3.8e7` at `delta = 3e-12` above |
| ceiling of `rho s^2` | 1e12 | the largest of a converged run is `1e6`; reached only when `s` is out of reach |
| tolerance of the subproblem | `1e-2`, divided by 10 per outer iteration down to the tolerance asked | the paper asks for a sequence that decreases to 0 |
| outer iterations | at most 50 | the most used by a converged run is 19 (first effort) and 34 (second) |
| stop | the descent's residual and the relative gap at the tolerance, or `is_kkt_point` at `kAugmentedStopFactor = 2` | below |
| stall | 3 outer iterations in a row at the ceiling of `rho s^2` without the gap halving: not converged (`kAugmentedStallIterations`) | below |

Measured on the 384 shapes with a target volatility of the grid below at `n` in {3, 5, 10, 20} and one seed,
from the start of the solver (`RiskParitySolveSweep.DISABLED_MeasureAugmentedLagrangian`), at the first effort
of the solver (tolerance `1e-6`, at most 200 iterations per descent), each followed by the polish:

| constraints | converged | KKT after the polish | outer iterations, mean (max) | descent iterations, mean (max) | polish steps, most | largest `rho s^2` | largest `abs(lambda s)` |
|---|---|---|---|---|---|---|---|
| Capped + return + volatility | 48 / 64 | 48 | 8.5 (13) | 41 (116) | 2 | 1e5 | 594 |
| Capped + volatility | 64 / 64 | 64 | 8.0 (15) | 40 (132) | 2 | 1e5 | 1.7e3 |
| LongOnly + return + volatility | 52 / 64 | 52 | 7.7 (16) | 59 (200) | 2 | 1e2 | 3.1 |
| LongOnly + volatility | 64 / 64 | 64 | 8.2 (17) | 121 (1228) | 2 | 1e3 | 1.9 |
| None + return + volatility | 60 / 64 | 60 | 9.2 (16) | 193 (1189) | 6 | 1e5 | 281 |
| None + volatility | 60 / 64 | 60 | 9.1 (19) | 378 (3634) | 3 | 1e4 | 90 |

At the second effort (tolerance `1e-11`, at most 2000 iterations per descent) the 128 shapes without a box all
converge (16.4 and 16.7 outer iterations in mean), and the 28 that do not (16 Capped + return + volatility, 12 LongOnly
+ return + volatility) are the shapes whose `s` is above the volatility the box allows (below). At the first effort the
first penalty decides which converge: with a first penalty of 1 all the 356 shapes that converge at the second effort
converge at the first, with more outer iterations (means 9.8 to 12.1); with 10, 348 do (8 without a box are left to
the second effort); with 100, 346, and the boxed shapes with a volatility alone need fewer outer iterations (5.0 to 6.3)
but 12 of the 128 without a box no longer converge; with 1000, 336, and 20 of the 128. On the earlier instances 10 was the
largest of the four that converged every reachable shape at the first effort; on these 1 is, and 10 costs the second
effort 8 shapes. The constant stays at 10 here: it was not part of this plan.

**The noise stop** (earlier instances). A stop on the tolerance alone never comes where the stationarity cannot go below it (the
table of the section above): on `n = 20`, `cond(Sigma) = 1e8`, no box, uneven budget, with a tolerance of
`1e-6` and 2000 iterations per descent, the iteration ran its 50 outer iterations and 94 406 descent iterations
and ended with an iterate at 0.04 of its noise. The iteration therefore also stops at the first iterate that
`is_kkt_point` accepts.

**The stop factor, the stall and the multipliers between efforts.** Measured on the hardest class of the sweep
(`RiskParitySolveSweep.DISABLED_MeasureHardestClass`: no box, target volatility, `cond(Sigma) = 1e8`, `n` in
{10, 20, 50}, seeds 1 to 3, both scales, both budgets, with and without a target return: 72 solves; total
seconds of the 72, longest, most iterations of one solve; the five rows are one run each, in one session, so the
seconds are indicative):

| augmented Lagrangian | OPTIMAL | seconds | longest | most iterations | largest residual / noise |
|---|---|---|---|---|---|
| stop at 16, no stall (before) | 72 | 13.6 | 0.8 | 54 391 | 15.37 |
| stop at 2, no stall | 72 | 15.0 | 0.8 | 58 390 | 1.95 |
| stop at 16, stall `k = 3` | 72 | 11.5 | 0.6 | 42 065 | 15.37 |
| stop at 2, stall `k = 3` (kept) | 72 | 13.1 | 0.6 | 46 064 | 1.95 |
| stop at 2, stall `k = 1` | 72 | 12.5 | 0.6 | 47 399 | 1.95 |

- **Stop at 2, gate at 16.** Every run that was `OPTIMAL` with the stop at 16 is `OPTIMAL` with the stop at
  2 (72 of 72 here; on the whole sweep below the stop at 2 gives 2770 of 2880, and the stop at 16 was not run on
  it again: 2788 both times on the earlier instances); the cost is about 7% more iterations.
- **The stall.** Runs that converge do reach the ceiling of the penalty in their first effort on this class;
  there the gap stops halving, and ending the effort after `k` such outer iterations lets the second effort
  start afresh sooner. `k = 1` and `k = 3` give the same statuses and nearly the same cost; 3 is kept as the
  less aggressive of the two. The stall does not cap iterations: the most used by an `OPTIMAL` is still
  46 064. The 110 `MAX_ITER` of the sweep were already cheap (0.15 s for all of them; `s` out of reach of a box,
  where the penalty reaches its ceiling within a few outer iterations of small descents), so the stall changes
  their cost by nothing measurable (`AugmentedLagrangian.UnreachableTargetDoesNotConverge`: at most 15 outer
  iterations, against the 60 allowed, on its two-asset case).
- **Not carried: the multiplier and the penalty between efforts** (earlier instances; not remeasured). Starting the
  second effort from the `(lambda, rho)` where the first stopped, instead of from `lambda = 0` and the first penalty,
  turned 9 or 10 of the 72 `OPTIMAL` into `MAX_ITER` (all `n = 50`; carrying only `lambda` did the same; 63, 62 and
  62 `OPTIMAL` in the three variants), because the second effort then begins near the multiplier and the penalty at
  which the first one stalled. The restart stays.
- **Above `n = 50`.** The same class at `n = 100`, one seed (8 solves): all `OPTIMAL`, largest residual 1.84 of
  the noise, 4.4 to 7.4 s and 13 266 to 23 362 iterations; about 18 times the longest solve at `n = 50` (0.40 s).

### The policy of the solve and its sweep

`solve_risk_parity` (`src/riskparity/risk_parity.cpp`): the Spinu point; the answer when it is inside the box
and no target is set. Otherwise the minimum-variance QP over the linear rows and the bounds (`INFEASIBLE` when
it is, or when `s < sigma_min (1 - 1e-12)`; the minimum-variance portfolio when `s = sigma_min`). Then, from
the point of the polyhedron nearest to the Spinu point (one QP), two efforts, each a descent (the augmented
Lagrangian with a target volatility) followed by the polish, the second from where the first stopped:

| effort | tolerance of the descent | iterations per descent |
|---|---|---|
| first | 1e-6 | 200 |
| second | 1e-11 | 2000 |

and `is_kkt_point` on the polished point decides: `OPTIMAL`, or, after the second effort, `MAX_ITER` with NaN
weights, multipliers and residual. The short first effort is what bounds the cost: on the earlier instances, with 2000 iterations per
descent at both efforts the grid below at `n <= 20` and one seed (768 runs, measured before the rule on the
budget row) gave the same count of each status in every row, with a longest solve of 1.00 s against 0.12 s
and a mean of 2212 against 584 iterations on the shapes with a target volatility alone and no box.

**Sweep** (`RiskParitySolveSweep.DISABLED_MeasureSweep`, run with `--gtest_also_run_disabled_tests`; 24 s of
solves, of which 17 s on the shapes with a target volatility and no box, nearly all of it at `n = 20` and 50; the
checks of the results come on top, 145 s in all). Every combination of `cond(Sigma)` in
{1e2, 1e4, 1e6, 1e8}, scale {1, 1e-4}, equal and uneven budgets, box None, LongOnly and Capped, with and without
`target_return` and `target_volatility` (`s` between `sigma_min` and `1.5 sigma_min`): 192 shapes, at `n` in
{3, 5, 10, 20, 50} and seeds 1 to 3, 2880 calls of `po::solve`. Every `OPTIMAL` is checked by
`expect_risk_parity_kkt` at the band of its `cond(Sigma)`; "above the band" counts those whose stationarity
passes by its rounding floor and not by the band (the section on the noise), "KKT failed" any failure of the
check:

| constraints | OPTIMAL | MAX_ITER | above the band | KKT failed | iterations, mean (max) | seconds, total (longest) |
|---|---|---|---|---|---|---|
| Capped | 240 / 240 | 0 | 0 | 0 | 25.3 (68) | 0.30 (0.020) |
| Capped + return | 240 / 240 | 0 | 0 | 0 | 22.9 (76) | 0.31 (0.014) |
| Capped + return + volatility | 174 / 240 | 66 | 0 | 0 | 65.7 (404) | 0.86 (0.086) |
| Capped + volatility | 240 / 240 | 0 | 0 | 0 | 64.1 (356) | 0.99 (0.082) |
| LongOnly | 240 / 240 | 0 | 0 | 0 | 13.5 (37) | 0.01 (0.001) |
| LongOnly + return | 240 / 240 | 0 | 0 | 0 | 23.5 (51) | 0.19 (0.005) |
| LongOnly + return + volatility | 224 / 240 | 16 | 2 | 0 | 114.4 (1217) | 2.13 (0.195) |
| LongOnly + volatility | 240 / 240 | 0 | 26 | 0 | 141.4 (1246) | 1.87 (0.129) |
| None | 240 / 240 | 0 | 0 | 0 | 13.5 (37) | 0.01 (0.000) |
| None + return | 212 / 240 | 28 | 0 | 0 | 30.6 (87) | 0.27 (0.032) |
| None + return + volatility | 240 / 240 | 0 | 72 | 0 | 1611.3 (46064) | 9.26 (0.670) |
| None + volatility | 240 / 240 | 0 | 118 | 0 | 1559.4 (43618) | 7.76 (0.584) |

| n | OPTIMAL | MAX_ITER | above the band | iterations, mean (max) | seconds, total (longest) |
|---|---|---|---|---|---|
| 3 | 520 / 576 | 56 | 17 | 52.8 (778) | 0.10 (0.002) |
| 5 | 548 / 576 | 28 | 36 | 73.2 (2220) | 0.17 (0.006) |
| 10 | 562 / 576 | 14 | 61 | 154.0 (3890) | 0.54 (0.026) |
| 20 | 566 / 576 | 10 | 52 | 1012.1 (46064) | 8.44 (0.670) |
| 50 | 574 / 576 | 2 | 52 | 262.4 (4618) | 14.73 (0.435) |

No result is `INFEASIBLE` (the shapes are built feasible) and none fails the check.
The 110 `MAX_ITER`:

- 82 with a target volatility and a target return, LongOnly (16) or Capped (66). On each of them the augmented
  Lagrangian, run again from the start of the solver with 2000 iterations per descent, ends not converged with
  `vol(w)` below `s` by at least `4.4e-3 s`: `s` is above the volatility of the vertex the iteration reaches. The
  largest volatility of a polyhedron has no cheap certificate, so this is `MAX_ITER` and not `INFEASIBLE`
  (`RiskParity.UnreachableVolatilityIsMaxIter`).
- 28 with a target return and no box, the unbounded line of the section on the noise.

The descents that did not converge in the sweep of the polyhedral descent above (66 of 2000, with a box or at
`cond(Sigma) = 1e8`) do not appear here: every shape without a target volatility is `OPTIMAL` except those 28.
The start is not the same (the point nearest to the Spinu point, not the minimum-variance point). Run directly,
the polish turns all 66 of those descents into points that `is_kkt_point` accepts
(`RiskParityDescentFuzz.DISABLED_MeasurePolishOfUnconverged`). On the earlier instances it turned 66 of 70: the other 4
were two shapes at both scales, `cond(Sigma) = 1e8`: `n = 10`, no box and no target (which the solver answers with
the Spinu point and never descends), and `n = 20`, Capped, equal budget, seed 31696, which through `po::solve` was
`MAX_ITER` with every number NaN at both scales. That is a bounded problem without a target that the solver did not
solve: **a known limitation**, documented and not fixed. On the instances of plan 8 the shape of seed 31696 is
`OPTIMAL` at both scales (`RiskParityFuzz.UnrepairedDescentShapesThroughSolve`), and no instance of the repository
reproduces the limitation now; the solver is a local method, and a bounded problem without a target can still end
`MAX_ITER` without a point.

**Cost.** One solve at `n = 20`, `cond(Sigma) = 1e8`, no box and a target volatility takes up to 0.67 s and
46 064 iterations (0.44 s and 4618 at `n = 50`; 7.4 s and 23 362 at `n = 100`): every subproblem of the augmented
Lagrangian is a descent that converges slowly there (306.2 iterations on average for the descent alone at
`cond(Sigma) = 1e8` without a box, in the table above), and the second effort starts the multiplier method again
(carrying it over loses correct answers, above). On the earlier instances the longest solve was 9.2 s and 55 323
iterations at `n = 50`. At `n <= 10` the longest solve of the sweep is 0.03 s.

**Fuzz battery** (`tests/unit/fuzz/risk_parity/risk_parity_fuzz_test.cpp`), the cases the sweep does not
cover:

- `RiskParityFuzz.AgainstScanReference`: `n = 3`, `cond(Sigma)` 1e2 to 1e8, the three boxes, a target return
  or a target volatility, equal and uneven budgets (96 shapes; 93 `OPTIMAL` with a scan to compare). `f` of the
  result is at most the scan's global minimum times `1 + 1e-9`, the scan's value first moved to the target by
  `sum_k lambda_k h_k(w_scan)` with the multipliers of the result: the scan's points lie on the ellipse only
  to its rounding, and on the earlier instances one of them, `1.8e-10 s` inside it, had an `f` lower than an exact
  solution 1.6e-10 away by `8.2e-9` of `f`. A result above that is accepted only as another local minimum: a KKT point at the
  band, more than `1e-6` from the scan's point, with no feasible point of the curve within 20 scan steps below
  its `f`. **18 of the 93** are such local minima, with `f` from 3.5% to 960 times above the global one (12 on the ellipse of a
  target volatility, with no box or long only; 6 on the line of a target return, with no box, long only or capped); the
  test fails if there are more.
- `RiskParityFuzz.UnevenBudgetsAndConditioning`: `b_i` log-spaced from `1e-6` to 1, `cond(Sigma)` 1e2 to 1e8,
  scales 1 and `1e-4`, the budget row alone and with the Capped box, `n` in {5, 20}, 3 seeds: 96 of 96
  `OPTIMAL`, all passing `expect_risk_parity_kkt` at the band; any other status must be `MAX_ITER` or
  `INFEASIBLE` with every number NaN. `Spinu.UnevenBudgetsConditioningSweep` covers Spinu alone; the solve sweep
  covers uneven budgets with the targets at `n <= 10` in the enabled grid; this adds the extreme budget through
  `po::solve` with a box at `n = 20`. The rule "no point without `OPTIMAL`, no `OPTIMAL` outside KKT" is checked
  by `RiskParitySolveSweep.OptimalRunsAreKkt` on 576 shapes and by both tests here.
- **Forward error, the thinnest class**: the weights at `cond(Sigma) = 1e8` with a target volatility are a
  conditioning limit, measured and pinned in the section on the noise above
  (`RiskParityFuzz.WeightsAtCondition1e8BelowTheirMeasuredCeiling`); the multipliers, with and without a box, in
  the same section (`RiskParityFuzz.TargetVolatilityMultipliersByCondition`,
  `RiskParityFuzz.TargetVolatilityBoundMultipliersByCondition`).

## Fuzz battery

`tests/unit/fuzz/` holds a generator of problem shapes and references computed by another route, used by the
`MaxReturnFuzz`, `LeverageFuzz` and `RatioFuzz` tests.

- **Shapes** (`random_problems.hpp`). `fuzz_shapes(count, seed)` cycles through every combination of `cond(Sigma)` in
  {1e2, 1e4, 1e6, 1e8}, the annual and daily scales (`Sigma` times 1 or 1e-4, `mu` times its square root), four mean
  patterns (spread in [0.03, 0.15]; equal; nearly equal, with a relative spread log-uniform in [1e-8, 1e-5]; below
  the risk-free rate: `rf` above every mean with cash, every mean negative without), four boxes (none, long-only,
  tight with `upper = 2 / n`, and long-only with a third of the weights fixed by `lower == upper`) and cash, so any
  256 consecutive shapes cover them all; `n` (3 to 24, or to a smaller `largest_n`) and the covariance seed are
  random. `with_leverage` gives each shape a leverage uniform in [0.8, 3] drawn from the shape's seed, so the shapes
  themselves are unchanged (below 1 without cash it is infeasible).
- **Reference** (`frontier_reference.hpp`). `max_return_reference` is the largest target return `r` with
  `MIN_VARIANCE(target_return = r)` feasible and of volatility at most `s`, by doubling and bisection on `r` to
  1e-13 relative; it returns the return of that portfolio, not `r`, because the solver accepts a target-return row to
  its feasibility tolerance, with a floor of 1e-12 that is 1e-9 of a daily mean. NaN when no feasible portfolio has
  volatility at most `s`.
- **Sign patterns** (`sign_pattern_reference.hpp`). `sign_pattern_optimum` is the optimum of a QP under
  `sum_i |x_i| <= L` by brute force for n <= 8: on the orthant of each sign pattern `s`, `s_i x_i >= 0`, the
  constraint is the linear row `s'x <= L`, and the optimum is the best `Optimal` value over the `2^n` patterns, with
  explicit rows only.
- **Ratio reference** (`ratio_reference.hpp`). `max_ratio_reference` maximizes the ratio along the efficient branch of
  `MIN_VARIANCE(target_return = r)`, never through the homogeneous QP. With cash fixed at 0 the ratio problem is the one
  without cash on the risky assets with `1'w = 1` and mean `v` (`mu - rf 1`, or `sigma` with `rf = 0`), so the
  reference solves that one, with the risky bounds and the leverage. From the minimum-variance return it doubles
  `r - r_mv` until the ratio stops increasing or no portfolio reaches `r`; there it bisects on feasibility for the
  largest return, and then it runs a golden-section search on `r / vol(r)` to 1e-13 relative, keeping the best ratio of
  a solved portfolio. Whether the ratio still increases at a point is the sign of
  `lambda_budget + beta'w + lambda_lev ||w||_1`, read from the frontier point's multipliers: by the multiplier
  convention `w' grad f = 2 vol^2 = lambda_budget + lambda_ret r + beta'w + lambda_lev ||w||_1`, so
  `d/dr (r / vol(r)) = (2 vol^2 - r lambda_ret) / (2 vol^3)` has that sign. It is a sum of terms of the order of `r`,
  where the difference of two ratios at a large `r` would be lost in their rounding. The result is `INFEASIBLE` when
  `MIN_VARIANCE` is, when the best ratio is not positive, or when the cash bounds exclude 0. It is `UNBOUNDED` when the
  ratio still increases after 30 doublings, about 1e9 times the scale of the returns, with weights of 7.5e8 to 5.0e9
  there. Measured on the battery, every peak was found by the second doubling (`r - r_mv` at most 3 steps) and every
  largest return by the first. Relative to `2 vol^2 / r`, the slope was at most -4.7e-4 where the reference stopped past
  a peak, and at least 5.9e-4 at the last doubling of the 28 `UNBOUNDED` cases. The weights of the `OPTIMAL` results
  are at most 3.2.
- **Bands.** KKT residuals (`expect_user_kkt`) and the agreement of the return (relative to the larger of the
  reference and `|mu|max`) are held to 1e-10 for `cond(Sigma) <= 1e4`, 1e-9 up to 1e6 and 1e-7 up to 1e8, by the
  condition number of the covariance actually solved; above 1e8 the 1e-7 band applies, except to the stationarity of
  problems with many interior weights, which is a conditioning limit and is not checked there (see
  [the numerical decisions](#numerical-decisions-of-max_sharpe-and-max_diversification)).
- **`MaxReturnFuzz.AgainstFrontierReference`**: `fuzz_shapes(400, 1)` without the fixed box (304 shapes), caps at 0.1,
  0.3, 0.5, 0.7 and 0.9 of the way from `sigma_min` to the largest asset volatility: `OPTIMAL` (or `INFEASIBLE` exactly
  when the reference is NaN), the return of the reference, `expect_user_kkt` and `vol <= s (1 + tol)`, all within the
  band. All 1520 cases pass, 779 of them with an inactive cap, with no case pinned; the worst agreement with the
  reference is 1.1e-12 at `cond <= 1e4`, 6.3e-12 at 1e6 and 1.8e-9 at 1e8.
- **`MaxReturnFuzz.NeverSilentlyWrong`**: the nearly-equal-mean shapes of `fuzz_shapes(400, 2)` (100) with a nearly
  riskless asset appended (variance log-uniform in [1e-10, 1e-6], mean 0.4 times the smallest), caps at 0.1, 0.5 and
  0.9 of the way to the largest asset volatility and at 1.01, 2 and 10 times `sigma_min`. The covariance solved has a
  condition number from about 6e2 to 2e16, so most cases take the 1e-7 band. The result must be `OPTIMAL` with the
  return of the reference, every weight within its bounds and `vol <= s (1 + tol)`, or `MAX_ITER`. Of the 528 cases,
  440 are `OPTIMAL` (worst agreement 1.1e-8, at a band of 1e-7) and 88 `MAX_ITER`.
- **`LeverageFuzz.AgainstSignPatterns`**: `fuzz_shapes(256, 4, 8)` (every combination once, n from 3 to 8) with
  `with_leverage`, for `MIN_VARIANCE`, `UTILITY`, `IMPORTANCE`, `KELLY` and `TRACKING_ERROR` (`MAX_RETURN` is checked by
  `MaxReturnFuzz.LeverageFuzz`), each at the shape's leverage and, where the optimum without leverage is above the least gross exposure the
  budget allows (1 without cash, 0 with cash), at a leverage a fraction in [0.2, 0.9] of the way down to it.
  `TRACKING_ERROR` runs over `T = n` periods `R = L' + rf` (`L` the factor of `Sigma`, so `R~'R~ = Sigma` and the
  band is that of the shape). `INFEASIBLE` exactly when no pattern is `Optimal`; otherwise `OPTIMAL` with the value of
  the reference (relative to the larger of the magnitude of its terms and `|g|max |w|max` over every weight, cash
  included, since the optimum at all cash is 0 to rounding) and `expect_user_kkt`, within the band. Of the 1979
  cases, 687 bind (59 to 180 per objective, 356 with two or more active generated rows, up to 8), 1057 do not and 235
  are `INFEASIBLE` (L < 1 without cash, or below what the cash cap of the tight box or the fixed weights force); the
  worst agreement is 4.2e-14 at `cond <= 1e4`, 3.9e-12 at 1e6 and 9.7e-11 at 1e8, at most 3.9e-3 of the band.
- **`MaxReturnFuzz.LeverageFuzz`**: `with_leverage(fuzz_shapes(400, 3))` without the fixed box (304 shapes), at the
  caps of `AgainstFrontierReference` from `sigma_min` under the leverage (the largest asset volatility alone where no
  portfolio is feasible under it), each at the shape's leverage and, where the `MAX_RETURN` optimum without leverage
  is above the least gross exposure the budget allows, at a leverage a fraction in [0.2, 0.9] of the way down to it.
  The reference is `max_return_reference` with the leverage in its `MIN_VARIANCE` (the frontier volatility still
  increases above the minimum-variance return, the leverage being convex); the checks are those of
  `AgainstFrontierReference`. Of the 2316 cases, 2170 are `OPTIMAL` (1131 with an inactive cap) and 146 `INFEASIBLE`;
  1066 bind the leverage, 760 of them with the cap. The worst agreement with the reference is 1.7e-12 at
  `cond <= 1e4`, 2.2e-12 at 1e6 and 1.6e-10 at 1e8. `gcc-debug` runs the first 256 shapes (every combination once) at
  the caps 0.1, 0.5 and 0.9: 412 cases bind the leverage, 277 with the cap (counted in `gcc-release` on those shapes;
  `gcc-debug` was not run in plan 8).
- **`RatioFuzz.AgainstFrontierReference`**: `fuzz_shapes(300, 3)` for `MAX_SHARPE` and `MAX_DIVERSIFICATION`, at each
  shape without a leverage and with `with_leverage`. The status must be that of `max_ratio_reference`. An `OPTIMAL`
  result must have the reference's ratio, the cash weight exactly 0 and `expect_user_kkt` with the gradient written in
  the test, `-v / vol + (v'w) Sigma w / vol^3`, and cash held to the bounds [0, 0]; all within the band. A result
  without a point has every weight and bound multiplier, the budget and leverage multipliers and the residual NaN. All 1200 cases
  pass, with no case pinned. 1039 are `OPTIMAL` (80 bind the leverage), 133 `INFEASIBLE` (`rf` above every mean with a
  box, and a leverage below 1) and 28 `UNBOUNDED` (24 `MAX_SHARPE`, 4 `MAX_DIVERSIFICATION`, all without a box or a
  leverage). The worst agreement of the ratio is 3.2e-12 at `cond <= 1e4`, 1.8e-11 at 1e6 and 2.1e-9 at 1e8 (at most
  0.032 of the band). The worst residual of `expect_user_kkt`, the stationarity, is 9.9e-13 at `cond <= 1e4`, 7.5e-11
  at 1e6 and 6.4e-9 at 1e8 (at most 0.075 of the band). Every other measure stays below 1e-5 of its band.
- **The known residual (history).** Until the tangent crossing, one case of `MaxReturnFuzz.LeverageFuzz` missed its
  band and was pinned by name: shape seed 3963545680699380072 (n 24, `cond(Sigma) = 1e8`, annual scale, nearly equal
  means, no box, cash, the shape's leverage 2.154), cap fraction 0.7 (cap 858.4). The cap is reached at `t = 6.7e15`
  on the leverage face, where `vol(Q(t))` varies by 1.5e-7 relative between `t` that differ in the last digits; the
  bisection ended at adjacent doubles, whose lower end was accepted with `vol = s (1 - 1.006e-7)`, a volatility
  complementarity of 1.006e-7 against the band of 1e-7. The tangent crossing solves it with nothing relaxed: `OPTIMAL`,
  `vol = s (1 - 4.4e-16)`, the return of the reference to 5.9e-14 relative, `lambda_vol = -1.2888e-13`, and the full
  band of `expect_user_kkt`. The shape's covariances on Windows and on Linux are frozen in
  `tests/unit/fuzz/known_shape_fixture.hpp` and checked at every cap of the fuzz
  (`MaxReturnFuzz.KnownShapeUnder{Windows,Linux}Covariance`): every cap is `OPTIMAL`, `|vol / s - 1| <= 5.6e-16`, the
  return within 6.6e-14 of the reference. Under the Linux covariance two caps needed one more fix: a solve of `Q(t)` of
  the climb (`t = 9.2e12`) took 556 iterations of the QP solver against a limit of 490, because generated rows leave
  and are generated again; the leverage family now counts as `2n` rows in the iteration limit
  (`GoldfarbIdnani.GrossExposureFamilyWithManyExchanges`). The case ran in `gcc-release` only (the debug run has no
  cap fraction
  0.7). Task 0's pinned case of `AgainstFrontierReference` (shape seed
  11394049047870518593, cap fraction 0.9, cash, long-only, nearly equal means, `cond(Sigma) = 1e8` daily, `MAX_ITER`
  where `Q(t)` was lost near `t = 1.2e12`) is solved by the shift along the cash bound and is no longer pinned.
- **Other platforms.** The generators and the library give the same bits on Windows (MinGW and MSVC), Linux and macOS,
  so a test instance is the same instance everywhere. Verified locally on Windows (MinGW and MSVC: the full suite, 371
  tests) and on Linux (the `manylinux_2_28` image, GCC 14.2.1, glibc 2.28, Release: the full suite, 371 tests, with the
  bits committed from the MinGW build); macOS is checked by the same tests in CI. In `src/` no libm function is called besides `sqrt` and `fma`
  (both rounded exactly): `po::portable::exp`, `log` and `hypot` (`src/core/portable_math.cpp`) replace them, the
  Givens of the QP and the Spinu objective being their users. The test generators use the same functions, Marsaglia's
  polar method for the normals and a rejection method for the integers, and no `<random>` distribution. The claim is
  checked by `Generators.HaveTheCommittedBits` (the instances, by their bits), `PortableMath.*` (the three functions
  against committed values) and `Determinism.*` (the solution bits of one problem of each objective), which run in
  CI on Windows, Linux and macOS. The change of the QP solution bits that the portable `hypot` made, and the run times, are in
  [docs/performance/portable-math.md](../performance/portable-math.md). Before plan 8 the generators built `Sigma` with
  `std::log`, `std::exp` and `std::cos`, which glibc rounds differently from the Windows UCRT: on Linux 561 of the 576
  covariance entries of the known residual's shape differed, by up to 1856 ulp (3.6e-13 relative; `mu` was the same),
  and with that covariance the case was `MAX_ITER` at cap fractions 0.3 to 0.9, on Linux and on Windows alike, before
  the tangent crossing; both covariances are now `OPTIMAL` at every fraction (the known residual above). The Linux
  covariance of the fixture comes from the `manylinux_2_28` image (GCC 14.2.1, glibc 2.28). Before the tangent crossing
  `MaxReturnFuzz.LeverageFuzz` failed on Linux and macOS; `RiskParityDescentFuzz.ConvergedRunsAreKkt`,
  `RiskParityReference.ForwardErrorMeasuresTheDistanceToTheRefinedPoint` and `MaxReturn.FaceShiftCashMultiplier` failed
  on Linux, and `RiskParity.TargetVolatilityEllipseAgainstScan` on macOS (its solve landed below the value of the scan,
  a reference of limited resolution). Moving to the portable generators changed every test instance; the tests that
  relied on the old ones were triaged one by one, and the counts and floors that moved (18 other local minima
  in `AgainstScanReference`, 540 `OPTIMAL` results in `OptimalRunsAreKkt`, 623 converged runs in `ConvergedRunsAreKkt`)
  are those of this document.
- **Runtime.** In `gcc-release` (`ctest`, 21 s for the 371 tests of the suite) `MaxReturnFuzz.AgainstFrontierReference`
  0.78 s, `NeverSilentlyWrong` 0.42 s, `MaxReturnFuzz.LeverageFuzz` 3.4 s, `LeverageFuzz.AgainstSignPatterns` 0.68 s and
  `RatioFuzz.AgainstFrontierReference` 1.9 s (1.3 s, 0.8 s, 6.5 s, 0.9 s and 3.1 s on the earlier instances, with the
  libm functions that the portable ones replaced). `gcc-debug` was not run in plan 8; its figures are those of plan 7 on
  the earlier instances: 25.0 s, 11.9 s, 37.7 s (117 s in full, hence the smaller debug run above), 11.6 s and 43.7 s
  (in full). The other tests run their full counts in both presets.

**Plan-2 scripts against the fixed DLL.** The Python probes of the plan-2 review (`MAX_RETURN` against Clarabel, an
SOCP reference, and an independent KKT check) were run against the DLL before and after these changes (plan 2; not repeated since). Their
thresholds are stricter than the bands: stationarity 1e-8, feasibility 1e-9, return gap 1e-9.

- The re-review cases (30 assets and a nearly riskless one of variance 1e-10, nearly equal means, long-only with
  `upper = 0.3`, `s = 10 sigma_min`): on HEAD, `OPTIMAL` with a weight of -0.031 below its bound of 0 and, in the other
  case, a return of 0.0800162 against 0.0800174; now both agree with the reference (lowest weight -1.4e-18, return
  within 2.4e-12 relative).
- The probes of fixed weights (`lower == upper`), slow motion of `w` and a nearly riskless asset at its cap: cases
  beyond the scripts' thresholds fell from 27 and 34 on HEAD to 5 and 5 (the second probe repeats the first). The 5:
  - four cases with n = 40, `cond(Sigma) = 1e8`, `s = 2 sigma_min` and `mu` spread from 1e-4 to 1e-7 relative, where
    `vol / s` is 1 + 1.5e-9, 1 + 2.2e-9, 1 + 1.7e-9 and 1 + 1.9e-9 (feasibility and volatility complementarity of the
    same size; stationarity at most 5.6e-13, return gap at most 7.3e-12). That is inside the band of 1e-7 for
    `cond(Sigma) = 1e8`, the rounding of the crossing. On HEAD these cases returned `lambda_vol = 0` with a
    stationarity of 1e-7 to 1e-4: silently wrong.
  - one case with fixed weights (n = 30, one weight fixed) at `s = sigma_min (1 + 1e-9)`: return gap 2.1e-5 against a
    Clarabel solution marked `AlmostSolved`, with stationarity 9.2e-12 and feasibility 3.5e-18; the same gap as on
    HEAD, from the reference.
- `MAX_ITER`: 1 case (n = 30, variance 1e-10, means spread by 1e-9, `s = 10 sigma_min`), also `MAX_ITER` on HEAD, where
  Clarabel reports `AlmostSolved`; HEAD had 3.
- The broader fuzz scripts are unchanged or better: 300 of 300 correct with a nearly riskless asset; 400 of 400
  `OPTIMAL` against Clarabel, worst stationarity 2.3e-9 (2.6e-9 on HEAD) and the same entries flagged against
  `AlmostSolved` references; 600 of 600 with cash; 240 of 240 `OPTIMAL` per condition number from 1e2 to 1e8; the
  equal-mean cases at `sigma_min`; 300-asset and 120-asset sample covariances all `OPTIMAL` (longest solve 0.82 s, 0.76 s
  on HEAD); `TRACKING_ERROR` identical. One script, `repro5.py`, stops with an error of its own: it passes
  `max_volatility = 0` for a case with cash (`sigma_min = 0`), which the setter rejects, on HEAD too.

## Sources of RISK_PARITY and of the earlier numerical decisions

Verified means the text was read for the claim it supports. The full list, with the pages, is in
`docs/research/2026-11-16-risk-parity-and-numerics.md`; [math.md](math.md#sources) repeats it by objective.

- Bruder and Roncalli (2012, section 2.3.3): the positive-orthant solution exists and is unique. Verified.
- Bai, Scheinberg and Tütüncü (2016), lemmas 2.1 to 2.3 and examples 3.1 and 3.2: at most `2^(n-1)` normalized solutions,
  and an indefinite Hessian. Verified on the 2013 preprint; the published version was not read.
- Andreani, Birgin, Martinez and Schuverdt (2007), algorithm 3.1: the augmented Lagrangian of the volatility. Verified.
  Its penalty rule (halving test, growth factor 10) is the one used; the plan had attributed a different rule
  (a factor 4 instead of 2) to Nocedal and Wright, section 17.4.
- Feng and Palomar (2015, SCRIP) and Spinu (2013): not read. The convergence theorem of SCRIP is not verified; Spinu's
  algorithm is taken from Griveau-Billion, Richard and Roncalli (2013), which was read.
- Ogita, Rump and Oishi (2005): the compensated dot product. Not verified against the paper; its effect is measured above.
- Tibshirani (1996, section 6): the family of the leverage rows. Verified. Niedermayer and Niedermayer (2007), proposition
  1 and corollary 1: the path of `MAX_RETURN`. Verified.
- Higham's rounding-error analysis of Cholesky: not verified in the book (above).

## Where the tests are

| Suite | What it checks |
|---|---|
| `tests/unit/linalg/cholesky_test.cpp` | The factor and the triangular solves; rejection of an indefinite matrix and of a pivot at rounding level |
| `tests/unit/linalg/qr_test.cpp` | Householder QR: `Q'a = R`, the lower factor of `a'a`, least squares, rank deficiency (duplicated, combined and zero columns), and a rank threshold relative to the columns (`cond = 1e10` full rank, `1e15` not) |
| `tests/unit/qp/goldfarb_idnani_test.cpp` | Unconstrained, equality-only and bounded problems; infeasible and inconsistent ones; the scale tests above; a provided factor of G gives the same solution; feasibility at the final x from a start of 3e17; the gross exposure family against the sign-pattern brute force (100 problems, 85 with the family active), inactive (bit for bit the problem without it) and in the homogeneous form, offset `L 1` and right-hand side 0; a face of the family held as an equality at `t = 1e10` (the family does not regenerate it) |
| `tests/unit/core/generators_test.cpp` | The test generators by their bits (the normals, `exp10`, the integers, `random_covariance`): a test instance is the same on every platform |
| `tests/unit/core/portable_math_test.cpp` | `po::portable::exp`, `log` and `hypot` against committed values, each checked against a 60-digit reference, branch by branch |
| `tests/unit/core/data_test.cpp` | Rejection of non-finite, asymmetric, indefinite, singular (a duplicated asset, a sample covariance of fewer periods than assets) and mismatched input; symmetrization within tolerance; the stored factor |
| `tests/unit/problems/min_variance_test.cpp` | KKT in the user's variables, the closed forms, the frontier and the capital market line, per-asset bounds, cash bounds, an ill-conditioned covariance |
| `tests/unit/problems/utility_test.cpp` | `UTILITY`, `IMPORTANCE` and `KELLY`: closed forms, defaults, KKT, and the identities between them and `MIN_VARIANCE` |
| `tests/unit/problems/tracking_error_test.cpp` | Exact replication, constrained least squares, cash, target return and long-only KKT, a rank-deficient series rejected, an ill-conditioned one (`cond(R) = 1e8`) accepted with the objective of an independent least squares, no covariance needed |
| `tests/unit/problems/max_return_test.cpp` | Frontier and capital market line, an infeasible cap, the cap at the minimum (including a vertex of the box), an inactive cap, equal means, the multiplier against finite differences, daily scale, random long-only problems, a conditioning sweep (`cond(Sigma)` 1e2 to 1e7), a risk-free rate above every mean (all cash), equal means with tight upper bounds or with cash, the budget multiplier under a common shift of the means, the leverage along the path (the vertex of the leverage in closed form, both caps binding against finite differences in s and L), and the multipliers mapped back from the shift along the leverage and the cash bound in closed form |
| `tests/unit/problems/leverage_test.cpp` | The leverage in `MIN_VARIANCE` and `UTILITY`: KKT with the subgradient term, a weight held at 0 by a subgradient inside [-1, 1], `L < 1` without cash infeasible, cash not counted, the multiplier against a finite difference, NaN without a leverage and 0 when inactive |
| `tests/unit/problems/max_sharpe_test.cpp` | `MAX_SHARPE` and `MAX_DIVERSIFICATION`: the tangency and the most diversified portfolio in closed form, the cash weight exactly 0 with `beta_0` (0 at the tangency, not under an active cap; bit for bit the no-cash weights for `MAX_DIVERSIFICATION`), cash bounds that exclude 0 and `rf` above every mean long-only `INFEASIBLE`, `1' Sigma^-1 mu~ < 0` without bounds `UNBOUNDED` with NaN weights and multipliers, the boundary `1' Sigma^-1 v = 0` exact in the data (`rf` at the minimum-variance return and a few ulps away for `MAX_SHARPE`, a covariance with `1' Sigma^-1 sigma = 0` for `MAX_DIVERSIFICATION`, over `cond(Sigma)` 1e2 to 1e8 and both scales: `UNBOUNDED`, or `OPTIMAL` against an extended-precision reference), the same boundary with three assets bounded and their lower rows active (no exact `kappa <= 0` missed, none decided by the data `UNBOUNDED`), bounded cones at `cond(Sigma)` 1e12 and 1e14 always `OPTIMAL` (status, budget, bounds and leverage; the whole KKT for long-only at 1e12), long-only and leverage KKT, and no sharper portfolio among 10 000 feasible ones |
| `tests/unit/fuzz/ratio_fuzz_test.cpp` | `MAX_SHARPE` and `MAX_DIVERSIFICATION` over the fuzz battery, without and with a leverage, against the frontier reference: the status, the ratio and KKT |
| `tests/unit/fuzz/leverage_fuzz_test.cpp` | The five QP objectives with a leverage over the fuzz battery (n <= 8) against the sign-pattern brute force |
| `tests/unit/fuzz/max_return_fuzz_test.cpp` | `MAX_RETURN` over the fuzz battery against the frontier reference, without and with a leverage, and nearly equal means with a nearly riskless asset: correct or `MAX_ITER`, never silently wrong |
| `tests/unit/riskparity/spinu_test.cpp` | Spinu against the diagonal and the two-asset closed forms (1e-12), the contributions against an equal budget (n up to 100, `cond(Sigma)` 1e2 to 1e6, 20 seeds) and an uneven one (`b_i` from 1e-6 to 1, `cond(Sigma)` up to 1e8, both scales), in the bands of the table above |
| `tests/unit/riskparity/risk_parity_test.cpp` | `RISK_PARITY` through `po::solve`: the Spinu point (multipliers 0 and NaN, a box that holds it, a negative lower bound); a target return and a target volatility against the scan of the feasible curve (`n = 3`); two assets with a target return (the 2 x 2 system); `s = sigma_min` (`lambda_vol = -inf`, the rest NaN) and just above it; `INFEASIBLE` (volatility below the minimum, return out of range, bounds); `MAX_ITER` without a point (a volatility out of reach, a line without a minimizer); equal means with a redundant and an inconsistent target; active bounds with both targets; the identity of the multipliers at `cond(Sigma)` 1e4 and 1e6; the change of units of `Sigma`; the multipliers against central differences of the optimal `f` in `r` and `s`; `s = sigma_min (1 + delta)` for `delta` from `1e-2` to `1e-10` with the safeguard of the multiplier correction; the multipliers of the target-return path never farther from the refined point than the least-squares fit, at `1e8` |
| `tests/unit/riskparity/newton_polish_test.cpp` | The polish from a point at `1e-6` to the rounding floor in at most two steps and no step from there; the multipliers settled without a step; the result independent of the equality multipliers of the start; a bound whose multiplier has the wrong sign released; the KKT residual in noise units sees an error of `1e-9` of a weight, a weight outside the box, a multiplier of the wrong sign and a NaN; weights that do not resolve the budget row are no KKT point; the Newton correction of the multipliers against the refined KKT point at `cond(Sigma)` 1e6 and 1e8, and a correction that does not contract or that changes the side of a held bound refused |
| `tests/unit/riskparity/augmented_lagrangian_test.cpp` | The augmented Lagrangian to a KKT point at the band with and without a target return; multiplier and penalty bounded at `s = sigma_min (1 + 1e-6)`; the outer and inner caps reported; a volatility out of reach not converged |
| `tests/unit/fuzz/risk_parity/risk_parity_fuzz_test.cpp` | `n = 3` against the global minimum of the scan (18 other local minima of 93, verified and counted); budgets from `1e-6` to 1 with `cond(Sigma)` up to `1e8` through `po::solve` (96 of 96 `OPTIMAL` and KKT); the two descent shapes the polish does not repair, through `po::solve`; the value and the weights of the thinnest class at `cond(Sigma) = 1e8` against the refined point (value at the band, weights within their measured limit); the multipliers of the target-volatility path against the refined point by `cond(Sigma)` (1e4 to 1e8), without a box and with LongOnly and Capped, the bound multipliers of the held weights by value, within measured ceilings |
| `tests/unit/fuzz/risk_parity/risk_parity_multiplier_measure_test.cpp` | Disabled measurement only: the multipliers of the target-volatility path at the point `po::solve` returns, least squares and after one Newton-KKT solve, against the refined KKT point, with an optional dump for the 80-digit reference; the multipliers of the QP path (`MIN_VARIANCE` with a target return) against its active-set KKT point in `long double` |
| `tests/unit/fuzz/risk_parity/risk_parity_solve_sweep_test.cpp` | `po::solve` over 576 shapes: no status other than `OPTIMAL` carries a point, every `OPTIMAL` passes `expect_risk_parity_kkt` at the band of its `cond(Sigma)` (the stationarity at the band or within twice its rounding floor), and every one above the band agrees in weights and in `f` with the KKT point of its active set refined in `long double` |
| `tests/unit/solve/requirements_test.cpp` | The objective x setter matrix (leverage accepted by the eight objectives that take it, `MAX_RETURN`, `MAX_SHARPE` and `MAX_DIVERSIFICATION` among them; `MAX_DIVERSIFICATION` needs no `mu`; `RISK_PARITY` takes `target_return`, `target_volatility` and `risk_budget` and no cash), its messages and the order of the checks |
| `tests/unit/solve/solve_test.cpp` | Argument checks of `po::solve` and their messages |
| `tests/unit/solve/determinism_test.cpp` | The solution bits of one problem of each objective, from literal inputs, committed: the same bits on every platform |
| `tests/unit/metrics/portfolio_test.cpp` | Portfolio return and volatility, utility and tracking error, with and without cash |
| `tests/unit/capi/` | The C boundary: status codes, `po_last_error()`, null handles, the new setters and `po_data_set_returns`, ABI constants, the metrics; `c_header_check.c` compiles the header as C11 |
| `python/tests/test_min_variance.py`, `test_qp_objectives.py` | The Python package against closed forms, `numpy` and the status contract; the frozen `Spec`, the objective identities, `MAX_RETURN` on the capital market line, the input checks of `Data` |
| `python/tests/test_ratio_objectives.py` | `MAX_SHARPE` and `MAX_DIVERSIFICATION` against the closed forms (`numpy`), the cash weight 0, `UNBOUNDED` and `INFEASIBLE`, the leverage and its multiplier, the hash of a `Spec` with a leverage and array bounds, `sharpe_ratio` and `diversification_ratio` against `numpy` |
