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
    their band. `MaxReturn.BudgetMultiplierFollowsCommonMean` checks the mapping through `expect_user_kkt` and
    against means lifted by 1.
  - Without a budget (cash) the common mean falls on the inequality face that holds the portfolio once it binds: the
    cash bound `w_cash >= 0`, or the leverage, which with cash acts as a budget. The path first tries the face it
    certified last, then the active face of largest multiplier in `Q(t)` itself, and keeps `Q(t)` unshifted when the
    certificate fails (next decision). Without the shift, in `MaxReturnFuzz.LeverageFuzz` the leverage binds at `t`
    near 1 and the spread of the means (1e-8 to 1e-5 relative) moves `w` along it until the cap at `t` from 1e4 to
    beyond 1e15: 41 cases of 15 shapes, all cash with nearly equal means, failed (25 `MAX_ITER`, 16 accepted at an
    adjacent end with `|vol - s| / s` from 3.4e-10 to 6e-6, above their band; at `t = 2.6e10` on one shape
    `vol(Q(t))` moved by 1.9e-5 between `t` 2e-9 apart). With it, one case remains (below).
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
  tests (27 732 shifted solves of an inequality face, 671 with `c = 0`): 25 048 accepted and 2013 rejected, with 366
  shifted solves not `Optimal` (fallback as well). Relative to `|c|` (`c` is negative for some faces, the leverage face among them, so the
  scale is its magnitude), the face multiplier was within rounding of 0 (`|u / t| <= 1.4e-14 |c|`, a face whose
  multiplier is 0 in exact arithmetic) in 2 accepted and 3 rejected solves; both outcomes are exact to rounding there,
  the rejected ones falling back to `Q(t)` unshifted. Outside them the smallest accepted multiplier was 9.0e-4 `|c|` and
  the rejected one closest to 0 was -5.0e-5 `|c|`, ten orders above the rounding.
- **The inactive cap: the limit certificate.** The cap is inactive at a point of the climb when the point satisfies the
  KKT conditions of `max mu~'w` over the rows (every row within the QP's feasibility tolerance, no active slope below
  `-16 n eps |v|max`, and `|v + N b|max <= 16 n eps` times its largest term). Those conditions prove the point optimal
  for every larger `t`, so no comparison of weights between two points is needed. The generated rows of the leverage
  enter as ordinary active rows (their normals `a - s_j` in `N`, and the family checked by the rule of its generated
  row at `sign(w)`), and the certificate is judged in the frame of the point's shift. On the fuzz battery and the
  `MaxReturn` unit tests (8976 evaluations of the stationarity) the residual of the accepted certificates was at most 0.64 `n eps` of the
  largest term (1033 accepted, 321 of them exactly 0, with equal means) and that of the rejected ones at least
  8.6e4 `n eps`: the factor 16 leaves 25 above the first and 5400 below the second. With the leverage and the face
  shift, over the `MaxReturn` unit tests and the three `MaxReturnFuzz` tests (22 201 evaluations, 1230 of them with
  every term 0, equal means): without generated rows, accepted at most 0.99 `n eps` (1567) and rejected at least
  3.1e5 `n eps` (17 641); with generated rows active, accepted at most 2.4 `n eps` (156) and rejected at least
  8.1e10 `n eps` (1607). The factor 16 leaves 6.7 above the largest accepted and 1.9e4 below the smallest rejected. The former decision, constant
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
  active normals, by Householder least squares; a slope below `-16 n eps |v|max` is a decrease. For the generated
  rows (52 234 slopes over the same runs), the decreasing ones were at most -4.2e6 `n eps |v|max` and the others at
  least -0.40: 2.6e5 and 40 from the threshold. Differencing two
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
  `vol(t)` is continuous, so at an end that is the answer the deviation is the rounding of `vol(Q(t))`. Measured
  with an instrumented copy of the policy that accepts every end and labels it by the checks of its case:
  - **Correct ends.** 1865 ends whose result passed every check, the volatility complementarity of
    `expect_user_kkt` included: the 101 of `MaxReturnFuzz.LeverageFuzz` in the normal search, and, with every
    bisection forced to adjacent doubles, 725 of `AgainstFrontierReference` and 1039 of `LeverageFuzz`. The largest
    deviation is 1.01e-7, the known residual below (`t = 6.7e15`, `cond(Sigma) = 1e8`); the largest of the forced
    ends of `AgainstFrontierReference` is 7.98e-8, in line with `1e-15 cond(Sigma)`. The two ends of the unit suite
    in the normal search (`MaxReturn.EqualMeansWithCash`) deviate by 1.6e-15 and 2.9e-15.
  - **Wrong ends.** The 9 ends of `NeverSilentlyWrong` in the normal search, reached where `Q(t)` is lost in its
    rounding (`t` from 1.8e8 to 6.1e10), deviate by 1.6e-3 to 0.23. Accepted, they are `OPTIMAL` with `lambda_vol < 0`
    and a volatility up to 23% below `s` (that test checks the return, the bounds and `vol <= s`, which nearly equal
    means let them pass, so the policy and not the test has to reject them).
  - **The constant.** 1e-5 is 99 times the largest correct deviation (125 times the largest forced one) and 160 times
    below the smallest wrong one, a separation of 1.6e4.
  - **Noise-relative forms tried, none of which separates.** Each form `F` was divided into the deviation, and the
    separation is the smallest `dev / F` of a wrong end over the largest of a correct one (it must be above 1, and
    about 2e4 for a margin of 140 on each side): `F1`, the weight noise of the path
    `n eps (max(1, |w|max) + t |v|max)` mapped through the factor, times `||(|L'| 1)||` over `s`: 0.087;
    `F2`, the residual noise of the stationarity mapped through the factor, `n eps ||L^-1 (|Sigma||w| + t |v|)|| / s`:
    0.0026; `F2b`, the same without the `t |v|` term: 36; `F3`, `n eps (max L_jj / min L_jj)^2`, a cheap estimate of
    `cond(Sigma)`: 0.0077. The forms with `t |v|` fail because the wrong ends sit at large `t`, where that term makes
    the noise as large as their deviation; `F3` fails because they come from covariances of condition number up to
    1e16 (the nearly riskless asset), where a noise proportional to the condition number exceeds them; `F2b` orders
    the two groups but leaves a margin of 6 on each side. `kAdjacentEnd` therefore stays a constant.
  - **What it guarantees.** An accepted adjacent end is guaranteed only to `|vol - s| <= 1e-5 s`, not to the band of
    its condition number: `kAdjacentEnd` is 1e5 times the band at `cond(Sigma) <= 1e4`, where the rounding of
    `vol(Q(t))` is about 1e-13, and looser than the band of 1e-7 at `cond(Sigma) = 1e8`. The battery holds one end
    accepted between the band and 1e-5, the known residual of `MaxReturnFuzz.LeverageFuzz` (in "Fuzz battery", at
    1.006e-7). No end of `AgainstFrontierReference` or of the unit suite lies there.
- **Breakdown.** Feasibility is decided at `t = 0`; a `Q(t)` that fails at `t > 0` is `MAX_ITER`, never
  `INFEASIBLE`.

## Numerical decisions of the QP solver

Five decisions were taken while validating the Goldfarb-Idnani implementation (`src/qp/`). Each is pinned by a test.

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
solver over `GoldfarbIdnani.GrossExposure*`, the `Leverage` suite and `LeverageFuzz.AgainstSignPatterns` (1992 solves
with the family), and over 54 larger solves outside the battery (n = 50, 150 and 300; `cond(Sigma)` = 1e2, 1e6 and
1e8; `MIN_VARIANCE` and `UTILITY`; L = 1.1, 1.5 and 3; all `OPTIMAL` and within their band by `expect_user_kkt`):

- **The generated row and its violation.** `s'x` is largest at `s = sign(x)`, so `(a - s)'x >= b` with that `s` is the
  most violated row of the family, and its slack `a'x - ||x||_1 - b` is the family's. It competes with the most
  violated explicit inequality by the rule of an explicit row, applied to the generated row: violated when the slack
  is below `-1e-12 max(1, |b|, ||a - s|| |x|max)`; the explicit row wins a tie. No new constant. Over the 1992 solves,
  every generated row came from a violation of at least 8.3e4 times that tolerance (5781 generations), while the
  slack that the rule tolerated as rounding was at most 0.16 of it in the middle of a solve and 6.2e-4 of it at the
  final `x`: the rule leaves 6 above the largest rounding and 8.3e4 below the smallest violation.
- **The tie at `x_i = 0`.** `s_i = +1` where `x_i = 0` exactly; the sign of a weight within the tolerance of 0 is that
  of its rounding. Both choices give the same slack, and a sign that is wrong for a weight held at 0 moves the slack
  of the generated row by `2 |x_i|`, which the same tolerance absorbs: no separate tolerance decides what is zero.
  208 generations had an exact zero in `x` and 698 a weight within the tolerance of 0. At the final `x` of
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
- **The iteration limit.** The family counts as n rows in the limit `10 (n + rows) + 10`, since at most n generated
  rows are active at once. The largest count measured was 49 iterations (0.15 of the limit) on the battery and 1627
  (0.27 of the limit) at n = 300, `cond(Sigma) = 1e8`, `L = 1.1`.
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
  `kappa = 1'y` at most 9.7e-17 of `sqrt(n) |y|max` (the scale of the row at the point); in the other 1039, `kappa` was
  at least 0.14 of it. The row alone does not decide the boundary, below.
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
    rule `kappa <= 0`, which it contains.
  - **Only on unbounded cones.** `kappa_noise` grows like `n^2 eps cond(Sigma) kappa`. For box-feasible points
    (`n = 150`), `16 kappa_noise / kappa` has a median of 2e-4 at `cond(Sigma) = 1e8`, 1.5e-2 at 1e10 and 2.0 at 1e12.
    Applied everywhere, the threshold returned `UNBOUNDED` on bounded feasible sets, where it is impossible:
    - at 1e12 and `n = 150`, 5 of 6 boxes [-1, 1];
    - at 1e13, some long-only and leveraged problems;
    - at 1e14, every bounded type.

    At `kappa = 0` the rows read `y_i >= 0` (finite lower), `y_i <= 0` (finite upper), `1'y = 0` and, with a leverage,
    `||y||_1 <= 0`. So the cone is bounded, and `kappa` stays away from 0, exactly when one of two things holds:
    - a leverage is set;
    - no pair `i != j` has `upper_i = +inf` and `lower_j = -inf` (`cone_bounded`).

    There only the row and `kappa <= 0` decide. On the reviewer's construction, with means of either sign, `n` 50 and
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
    - every `OPTIMAL` result has an exact `kappa` of at least 19 `kappa_noise`, the lower rows of the bounded assets
      active, and weights whose relative error bounds the error of `kappa` by 5.8e-3 `kappa_noise`;
    - the 58 `UNBOUNDED` results with an exact `kappa > 0` have at most 8.0 `kappa_noise`;
    - none above the threshold is `UNBOUNDED`.
  - **The boundary side.** The error of the computed `kappa` against its exact value, in units of `kappa_noise`, is at
    most 0.040: on `MaxSharpe.RiskFreeAtMinimumVarianceReturn` at `rf = r_gmv`, where the data make `kappa` exactly 0
    (16 problems, `cond(Sigma)` 1e2 to 1e8, both scales; at most 0.040), on
    `MaxDiversification.InverseCovarianceTimesVolatilitiesSumsToZero`, also exactly 0 (24 problems; at most 0.016), and
    on the reviewer's probes (`n = 6`, `cond(Sigma)` 1e2 to 1e8, both scales, 3 draws each, `rf = r_gmv` computed in
    floating point, at 0 and 1 ulp from it and moved by 3 to 3e4 times `cond(Sigma) eps` relative, 312 solves against an
    exact rational `kappa`; 187 with the row inactive, at most 0.029). The factor 16 is 400 times the largest error.
  - **The battery side.** The 1039 `OPTIMAL` solves of `RatioFuzz` have `kappa` at least 2.5e6 `kappa_noise`, 1.6e5
    times the threshold, and the factor changes no status of the battery. Outside the boundary tests, the unbounded
    cones of the unit suites of `MAX_SHARPE` and `MAX_DIVERSIFICATION` are at least 1.1e11.
  - **What it gives.** On the probes no result is `OPTIMAL` with an exact `kappa <= 0` (several were before). The 99
    `OPTIMAL` results have an exact `kappa` of at least 23 `kappa_noise`, weights from 11 to 1.3e12 and a relative error
    of the weights of at most 6.5e-4, consistent with an error of `kappa` of at most 0.029 `kappa_noise` over a `kappa`
    above 16: backward-stable results whose band is relative to weights far from order 1 (see "Tolerances"). The other
    side of the same honesty: 54 problems whose exact `kappa` is positive but at most 16 `kappa_noise` are
    `UNBOUNDED`, with exact weights of at least 3.2e6, since their sign is within the rounding of the data.
- **Stationarity above `cond(Sigma)` = 1e8: a conditioning limit.** Where many weights are interior, the relative
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
  budget, the bounds and the leverage at the 1e-7 band, and the whole KKT only for long-only at 1e12 (margin 74 at `n = 150`; 3.9e-8 at 1e13 would leave 2.6). No
  band is added for the stationarity there.
- **Cash bounds.** With cash the cash weight is held at exactly 0, and the request is `INFEASIBLE` when the cash bounds
  exclude 0: `lower_0 > 0` or `upper_0 < 0`, compared with 0 exactly.
- **The multipliers.** The mapping uses `nu = 2 sigma_y^2`, which complementarity gives in exact arithmetic, rather
  than the QP's `nu`, and computes `sigma_y = |L'y|` from the factor of `Data`. The difference between the two is the
  rounding of the active rows, and enters the user's stationarity as a term along `v`: the worst stationarity of the
  battery is at most 0.087 of its band (below).

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
  band. All 1520 cases pass, 794 of them with an inactive cap, with no case pinned; the worst agreement with the
  reference is 1.1e-12 at `cond <= 1e4`, 8.2e-12 at 1e6 and 8.1e-10 at 1e8.
- **`MaxReturnFuzz.NeverSilentlyWrong`**: the nearly-equal-mean shapes of `fuzz_shapes(400, 2)` (100) with a nearly
  riskless asset appended (variance log-uniform in [1e-10, 1e-6], mean 0.4 times the smallest), caps at 0.1, 0.5 and
  0.9 of the way to the largest asset volatility and at 1.01, 2 and 10 times `sigma_min`. The covariance solved has a
  condition number from about 4e2 to 4e16, so most cases take the 1e-7 band. The result must be `OPTIMAL` with the
  return of the reference, every weight within its bounds and `vol <= s (1 + tol)`, or `MAX_ITER`. Of the 528 cases,
  413 are `OPTIMAL` (worst agreement 7.5e-9, at a band of 1e-7) and 115 `MAX_ITER`.
- **`LeverageFuzz.AgainstSignPatterns`**: `fuzz_shapes(256, 4, 8)` (every combination once, n from 3 to 8) with
  `with_leverage`, for `MIN_VARIANCE`, `UTILITY`, `IMPORTANCE`, `KELLY` and `TRACKING_ERROR` (`MAX_RETURN` is checked by
  `MaxReturnFuzz.LeverageFuzz`), each at the shape's leverage and, where the optimum without leverage is above the least gross exposure the
  budget allows (1 without cash, 0 with cash), at a leverage a fraction in [0.2, 0.9] of the way down to it.
  `TRACKING_ERROR` runs over `T = n` periods `R = L' + rf` (`L` the factor of `Sigma`, so `R~'R~ = Sigma` and the
  band is that of the shape). `INFEASIBLE` exactly when no pattern is `Optimal`; otherwise `OPTIMAL` with the value of
  the reference (relative to the larger of the magnitude of its terms and `|g|max |w|max` over every weight, cash
  included, since the optimum at all cash is 0 to rounding) and `expect_user_kkt`, within the band. Of the 1979
  cases, 699 bind (61 to 179 per objective, 404 with two or more active generated rows, up to 8), 1056 do not and 224
  are `INFEASIBLE` (L < 1 without cash, or below what the cash cap of the tight box or the fixed weights force); the
  worst agreement is 3.1e-14 at `cond <= 1e4`, 6.3e-13 at 1e6 and 1.1e-10 at 1e8, at most 1.1e-3 of the band.
- **`MaxReturnFuzz.LeverageFuzz`**: `with_leverage(fuzz_shapes(400, 3))` without the fixed box (304 shapes), at the
  caps of `AgainstFrontierReference` from `sigma_min` under the leverage (the largest asset volatility alone where no
  portfolio is feasible under it), each at the shape's leverage and, where the `MAX_RETURN` optimum without leverage
  is above the least gross exposure the budget allows, at a leverage a fraction in [0.2, 0.9] of the way down to it.
  The reference is `max_return_reference` with the leverage in its `MIN_VARIANCE` (the frontier volatility still
  increases above the minimum-variance return, the leverage being convex); the checks are those of
  `AgainstFrontierReference`. Of the 2314 cases, 2168 are `OPTIMAL` (1129 with an inactive cap) and 146 `INFEASIBLE`;
  1057 bind the leverage, 763 of them with the cap. The worst agreement with the reference is 1.5e-12 at
  `cond <= 1e4`, 4.2e-12 at 1e6 and 5.2e-10 at 1e8. `gcc-debug` runs the first 256 shapes (every combination once) at
  the caps 0.1, 0.5 and 0.9: 407 cases bind the leverage, 278 with the cap.
- **`RatioFuzz.AgainstFrontierReference`**: `fuzz_shapes(300, 3)` for `MAX_SHARPE` and `MAX_DIVERSIFICATION`, at each
  shape without a leverage and with `with_leverage`. The status must be that of `max_ratio_reference`. An `OPTIMAL`
  result must have the reference's ratio, the cash weight exactly 0 and `expect_user_kkt` with the gradient written in
  the test, `-v / vol + (v'w) Sigma w / vol^3`, and cash held to the bounds [0, 0]; all within the band. A result
  without a point has every weight and bound multiplier, the budget and leverage multipliers and the residual NaN. All 1200 cases
  pass, with no case pinned. 1039 are `OPTIMAL` (91 bind the leverage), 133 `INFEASIBLE` (`rf` above every mean with a
  box, and a leverage below 1) and 28 `UNBOUNDED` (24 `MAX_SHARPE`, 4 `MAX_DIVERSIFICATION`, all without a box or a
  leverage). The worst agreement of the ratio is 1.4e-12 at `cond <= 1e4`, 2.0e-11 at 1e6 and 2.5e-9 at 1e8 (at most
  0.025 of the band). The worst residual of `expect_user_kkt`, the stationarity, is 8.2e-13 at `cond <= 1e4`, 8.7e-11
  at 1e6 and 8.4e-9 at 1e8 (at most 0.087 of the band). Every other measure stays below 1e-3 of its band.
- **The known residual.** One case of `MaxReturnFuzz.LeverageFuzz` misses its band and is pinned by name: shape seed
  3963545680699380072 (n 24, `cond(Sigma) = 1e8`, annual scale, nearly equal means, no box, cash, the shape's leverage
  2.154), cap fraction 0.7 (cap 858.4). It is an acceptance at an adjacent end by `kAdjacentEnd`, not an error of the
  face shift. The pin asserts `OPTIMAL`, the return of the reference, `expect_user_equalities`,
  `expect_user_bounds`, `expect_user_leverage` and `expect_user_stationarity` within the band, `vol <= s (1 + tol)`
  and `lambda_vol <= 0`, and relaxes only the volatility complementarity, to at most 1.1e-7; every other case keeps
  its full band, and any other miss fails. The cap is reached at `t = 6.7e15` on the leverage face, shifted, with 23 generated
  rows active; there `vol(Q(t))` varies by 1.5e-7 relative between `t` that differ in the last digits (858.428732 to
  858.428857), the rounding of `1e-15 cond(Sigma)`. The bisection never holds one active set across a bracket and ends
  at adjacent doubles, whose lower end is accepted (`kAdjacentEnd`) with `vol = s (1 - 1.006e-7)`: the volatility
  complementarity of `expect_user_kkt` is 1.006e-7 against the band of 1e-7. The return agrees with the reference to
  6e-14 relative and `lambda_vol = -1.3e-13`. The case runs in `gcc-release` only (the debug run has no cap fraction
  0.7). Task 0's pinned case of `AgainstFrontierReference` (shape seed
  11394049047870518593, cap fraction 0.9, cash, long-only, nearly equal means, `cond(Sigma) = 1e8` daily, `MAX_ITER`
  where `Q(t)` was lost near `t = 1.2e12`) is solved by the shift along the cash bound and is no longer pinned.
- **Runtime.** In `gcc-release` `AgainstFrontierReference` 1.3 s, `NeverSilentlyWrong` 0.8 s, `MaxReturnFuzz.LeverageFuzz`
  6.5 s, `LeverageFuzz.AgainstSignPatterns` 0.9 s and `RatioFuzz.AgainstFrontierReference` 3.1 s; in `gcc-debug` 25.0 s,
  11.9 s, 37.7 s (117 s in full, hence the smaller debug run above), 11.6 s and 43.7 s (in full). The other tests
  run their full counts in both presets.

**Plan-2 scripts against the fixed DLL.** The Python probes of the plan-2 review (`MAX_RETURN` against Clarabel, an
SOCP reference, and an independent KKT check) were run against the DLL before and after these changes. Their
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

## Where the tests are

| Suite | What it checks |
|---|---|
| `tests/unit/linalg/cholesky_test.cpp` | The factor and the triangular solves; rejection of an indefinite matrix and of a pivot at rounding level |
| `tests/unit/linalg/qr_test.cpp` | Householder QR: `Q'a = R`, the lower factor of `a'a`, least squares, rank deficiency (duplicated, combined and zero columns), and a rank threshold relative to the columns (`cond = 1e10` full rank, `1e15` not) |
| `tests/unit/qp/goldfarb_idnani_test.cpp` | Unconstrained, equality-only and bounded problems; infeasible and inconsistent ones; the scale tests above; a provided factor of G gives the same solution; feasibility at the final x from a start of 3e17; the gross exposure family against the sign-pattern brute force (100 problems, 85 with the family active), inactive (bit for bit the problem without it) and in the homogeneous form, offset `L 1` and right-hand side 0; a face of the family held as an equality at `t = 1e10` (the family does not regenerate it) |
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
| `tests/unit/problems/requirements_test.cpp` | The objective x setter matrix (leverage accepted by the eight objectives solved, `MAX_RETURN`, `MAX_SHARPE` and `MAX_DIVERSIFICATION` among them; `MAX_DIVERSIFICATION` needs no `mu`), its messages and the order of the checks |
| `tests/unit/problems/solve_test.cpp` | Argument checks of `po::solve` and their messages |
| `tests/unit/metrics/portfolio_test.cpp` | Portfolio return and volatility, utility and tracking error, with and without cash |
| `tests/unit/capi/` | The C boundary: status codes, `po_last_error()`, null handles, the new setters and `po_data_set_returns`, ABI constants, the metrics; `c_header_check.c` compiles the header as C11 |
| `python/tests/test_min_variance.py`, `test_qp_objectives.py` | The Python package against closed forms, `numpy` and the status contract; the frozen `Spec`, the objective identities, `MAX_RETURN` on the capital market line, the input checks of `Data` |
| `python/tests/test_ratio_objectives.py` | `MAX_SHARPE` and `MAX_DIVERSIFICATION` against the closed forms (`numpy`), the cash weight 0, `UNBOUNDED` and `INFEASIBLE`, the leverage and its multiplier, the hash of a `Spec` with a leverage and array bounds, `sharpe_ratio` and `diversification_ratio` against `numpy` |
