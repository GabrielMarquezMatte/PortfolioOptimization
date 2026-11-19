# Mathematics

The equations PortOpt solves, as the code writes them. This version covers `MIN_VARIANCE`, `UTILITY`,
`IMPORTANCE`, `KELLY`, `TRACKING_ERROR`, `MAX_RETURN`, `MAX_SHARPE`, `MAX_DIVERSIFICATION` and `RISK_PARITY`: every
objective is solved. The sources are listed at the end, with what was and was not verified.

## Notation

- N risky assets. `w` in R^N are the weights of the risky assets, `mu` in R^N the expected returns and `Sigma` the
  N x N covariance matrix, symmetric positive definite. `sigma_i = sqrt(Sigma_ii)`.
- **Cash (the risk-free asset).** When the data carries a risk-free rate `rf` (`po_data_set_risk_free`), there is a
  further asset N + 1 with return `rf` and zero variance and covariances. Its weight is `w_0 = 1 - 1'w`, and it is
  always the last entry of every weight vector. `Sigma` stays that of the risky assets only and is never
  regularized. The return of a portfolio is `rf + (mu - rf 1)'w` and its variance is `w' Sigma w`. Without cash the
  budget constraint `1'w = 1` holds.
- `mu~ = mu - rf 1` is the excess return, with `rf = 0` without cash.
- A weight vector has N entries without cash and N + 1 with cash.
- Units: every input is in the user's units (mu, Sigma, rf, the target return). The library converts nothing.

## MIN_VARIANCE

Without a box, a target return or cash:

    minimize    w' Sigma w
    subject to  1'w = 1

The solution is `w = Sigma^-1 1 / (1' Sigma^-1 1)`. With A = 1' Sigma^-1 1, B = 1' Sigma^-1 mu and
C = mu' Sigma^-1 mu, a target return `r` (the constraint `mu'w = r`) gives the frontier

    sigma^2(r) = (A r^2 - 2 B r + C) / (A C - B^2)

a parabola in `r`, minimal at `r = B / A`.

With cash, the budget disappears (the cash weight absorbs it) and the target return `r` becomes
`(mu - rf 1)'w = r - rf`. The frontier is then the capital market line

    sigma(r) = |r - rf| / sqrt((mu - rf 1)' Sigma^-1 (mu - rf 1))

and without a target return the minimum is all cash, `w = 0`.

With a box, each weight (cash included) has bounds `lower_k <= w_k <= upper_k`, either one pair for all weights
(`po_spec_set_bounds`) or one pair per weight (`po_spec_set_asset_bounds`, cash last); an infinite bound is no
bound. Infeasibility, such as an upper bound of 0.2 on three assets that must sum to 1, is a result with status
`PO_SOLVE_INFEASIBLE`, not an error.

## Objectives built on the covariance

`UTILITY`, `IMPORTANCE` and `KELLY` share `MIN_VARIANCE`'s constraints (budget or cash, the box and, for
`MIN_VARIANCE` only, a target return) and differ in the objective. Each is the minimization of `f` over the weights,
with `mu~ = mu - rf 1`:

| Objective | f | Setter |
|---|---|---|
| `MIN_VARIANCE` | `w' Sigma w` | `target_return` |
| `UTILITY` | `(delta / 2) w' Sigma w - mu~'w` | `risk_aversion` (delta > 0, default 5) |
| `IMPORTANCE` | `(1 - lambda) w' Sigma w - lambda mu~'w` | `importance` (0 <= lambda < 1, default 0.5) |
| `KELLY` | `w' Sigma w / (2 (1 + rf)) - mu~'w` | none; needs `1 + rf > 0` |

One builder (`src/problems/mean_variance.cpp`) serves all four. It writes the QP as `G = a Sigma`,
`g = -b mu~` and gives the solver the Hessian factor `sqrt(a) L`, with `L` the Cholesky factor of `Sigma` that
`Data` already holds, so no objective factors a scaled covariance again. `0.5 w'Gw + g'w` equals `f` and the
multipliers of the QP are the user's:

| Objective | a | b |
|---|---|---|
| `MIN_VARIANCE` | 2 | 0 |
| `UTILITY` | delta | 1 |
| `IMPORTANCE` | 2 (1 - lambda) | lambda |
| `KELLY` | 1 / (1 + rf) | 1 |

Two identities follow and are tested on the weights. `KELLY` is `UTILITY` with `delta = 1 / (1 + rf)` (delta = 1 without cash).
`IMPORTANCE` with `lambda` in (0, 1) is `UTILITY` with `delta = 2 (1 - lambda) / lambda`, scaled by `lambda`: the
same weights and multipliers `lambda` times those of `UTILITY`; at `lambda = 0` it is `MIN_VARIANCE`.

Without a box, a target return or cash, `UTILITY` has the closed form `w = Sigma^-1 (mu - nu 1) / delta` with
`nu = (B - delta) / A`, so that `lambda_budget = -nu = (delta - B) / A` (A and B as above). With cash and no box,
`w = Sigma^-1 mu~ / delta`.

## TRACKING_ERROR

With `R` the returns (T x N, one row per period, risky assets only) and `b` the benchmark (T entries), the
portfolio return of period t is `R_t w + rf w_0` with `w_0 = 1 - 1'w`. Writing `R~ = R - rf 1 1'` (without cash
`R~ = R`),

    f = |R~ w + rf 1 - b|^2

so the QP is

    G = 2 R~'R~,    g = -2 R~'(b - rf 1)

and `0.5 w'Gw + g'w` is `f` minus the constant `|b - rf 1|^2`. It needs the returns and not the covariance; a target
return also needs `mu`. `R~` is factored by Householder QR, `R~ = Q R_qr` (`src/linalg/qr.cpp`), and never through
`R~'R~`, whose condition number is the square of that of `R~`: since `R~'R~ = R_qr'R_qr`, the Hessian factor is
`sqrt(2) R_qr'` (rows signed so that its diagonal is nonnegative). That factorization is the one decision on the
column rank of `R~`: with fewer periods than assets, or a diagonal entry `|R_jj| <= c max(T, N) eps max_j |R~_j|`
(`c = 16`, `kRankFactor` in `src/linalg/qr.hpp`; `|R~_j|` the column norms), the request is rejected with
`PO_INVALID_ARGUMENT`, "returns minus the risk-free rate do not have full column rank numerically" (a duplicated or
constant column, or one that depends on the others to working precision). `po_tracking_error` reports the root mean
square `sqrt(f / T)`.

## Leverage: an implicit family of rows

The gross exposure limit

    sum_i |w_i| <= L

over the risky weights (cash is not counted) is set by `leverage` (`po_spec_set_leverage`, `L` positive and finite) and
applies to every objective but `RISK_PARITY`. Without cash the budget gives `sum_i |w_i| >= |1'w| = 1`, so `L < 1` is
`PO_SOLVE_INFEASIBLE`. With cash the cash weight absorbs the budget, so a limit below 1 can be feasible
for the objectives that keep the cash weight free (all cash is feasible when the bounds allow it); `MAX_SHARPE` and
`MAX_DIVERSIFICATION` fix it at 0, so for them `L < 1` is infeasible with cash too.

The constraint is not linear, but it is the intersection of `2^n` linear rows: `sum_i |w_i| <= L` holds exactly when
`s'w <= L` for every `s` in `{-1, 1}^n`. `qp::GrossExposure` states the family `a'x - ||x||_1 >= b`, that is the rows
`(a - s)'x >= b` over every `s`; the user's leverage is `a = 0`, `b = -L` and the homogeneous form of `MAX_SHARPE` is
`a = L 1`, `b = 0`. The solver never lists the rows. At a point `x`, `s'x` is largest at `s = sign(x)` (taking `+1` where
`x_i = 0`), so the row with that `s` is the most violated one and competes with the explicit rows like any of them; when
it wins, it is generated and added to the active set, and the active generated rows are numbered after the explicit ones.
At most `n` of them are active at once. The details and measurements are in
[validation.md](validation.md#the-gross-exposure-family).

Stationarity then reads `G x + g = E'y + C'z + sum_j zeta_j (a - s_j)`, with `zeta_j >= 0` the multipliers of the active
generated rows (`src/qp/problem.hpp`). For the user's leverage `a = 0`, and `-sum_j zeta_j s_j = lambda_lev g` with

    lambda_lev = -sum_j zeta_j <= 0,    g = sum_j zeta_j s_j / sum_j zeta_j

The family is generated lazily by the dual active-set method; the same construction (the constraint `s'w <= L` for every
sign vector `s`, generated as it is violated) is the one Tibshirani (1996, section 6) gives for the lasso, and with rows
that are only added it ends in finitely many steps. (The termination argument of Goldfarb-Idnani with removals was not
verified in the 1983 paper.)

`g` lies in the subdifferential of `||w||_1`: on an active row `s_j'w = ||w||_1`, so `s_ji = sign(w_i)` where `w_i` is not
0, and where `w_i = 0` the active rows may disagree, which gives any `g_i` in [-1, 1]. `lambda_lev` is the derivative of
the optimal `f` with respect to `L`; it is 0 when the leverage is set but no generated row is active, and NaN when no
leverage was set.

## MAX_RETURN

The path that follows is piecewise affine in the parameter, and between two breakpoints the optimum is the
interpolation of the two ends: that is Proposition 1 and Corollary 1 of Niedermayer and Niedermayer (2007) for the
critical line algorithm. The solver doubles and bisects on `t` instead of running their exact ratio test.

    minimize    -mu~'w
    subject to  sqrt(w' Sigma w) <= s

over the same budget, box and cash constraints as the other objectives, with `s` from `max_volatility`. It needs
`mu`, the covariance and `s`, and takes no target return.

The volatility constraint is not linear, so the QP solver cannot take it. `MAX_RETURN` instead solves the
parametric QP

    Q(t):  minimize  0.5 w' Sigma w - t v'w    over the budget, box and cash rows,    v = mu~ - c 1

(`G = Sigma` with the factor `L`, `g = -t v`; only `g` depends on `t`). Where the rows hold `1'w = 1` (no cash), `c`
is the midrange of `mu~`, and `-t v'w = -t mu~'w + t c` on every feasible `w`: the same optimum as with `mu~`. The
solver starts from the unconstrained minimizer `t Sigma^-1 v`, and the shift keeps out of it the common mean that the
budget absorbs; with nearly equal means it would otherwise swamp the weights at the large `t` where the cap is reached.
With cash, `c = 0`. The multipliers of `Q(t)` are those of `mu~` except the budget's, which carries `t c` more, so the
reported budget multiplier is the QP's over `t` (or the slope's) minus `c`. Below, `v` is `mu~` whenever `c = 0`.
Its optimum `w(t)` is piecewise affine in
`t` and its volatility is nondecreasing, from the minimum volatility `sigma_min` at `t = 0`. The solver finds the
`t*` where the volatility reaches `s`:

1. At `t = 0` it solves the minimum-variance point. `s` below `sigma_min` is `PO_SOLVE_INFEASIBLE`. Feasibility is
   decided there, once: a `Q(t)` that fails at `t > 0` has broken down along the path and is `PO_SOLVE_MAX_ITER`,
   never `PO_SOLVE_INFEASIBLE`.
2. It doubles `t`, from `t_1 = sigma_min^2 / |v|max` (with `sigma_min = 0`, all cash, from
   `min(s^2, max_i Sigma_ii) / |v|max`, so that a large cap does not start far past the last breakpoint; `|v|max` is
   floored at `eps`), until the volatility reaches `s`, or until the limit certificate below holds (the inactive
   cap).
3. It bisects the bracket. When both ends are on the same face (the same active explicit rows and, with the
   leverage active, the same support with the same signs and the same zeros; not the same sign vectors of the
   generated rows, of which the solver may pick several sets for one face), it computes the crossing on the face of
   the lower end:
   - the weights off the face's free coordinates are fixed exactly: 0 where the active sign vectors disagree, and the
     bound of each active bound row;
   - the tangent `dw = Z H_r^-1 Z'v` of the face, by the null space `Z` of the remaining active rows over the free
     coordinates (Householder QR) and the QR of `L'Z` for the reduced Hessian `H_r = Z'Sigma Z`, in the frame of the
     point's shift `v`, anchored at a solve of the face at `t_a` (`w_a`), never at the QP solution itself, whose
     rounding carries `eps |t Sigma^-1 v|`;
   - when no weight is fixed (the budget alone, for example) and the `m + 4` triangular solves with `L` of the range
     space cost fewer flops than the null space (`(m + 4) n^2` against `4 n^2 m + 2 n (n - m)^2`), the same face solve
     in the range space instead: in `y = L'w`, with `L^-1 N = Q R`, `y = Q [R'^-1 b; (Q' t L^-1 v)_2]` and
     `dy = Q [0; (Q' L^-1 v)_2]`. Its weights carry `eps cond(Sigma)` rather than `eps cond(Z'Sigma Z)`, and its
     stationarity residual `eps |L| |t L^-1 v|`, so a crossing it does not certify is computed again by the null
     space;
   - the root `h` of `|L'(w_a + h dw)|^2 = s^2`, written without cancellation, `t* = t_a + h`, `w* = w_a + h dw`, and
     the multipliers `u` by least squares in `N u = Sigma w* - t* v`.

   The crossing is accepted when the KKT conditions of `Q(t*)` hold at that one point (`Q(t*)` is convex, so they
   suffice): every row and the leverage hold to the QP solver's feasibility tolerance, each free weight keeps the sign
   of the face, the active inequality multipliers are nonnegative, and the leverage term is `zeta g` with `zeta >= 0`
   and `g` in the subdifferential of `||w||_1` (`g_i = sign(w_i)` on the support, `|g_i| <= 1` on the zeros), a form
   that does not depend on the sign vectors chosen. The dual tolerance is `16 n eps (t* |v|max + |Sigma|max |w|max)`.
   A crossing that does not certify (the face ends before `vol = s`) leaves the bisection to go on. A bracket whose
   ends are adjacent doubles takes the certified crossing of the face of either end, the lower first; when neither
   certifies the path has broken down (`PO_SOLVE_MAX_ITER`).

The reported point is the crossing: `vol(w*) = s` to the rounding of its evaluation, and `(w*, u)` solves `Q(t*)` to
the tolerances of the certificate.

At `t*`, `Sigma w - t* v = N u`, where `N u` collects the QP multipliers `u` of the rows. Dividing by `t*`:

    -v = N (u / t*) - (1 / t*) Sigma w = N (u / t*) + lambda_vol Sigma w / vol(w),    lambda_vol = -s / t*

since `vol(w) = s`. The multipliers of the budget and bound rows are the QP's divided by `t*`, and the volatility
multiplier is `lambda_vol = -s / t*`, the derivative of the optimal `f = -mu~'w` with respect to `s` (<= 0).

**When `w(t)` is constant.** Where `w` is constant, `Sigma w - t v = N u(t)` for every `t`, so
`u(t) = u(t_b) + (t - t_b) b` with `N b = -v`. The slope `b` is computed from that equation by Householder least
squares over the active normals, not by differencing two QP solutions (whose rounding grows with `t`); a slope below
`-16 n eps |v|max` is a decrease, and an active multiplier that decreases reaches zero at the vertex exit
`t_b + u_k(t_b) / |b_k|`. A combination that comes out negative by rounding (within those bounds) is set to 0
(`16` is `kNoiseFactor` in `src/problems/parametric_path.cpp`).

The **limit certificate** at a point of the path: every row holds to the QP solver's feasibility tolerance, no active
slope decreases, and `|v + N b|max <= 16 n eps` times the largest of its terms (`limit_holds`). These are the KKT
conditions of the linear program `max mu~'w` over the rows, with `b` as its multipliers, so the point is that
program's optimum; and `w` with `u(t) + (t' - t) b` satisfies the KKT conditions of `Q(t')` for every `t' >= t`.

At `s = sigma_min`, two points of the path on the same active set have the same weights when the equality rows and the
active rows fix `w` (a vertex, decided combinatorially: the solver keeps the active normals independent, so their rank
is their count, and a count of `n` is a vertex), or else when the weights differ by no more than the noise of the
stationarity equation `Sigma w - t v = N u`,

    16 n eps (max(1, |w|max) + t |v|max).

Special cases:

- **Inactive cap.** When the limit certificate holds at a point of the climb whose volatility is below `s`, `w` stays
  optimal for every larger `t` and the volatility never reaches `s`. The weights are the maximum-return point of the
  box and `lambda_vol = 0`; the other multipliers are the slopes `b` (the budget's minus `c`).
- **`s = sigma_min`** (to `1e-12` relative). Let `t_b = sup{t : w(t) = w(0)}`; `t` is halved from `t_1` until `w(t)` is
  constant with `w(0)`. If `t_b` is infinite (no decrease), `lambda_vol = 0` when the limit certificate holds there,
  and `PO_SOLVE_MAX_ITER` when it does not. If `t_b` is finite (a vertex exit),
  `lambda_vol = -s / t_b`. The efficient frontier is vertical at `sigma_min` (`t_b = 0`) only when the weights moved
  beyond their rounding at a larger `t` on the active set of `w(0)` and no exit comes before that `t`: then
  `lambda_vol = -inf` and the budget multiplier (when present), the bound multipliers and the leverage multiplier (when a
  leverage was set, active or not) are NaN, although those constraints exist. This and RISK_PARITY at
  `s = sigma_min` are the two cases where a multiplier is NaN for an existing constraint. When the weights and the multipliers disagree, the multipliers decide (`0` or
  `-s / t_b`), never `-inf`.

**The leverage along the path.** With a leverage the rows of `Q(t)` include the family above. Its generated rows
`(a - s_j)'w >= b` enter the active set like explicit rows: their normals are columns of `N`, their sign vectors `s_j` are
part of the identity of the active set, and the slope `b`, the limit certificate and the segment decisions treat them as
any other active row. When a weight crosses 0 the row generated at `sign(w)` changes, so the active set changes there,
which is a breakpoint of the path. The user's multipliers are those of `Q(t*)` over `t*` as before, with
`lambda_lev = -sum_j zeta_j / t*`; for an inactive cap, like the others, they come from the slopes `b`.

**The shift along a face.** The shift along the budget generalizes. A *face* is a row `r'w >= rho` (or `= rho`) whose
normal `r` has every entry -1 or 1, so that it holds the whole portfolio: the budget, a bound of the cash weight, or a
generated row of the leverage (`r = -s`). With the face held as an equality, every feasible `w` has `r'w = rho`, so with
`v = mu~ - c r`

    -t v'w = -t mu~'w + t c rho

and the optimum on the face is that of `mu~`; `c` is the midrange of the `r_i mu~_i`, which minimizes `|v|max`, and it may
be negative (for the leverage face `r = -s`). The multiplier of the face in `Q(t)` is that of `mu~` less `t c`. Where the
rows hold the budget (no cash) the budget is the face of every point. With cash there is no budget, and the shift is
along the inequality face that holds the portfolio once it binds, the cash bound or the leverage; the point of the path
is shifted only when the optimum restricted to that face solves `Q(t)`, that is when the multiplier of the face in `Q(t)`
is nonnegative, and `Q(t)` is solved unshifted otherwise. At the vertex `L e_best` of the leverage with cash,
`lambda_lev = -mu~_best`. The decision and its measurements are in
[validation.md](validation.md#numerical-decisions-of-max_return).

Without a box, `MAX_RETURN` lies on the efficient frontier. With cash it is the capital market line,
`w = s Sigma^-1 mu~ / sqrt(mu~' Sigma^-1 mu~)` over the risky assets, and `lambda_vol = -sqrt(mu~' Sigma^-1 mu~)`,
the negative of the slope of the line.

## MAX_SHARPE and MAX_DIVERSIFICATION

    minimize    f(w) = -v'w / sqrt(w' Sigma w)
    subject to  1'w = 1,  lower_i <= w_i <= upper_i,  sum_i |w_i| <= L (with a leverage)

with `v = mu~` for `MAX_SHARPE` (the Sharpe ratio `(return - rf) / vol` with `rf = 0` without cash) and
`v = sigma = sqrt(diag Sigma)` for `MAX_DIVERSIFICATION` (the diversification ratio, which needs no `mu`). Both ratios
are of risky assets only: with cash the cash weight is fixed at 0, so `1'w = 1` holds over the risky weights, and cash
bounds that exclude 0 are `PO_SOLVE_INFEASIBLE` before anything is solved. The bounds are those of the risky weights.

**Homogenization** (Cornuejols and Tütüncü; Choueifaty for the diversification ratio). The ratio is homogeneous of
degree 0 in `w`, so with `y = kappa w`, `kappa = 1'y > 0`, and the scale fixed by `v'y = 1`, maximizing the ratio is
minimizing `y' Sigma y`. Each constraint on `w` becomes homogeneous in `y`, and `kappa > 0` is relaxed to `1'y >= 0`:

    minimize    y' Sigma y
    subject to  v'y = 1
                y_i - lower_i 1'y >= 0      (finite lower bounds)
                upper_i 1'y - y_i >= 0      (finite upper bounds)
                1'y >= 0
                L 1'y - ||y||_1 >= 0        (with a leverage: the gross exposure family with offset L 1, rhs 0)

`src/problems/homogeneous.cpp` builds it with `G = 2 Sigma` (factor `sqrt(2) L`), `g = 0`, one equality row and the
rows in this order: the bound rows of each risky asset (lower then upper), then `1'y >= 0`.

- The QP `Infeasible`: no feasible `w` has `v'w > 0` (`v'y = 1` cannot be met). `PO_SOLVE_INFEASIBLE`.
- The row `1'y >= 0` in the active set of the optimum: the optimum is at `kappa = 0`, a limit of feasible portfolios
  whose weights grow without bound. The supremum of the ratio exists but is not attained: `PO_SOLVE_UNBOUNDED`, with
  NaN weights and multipliers. It needs an unbounded feasible set. At `kappa = 0` the rows read `y_i >= 0` (finite
  lower), `y_i <= 0` (finite upper), `1'y = 0` and, with a leverage, `||y||_1 <= 0`. So `kappa = 0` forces `y = 0`,
  which contradicts `v'y = 1`, exactly when a leverage is set or no pair `i != j` has `upper_i = +inf` and
  `lower_j = -inf` (long-only without a cap included). Without bounds it happens when `1' Sigma^-1 v <= 0`.
- On an unbounded feasible set only: at `1' Sigma^-1 v` of order rounding the sign of `kappa` is not decided by the
  data, so a `kappa = 1'y` at most 16 times its rounding (the row inactive) is `PO_SOLVE_UNBOUNDED` too, even where the
  exact `1' Sigma^-1 v` is positive and tiny and the exact weights are of the order of the inverse of that rounding;
  see [the boundary](validation.md#numerical-decisions-of-max_sharpe-and-max_diversification) in `validation.md`. On a
  bounded one `kappa` stays away from 0 and only `kappa <= 0` counts.
- Otherwise `PO_SOLVE_OPTIMAL` with `w = y / kappa` (and cash exactly 0).

**Why a leverage excludes `PO_SOLVE_UNBOUNDED`.** The generated rows give `||y||_1 <= L kappa`, and `v'y = 1` gives
`1 = v'y <= ||v||_inf ||y||_1 <= ||v||_inf L kappa`, so `kappa >= 1 / (L ||v||_inf) > 0`: the row `1'y >= 0` is never
active (this is derived here, not taken from the sources).

Without bounds or a leverage, the optimum is `y = Sigma^-1 v / (v' Sigma^-1 v)`, so `w = Sigma^-1 v / (1' Sigma^-1 v)`
when `1' Sigma^-1 v > 0`: the tangency portfolio for `MAX_SHARPE`, the most diversified one for `MAX_DIVERSIFICATION`.

**Multipliers in the user's `w`.** With `sigma_y = sqrt(y' Sigma y)`, `w = y / kappa` gives `sigma_w = sigma_y / kappa`
and `v'w = 1 / kappa`, so

    grad f(w) = -v / sigma_w + (v'w) Sigma w / sigma_w^3 = -kappa v / sigma_y + kappa Sigma y / sigma_y^3.

The stationarity of the QP (`src/qp/problem.hpp`) is

    2 Sigma y = nu v + C'z + eta 1 + sum_j zeta_j (L 1 - s_j)

with `nu` the multiplier of `v'y = 1`, `z` those of the bound rows, `eta` that of `1'y >= 0` and `zeta_j` those of the
active generated rows. Every inequality row has right-hand side 0, so multiplying by `y` and using complementarity and
`v'y = 1` leaves `nu = 2 sigma_y^2`. Substituting `Sigma y` into the gradient, the `v` terms cancel:

    grad f(w) = c (C'z + eta 1 + sum_j zeta_j (L 1 - s_j)),    c = kappa / (2 sigma_y^3) > 0.

The normal of a lower row is `e_i - lower_i 1` and that of an upper row `upper_j 1 - e_j`, so in the convention below

- `beta_i += c z_i` for an active lower row of weight i, and `beta_j -= c z_j` for an upper row;
- `lambda_lev = -c sum_j zeta_j`, since `-c sum_j zeta_j s_j = lambda_lev g` with `g = sum_j zeta_j s_j / sum_j zeta_j`
  in the subdifferential of `||w||_1`;
- the budget term is `lambda = c (eta - sum_lower z_i lower_i + sum_upper z_j upper_j + L sum_j zeta_j)`, where `eta` is
  0 at every `OPTIMAL` result, since an active `1'y >= 0` is `UNBOUNDED`. Without cash `lambda` is the budget
  multiplier. With cash the budget multiplier is NaN and the term is `beta_0 = -lambda`, the multiplier of the cash
  weight held at 0, which enters the gradient as `-beta_0 1`.

`f` is homogeneous of degree 0, so `w' grad f = 0`, and the budget term balances the others:
`lambda = -(sum_i beta_i w_i + lambda_lev ||w||_1)`. At an unconstrained tangency it is 0.

## RISK_PARITY

    minimize    f(w) = sum_i (c_i - b_i)^2,    c_i = w_i (Sigma w)_i / w' Sigma w
    subject to  1'w = 1,  mu'w = r (with a target return),  sqrt(w' Sigma w) = s (with a target volatility),
                lower_i <= w_i <= upper_i

`c` is the vector of risk contributions (it sums to 1 when `w' Sigma w > 0`) and `b` is the risk budget
(`po_spec_set_risk_budget`, default `1 / N`; every `b_i` positive and finite, `|sum_i b_i - 1| <= N eps`, the rounding of
a sum of `N` terms). `RISK_PARITY` has no cash and no leverage: `risk_free` in the data is `PO_INVALID_ARGUMENT`
(`"risk_free does not apply to RISK_PARITY"`), and so is `leverage` (`"leverage does not apply to RISK_PARITY"`), and
the same for `max_volatility`, `risk_aversion` and `importance`. `target_return` needs `mu`. `target_volatility` is for
`RISK_PARITY` only (`"target_volatility does not apply to <OBJECTIVE>"` elsewhere). A budget with `k` entries on data of
`N` assets is `"risk_budget has k entries; expected N"`.

**The derivatives.** With `p = Sigma w` and `v = w' Sigma w`, `c = w o p / v` and

    J = dc/dw = (diag(p) + diag(w) Sigma - 2 c p') / v,        grad f = 2 J'(c - b)

(the first two terms differentiate the numerator `w_i p_i`, the last the denominator, `dv/dw = 2 p`). `c` is homogeneous
of degree 0, so `J w = 0` and `w' grad f = 0`, and `J'J` is singular along `w`. `f` is not convex, and the exact Hessian
is indefinite in general (Bai, Scheinberg and Tütüncü, 2016, examples 3.1 and 3.2); it is used only in the final Newton
step on the KKT system.

**The local character.** The positive-orthant solution of `c = b` (`f = 0`) exists and is unique for `Sigma` positive
definite (Bruder and Roncalli, 2012, section 2.3.3). That is not the only one in general: for each sign pattern
`beta` in `{-1, 1}^N` the problem `min 1/2 x' Sigma x - sum_i b_i ln(beta_i x_i)` over `beta_i x_i > 0` has a unique
solution `x^beta`, with `x^(-beta) = -x^beta`, and the normalization `x^beta / 1'x^beta` exists exactly when
`1'x^beta != 0`. So there are at most `2^(N-1)` normalized solutions, one per pair of opposite orthants, and none for a
pair whose weight sum is 0 (Bai, Scheinberg and Tütüncü, 2016, lemma 2.2). With bounds that admit negative weights and
the positive solution outside the box, an exact solution in another orthant can lie inside it, while the iteration
below, started from the positive one, ends at a KKT point with `f > 0`. `PO_SOLVE_OPTIMAL` therefore means a KKT point
reached from the Spinu solution, not the global minimum (`RiskParityFuzz.AgainstScanReference` finds 6 local minima
among 91 shapes of `N = 3`).

**Spinu (the fast path).** The minimizer `y` of `F(y) = 1/2 y' Sigma y - b' log y` over `y > 0` has `y_i (Sigma y)_i =
b_i`, so `w = y / 1'y` has `c = b`. `F` is strictly convex and coercive, `grad F = Sigma y - b / y` and
`hess F = Sigma + diag(b / y^2)`. Newton with Armijo backtracking that keeps `y > 0` by a fraction of the way to the
boundary. `F` is not standard self-concordant when `b_i < 1` (the third derivative of `-b ln y` is `2 b / y^3` against
`(b / y^2)^(3/2)`), but `F / min_i b_i` is, which gives a damped step whose iteration bound scales with `1 / min_i b_i`
(a derivation of this work); the solver uses Armijo instead. It stops on the size of the step just taken,
`max_i |dy_i / y_i| <= 16 eps` or when that size stops halving at `<= 1e-8`, because the Newton decrement bounds
`F - F*` and not the error of the weights (measurements in [validation.md](validation.md#numerical-decisions-of-risk_parity)).
Without a target return, a target volatility or an active bound, that point is the answer: every multiplier is 0 and `f`
is at its rounding floor. The same holds when the box contains it.

**The polyhedron and the descent.** With a target or a bound that cuts the Spinu point off, the feasible set of the
linear rows is the polyhedron `Omega = {1'w = 1, mu'w = r, lower <= w <= upper}`. From the point of `Omega` nearest to the
Spinu point (one QP) the solver runs a Gauss-Newton sequence of QPs on `Omega`: with `J` and `c - b` at `w_k`, the step
`d` minimizes

    1/2 d' (2 J'J + tau I) d + (2 J'(c - b))' d       over   w_k + d in Omega

with the existing Goldfarb-Idnani solver, followed by an Armijo line search along `d`. The proximal term `tau I` is there
because `J'J` is singular. It is Levenberg-Marquardt: `tau_0 = 0.05 tr(2 J'J) / N`, divided by 3 after a full step and
multiplied by 4 after a shortened one. This is scale free (`J` does not change when `Sigma` is scaled); the formula
`0.05 tr(Sigma) / (2N)` that the first plan had is not, and it was replaced after measurement
([validation.md](validation.md#the-polyhedral-descent)). The idea of a convex-QP sequence on the linear constraints is that
of SCRIP (Feng and Palomar, 2015) and of Bai, Scheinberg and Tütüncü.

**The volatility as an equality.** `h(w) = sqrt(w' Sigma w) - s = 0` is not linear and no source treats it, so it goes
in an augmented Lagrangian over `Omega` (Andreani, Birgin, Martinez and Schuverdt, 2007, algorithm 3.1):
`L(w) = f(w) - lambda h(w) + (rho / 2) h(w)^2`, each outer iteration a descent on `L` over `Omega`, then
`lambda <- lambda - rho h`, `rho` times 10 when `|h|` does not fall by half, a safeguarded `lambda` and a decreasing inner
tolerance. `grad h = Sigma w / vol(w)`.

**The polish and the gate.** The point is polished by Newton on the KKT system of its active set (the exact Hessian of the
Lagrangian, the active rows as equalities), and a final test decides the status: the KKT residual, feasibility included,
must be within 16 times the rounding noise of the problem, a bound derived from the data and not a fixed number (the
stopping factor inside the solver is 2). The noise is that of the compensated dot product (Dot2, Ogita, Rump and Oishi,
2005) used for `Sigma w`.

**Status.**

- `PO_SOLVE_OPTIMAL`: the gate passed.
- `PO_SOLVE_INFEASIBLE`, only when certified: the linear rows are incompatible (the minimum-variance QP over them is
  infeasible), `target_return` is outside the achievable range, or `s < sigma_min (1 - 1e-12)`.
- `s = sigma_min` within `1e-12` relative (`qp::kRelativeTolerance`, the feasibility tolerance of the QP): the only feasible
  point is the minimum-variance portfolio, and the result is `PO_SOLVE_OPTIMAL` there with `lambda_vol = -inf` and every
  other multiplier NaN. The gradients of the constraints are dependent, the LICQ fails, and in general no KKT multipliers
  exist (`df*/ds = -inf`); the same convention as `MAX_RETURN` on a vertical frontier. A target return that fixes `w` in
  two assets is another case: there the LICQ holds and the multipliers exist.
- `PO_SOLVE_MAX_ITER` in every other case, and it carries **no point**: weights, multipliers and `po_result_residual` are
  NaN (a rule that is now uniform for every objective). Its causes are a target volatility the iteration cannot reach (above
  the volatility of every vertex it reaches; the largest volatility of a polyhedron has no cheap certificate), a target
  return with no bounds (`f` is homogeneous of degree 0 and has no minimizer along the line the two equalities leave), a
  solution whose gross exposure makes the budget row unresolvable in double (`eps (sum_i |w_i| + 1) > 1e-12`, a gross
  exposure above about 4500, which is `MAX_ITER` even when it is stationary), and a bounded budget-only problem the
  solver does not solve (Capped box, equal budget, `N = 20`, `cond(Sigma) = 1e8`, seed 31696: documented, not fixed).

**Multipliers.** The stationarity condition is that of the [convention below](#multiplier-convention),

    grad f = lambda_budget 1 + lambda_return mu + lambda_vol Sigma w / vol(w) + sum_i beta_i e_i

with `beta_i >= 0` at an active lower bound and `<= 0` at an active upper bound. The budget multiplier is finite (no
cash); `lambda_return` and `lambda_vol` are NaN when their setter was not called; `lambda_lev` is always NaN; `beta_i` is 0
for a weight with no finite bound. Each is the derivative of the optimal `f` with respect to its right-hand side, so
`lambda_vol` is `df*/ds`, of either sign. Because `w' grad f = 0`, multiplying the condition by `w` gives the identity

    lambda_budget = -lambda_return r - lambda_vol s - sum_i beta_i w_i

(`w' Sigma w / vol(w) = s`), which holds to rounding at every `OPTIMAL` result and is checked by
`RiskParity.HomogeneityIdentity`, and `RiskParity.MultipliersAreDerivatives` checks them against central differences. At the Spinu point all terms are 0.

**What the weights are good to.** With a target volatility and few or no active bounds at `cond(Sigma)` around 1e8, the
weights are determined only to the rounding of the KKT residual amplified by the conditioning of the problem: an `OPTIMAL`
result can be up to about 1e-4 relative from the exact KKT point (measured up to 3.7e-5 over 229 results), while its
objective value is right to its band. **The multipliers** are not settled by least squares at the rounded weights
alone, which loses the `H dw` term of the Newton step (relative errors of 3.7e-2 at `cond(Sigma) = 1e8` were measured
that way): one more Newton solve at the returned weights corrects them, `lambda + dlambda` for the equalities and
`beta_j = (grad f - A'lambda)_j + (H dw)_j` for a held weight, kept only when a safeguard accepts it (a finite
solution, a second solve that contracts, a stationarity within the noise of the gate and the signs of the held
bounds; not the rank test of the QR, which the scaled KKT matrix fails at `cond(Sigma) = 1e8` and `N >= 12` on about
half the shapes while the correction is right). With a target volatility at least `1e-3` relative above `sigma_min`
their measured relative error is then up to 2.7e-13 at `cond(Sigma) = 1e4`, 5.2e-11 at 1e6 and 5.3e-10 at 1e8 (bands
1e-10, 1e-9, 1e-7; `N` from 5 to 50, bound multipliers on `N` = 5 and 8), at the weights those solves returned;
measured, not derived. Where the weights are further off the multipliers follow (8.7e-5 at a weight error of 2.3e-5).
Nearer the minimum, `s = sigma_min (1 + delta)`, the error grows roughly like `1 / delta`: at `cond(Sigma) = 1e4` up to
2.7e-12 at `delta = 1e-4`, 5.9e-11 at `1e-6`, 7.8e-9 at `1e-8` and 4.5e-7 at `1e-10` (1.3e-6 at `1e-10` over
`cond(Sigma) <= 1e8`). At 1e10 an `OPTIMAL` result is not validated: the multipliers were measured up to 0.25 off,
with weights up to 1.5e-2 off. `RiskParity.MultipliersAreDerivatives` runs on
well-conditioned data only. A solve with a target volatility, no box and `N = 50` at `cond(Sigma) = 1e8` can take
about 10 s; at `N <= 20` it takes at most about 0.15 s. See
[validation.md](validation.md#numerical-decisions-of-risk_parity).

## Which setters and data each objective takes

Any other combination is `PO_INVALID_ARGUMENT`:

| Objective | requires | also accepts, besides the bounds | `leverage` |
|---|---|---|---|
| `MIN_VARIANCE` | covariance (mu with a target) | `target_return` | yes |
| `MAX_RETURN` | mu, covariance, `max_volatility` | none | yes |
| `UTILITY` | mu, covariance | `risk_aversion` | yes |
| `IMPORTANCE` | mu, covariance | `importance` | yes |
| `KELLY` | mu, covariance, `1 + rf > 0` | none | yes |
| `TRACKING_ERROR` | returns (mu with a target) | `target_return` | yes |
| `MAX_SHARPE` | mu, covariance | none | yes |
| `MAX_DIVERSIFICATION` | covariance | none | yes |
| `RISK_PARITY` | covariance (mu with a target) | `target_return`, `target_volatility`, `risk_budget` | no |

`leverage` is accepted by every objective except `RISK_PARITY` (`"leverage does not apply to RISK_PARITY"`), which also
takes no cash (`"risk_free does not apply to RISK_PARITY"`). Every objective is solved, so `PO_UNSUPPORTED` no longer
occurs for a valid objective. `po_solve` checks, in this order: every set setter applies to the objective, the data it requires is present, and per-asset bounds match the weight count. The
messages are `"<setter> does not apply to <OBJECTIVE>"`, `"<OBJECTIVE> requires <item>"` (the items are `mu`,
`covariance`, `returns` and `max_volatility`), `"KELLY requires risk_free > -1"` and `"risk_budget has k entries;
expected n"`. The setters are the trust
boundary of the values: `risk_aversion` is positive and finite, `importance` is in [0, 1), `max_volatility` is
positive and finite, `target_return` is finite, `leverage` is positive and finite, `target_volatility` is positive and
finite and `risk_budget` has positive finite entries that sum to 1 within `n eps`.

## The QP rows

`qp::Problem` is `min 0.5 x'Gx + g'x` subject to `E x = e` and `C x >= c`. Every objective sets the variables to
the risky weights, `x = w`. MIN_VARIANCE has

    G = 2 Sigma,    g = 0

so that `0.5 w'Gw = w' Sigma w` exactly and the multipliers of the QP are the user's (see below); the `G` and `g` of
the other objectives are given above, and the rows are the same for all of them. The rows come from
`src/problems/constraints.cpp`, in this order; a row exists only if its constraint does.

| Constraint | Row | Kind | Right-hand side |
|---|---|---|---|
| Budget, no cash | `1'w` | equality | `1` |
| Target return `r`, no cash | `mu'w` | equality | `r` |
| Target return `r`, with cash | `(mu - rf 1)'w` | equality | `r - rf` |
| Lower bound `l_i`, risky asset i (finite) | `e_i'w` | inequality `>=` | `l_i` |
| Upper bound `u_i`, risky asset i (finite) | `-e_i'w` | inequality `>=` | `-u_i` |
| Lower bound `l_0` of cash (finite) | `-1'w` | inequality `>=` | `l_0 - 1` |
| Upper bound `u_0` of cash (finite) | `1'w` | inequality `>=` | `1 - u_0` |

The cash rows follow from `w_0 = 1 - 1'w`: `w_0 >= l_0` is `-1'w >= l_0 - 1`, and `w_0 <= u_0` is `1'w >= 1 - u_0`.
Bounds of the risky assets come first, one lower then one upper per asset, then those of cash. The leverage is not a
row of this table: it is the implicit family of [the leverage section](#leverage-an-implicit-family-of-rows), generated
as the solver needs it.

## Multiplier convention

The solver uses `G x + g = E'y + C'z` with `z >= 0` (`src/qp/problem.hpp`). Writing `f(w)` for the objective in the minimization forms above (`w' Sigma w` for MIN_VARIANCE, without the
factor 1/2; a maximized objective enters negated), the stationarity condition at the optimum, over the risky weights,
is

    grad f(w) = lambda_budget 1 + lambda_return grad(return)(w) + lambda_vol Sigma w / vol(w)
                + lambda_lev g + sum_i beta_i e_i - beta_0 1

- `lambda_budget`, `lambda_return`, `lambda_vol` and `lambda_lev` are `po_result_multiplier` for
  `PO_CONSTRAINT_BUDGET`, `PO_CONSTRAINT_TARGET_RETURN`, `PO_CONSTRAINT_VOLATILITY` and `PO_CONSTRAINT_LEVERAGE`.
  `grad(return)` is `mu`, or `mu - rf 1` with cash, `Sigma w / vol(w)` is the gradient of `vol(w) = sqrt(w' Sigma w)`
  and `g` is in the subdifferential of `||w||_1` over the risky weights: `g_i = sign(w_i)` where `w_i` is not 0 and any
  value in [-1, 1] where `w_i = 0`.
- The term of an absent constraint is omitted: there is no budget term with cash, no return term without a target
  return, no volatility term unless the objective is MAX_RETURN or RISK_PARITY with a target volatility, no leverage term without a leverage and no `beta_0`
  without cash. With cash, `MAX_SHARPE` and `MAX_DIVERSIFICATION` fix the cash weight at 0, so their budget term is
  carried by `beta_0` and the budget multiplier is NaN.
- `beta_k` is the lower-bound multiplier minus the upper-bound multiplier of weight k (`po_result_bound_multipliers`):
  `>= 0` at an active lower bound and `<= 0` at an active upper bound. It is 0, not NaN, for a weight without a
  finite bound.
- `beta_0` is the entry of cash, the last one. It enters with a minus sign because the cash weight is `1 - 1'w`: the
  cash rows have the opposite signs of those of an asset (see the table).
- A multiplier of an absent constraint is NaN: the budget when there is cash, the target return when none was set
  (and for the objectives that take none), the volatility constraint unless the objective is MAX_RETURN or RISK_PARITY with a target volatility, and the
  leverage constraint when no leverage was set (always, for RISK_PARITY). There are two documented exceptions, where
  existing constraints have NaN multipliers: the vertical `s = sigma_min` case of MAX_RETURN above (the budget
  multiplier when present, the bound multipliers and the leverage multiplier when a leverage was set), and RISK_PARITY
  at `s = sigma_min` (the budget and bound multipliers; see its section).
- Each multiplier is the derivative of the optimal `f` with respect to the right-hand side of its constraint. A `<=`
  constraint therefore has a multiplier `<= 0`, like `beta` at an upper bound, and a `>=` constraint one `>= 0`. For
  MAX_RETURN, `lambda_vol <= 0` is the derivative with respect to the cap `s`, and 0 when the cap is inactive. For
  RISK_PARITY it is the derivative of the optimal `f` with respect to the target `s`, of either sign.
  `lambda_lev <= 0` is the derivative with respect to the limit `L`, and 0 when the leverage is set but inactive.
- When the constraints are dependent (a target-return row parallel to the budget row, two equal bounds), the
  multipliers are not unique; one valid set is returned.
- The weights and all multipliers are NaN when the status provides no point: always for `PO_SOLVE_INFEASIBLE` and
  `PO_SOLVE_UNBOUNDED`, and for `PO_SOLVE_MAX_ITER` of every objective (the dual iterate of the QP solver is not primal
  feasible, and the iterate of `RISK_PARITY` is not a point the solver can vouch for). `po_result_residual` is then NaN
  too; otherwise it is 0 for every objective but `PO_OBJECTIVE_RISK_PARITY`, where it is `sum_i (c_i - b_i)^2`.
- For `RISK_PARITY` with `s = sigma_min` the volatility multiplier is `-inf` and the others NaN (see its section).

For MIN_VARIANCE without a box or cash, `lambda_budget = 2 / A`: it is the derivative of the optimal variance with respect to the
budget `b`: the optimal variance is `b^2 / A`, whose derivative at `b = 1` is 2 / A.

## Sources

Verified means the text was read for the claim it supports; the others are marked.

- Tibshirani (1996), "Regression shrinkage and selection via the lasso", JRSS-B 58(1), section 6: the family of
  rows of the leverage. Verified.
- Niedermayer and Niedermayer (2007), "Applying Markowitz's critical line algorithm", Univ. Bern DP 07-01, proposition 1
  and corollary 1: the piecewise affine path of `MAX_RETURN`. Verified.
- Cornuejols and T�t�nc� (2007, section 8.2) and Choueifaty and Coignard (2008): the homogenization of `MAX_SHARPE` and
  `MAX_DIVERSIFICATION`. Not read; cited at second hand. The bound `kappa >= 1 / (L ||v||_inf)` is derived here.
- Bai, Scheinberg and T�t�nc� (2016), Quantitative Finance 16(3): lemmas 2.1 to 2.3, examples 3.1 and 3.2. Verified on
  the preprint (2013); the published version was not read. Lemma 2.1 says positive semidefinite, but its proof needs
  positive definite `Sigma` (a derivation of this work).
- Bruder and Roncalli (2012), "Managing risk exposures using the risk budgeting approach", section 2.3.3: existence and
  uniqueness of the positive-orthant solution. Verified.
- Andreani, Birgin, Martinez and Schuverdt (2007), SIAM J. Optim. 18(4), algorithm 3.1: the augmented Lagrangian.
  Verified on the authors' preprint.
- Feng and Palomar (2015), SCRIP: not read (the certificate of the site was refused); known through the package
  `riskParityPortfolio` and Palomar (2025), chapter 11.
- Spinu (2013): not read (SSRN refused the request); the algorithm is taken from Griveau-Billion, Richard and Roncalli
  (2013), appendix A.3, which was read.
- Ogita, Rump and Oishi (2005), "Accurate sum and dot product", SIAM J. Sci. Comput. 26(6): the compensated dot product
  (Dot2). Not in the research note and not verified against the paper: the implementation follows the algorithm as the
  authors of this work know it, and its accuracy is measured (`validation.md`), not taken from the paper's bound.
- Higham, *Accuracy and Stability of Numerical Algorithms*: the rigorous backward error of Cholesky is stated with
  `|R'| |R|` and not `|Sigma|`, which is why the noise model of the `MAX_SHARPE` guard is a first-order bound and not a
  proof. Not verified in the book.
- Goldfarb and Idnani (1983): the termination argument with removals was not verified in the paper.
