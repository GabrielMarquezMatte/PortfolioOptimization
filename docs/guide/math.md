# Mathematics

The equations PortOpt solves, as the code writes them. This version covers `MIN_VARIANCE`, `UTILITY`,
`IMPORTANCE`, `KELLY`, `TRACKING_ERROR`, `MAX_RETURN`, `MAX_SHARPE` and `MAX_DIVERSIFICATION`. Later objectives add
their own sections.

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

## MAX_RETURN

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
3. It bisects until both ends of the bracket have the same active set, takes the exact root `t*` of
   `|L'w|^2 = s^2` on that affine segment and solves `Q(t*)`. When `Q(t*)` has the same active set (it lies on the
   segment) and `|vol(w(t*)) - s|` is at most `1e-6` of the bracket's volatility range (`kCrossingSanity`), the
   crossing is accepted; otherwise `Q(t*)` narrows the bracket like any other point. A bracket whose ends are
   adjacent doubles is accepted as it is: on one segment by its root, across a breakpoint by its lower end when that
   end's volatility is within `1e-5 s` of `s` (`kAdjacentEnd`); `vol(t)` is continuous, so an end farther away means
   that `Q(t)` is lost in its rounding there, and the path has broken down (`PO_SOLVE_MAX_ITER`).

The reported point is the point of the segment at the root: the weights and the multipliers are the interpolation
of the two ends of the bracket at the fraction of the root. The KKT conditions of `Q(t)` on one active set are linear
in `(w, u, t)`, and the interpolation keeps the nonnegative multipliers nonnegative and the inactive rows feasible,
so that point is an exact solution of `Q(t*)` up to the rounding of the ends, and `vol(w) = s` to rounding. The
rounding of `vol(w(t*))` of a single QP solution, about `1e-15 cond(Sigma) s`, does not enter. (At adjacent ends
across a breakpoint, the lower end is reported, with `lambda_vol = -vol / t`.)

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
  `lambda_vol = -inf` and the budget multiplier (without cash) and the bound multipliers are NaN, although those
  constraints exist. This is the one case where a multiplier is NaN for an existing constraint. When the weights and
  the multipliers disagree, the multipliers decide (`0` or `-s / t_b`), never `-inf`.

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
  NaN weights and multipliers. It needs an unbounded feasible set: with every weight bounded on both sides, or with a
  leverage, `kappa = 0` forces `y = 0`, which contradicts `v'y = 1`. Without bounds it happens when
  `1' Sigma^-1 v <= 0`.
- Otherwise `PO_SOLVE_OPTIMAL` with `w = y / kappa` (and cash exactly 0).

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

## Which setters and data each objective takes

Any other combination is `PO_INVALID_ARGUMENT`:

| Objective | requires | also accepts, besides the bounds |
|---|---|---|
| `MIN_VARIANCE` | covariance (mu with a target) | `target_return` |
| `MAX_RETURN` | mu, covariance, `max_volatility` | none |
| `UTILITY` | mu, covariance | `risk_aversion` |
| `IMPORTANCE` | mu, covariance | `importance` |
| `KELLY` | mu, covariance, `1 + rf > 0` | none |
| `TRACKING_ERROR` | returns (mu with a target) | `target_return` |
| `MAX_SHARPE` | mu, covariance | `leverage` |
| `MAX_DIVERSIFICATION` | covariance | `leverage` |

`po_solve` checks, in this order: the objective is solved by this version (`PO_UNSUPPORTED` otherwise), every set
setter applies to the objective, the data it requires is present, and per-asset bounds match the weight count. The
messages are `"<setter> does not apply to <OBJECTIVE>"`, `"<OBJECTIVE> requires <item>"` (the items are `mu`,
`covariance`, `returns` and `max_volatility`) and `"KELLY requires risk_free > -1"`. The setters are the trust
boundary of the values: `risk_aversion` is positive and finite, `importance` is in [0, 1), `max_volatility` is
positive and finite and `target_return` is finite.

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
Bounds of the risky assets come first, one lower then one upper per asset, then those of cash.

## Multiplier convention

The solver uses `G x + g = E'y + C'z` with `z >= 0` (`src/qp/problem.hpp`). Writing `f(w)` for the objective in the minimization forms above (`w' Sigma w` for MIN_VARIANCE, without the
factor 1/2; a maximized objective enters negated), the stationarity condition at the optimum, over the risky weights,
is

    grad f(w) = lambda_budget 1 + lambda_return grad(return)(w) + lambda_vol Sigma w / vol(w)
                + sum_i beta_i e_i - beta_0 1

- `lambda_budget`, `lambda_return` and `lambda_vol` are `po_result_multiplier` for `PO_CONSTRAINT_BUDGET`,
  `PO_CONSTRAINT_TARGET_RETURN` and `PO_CONSTRAINT_VOLATILITY`. `grad(return)` is `mu`, or `mu - rf 1` with cash,
  and `Sigma w / vol(w)` is the gradient of `vol(w) = sqrt(w' Sigma w)`.
- The term of an absent constraint is omitted: there is no budget term with cash, no return term without a target
  return, no volatility term unless the objective is MAX_RETURN and no `beta_0` without cash.
- `beta_k` is the lower-bound multiplier minus the upper-bound multiplier of weight k (`po_result_bound_multipliers`):
  `>= 0` at an active lower bound and `<= 0` at an active upper bound. It is 0, not NaN, for a weight without a
  finite bound.
- `beta_0` is the entry of cash, the last one. It enters with a minus sign because the cash weight is `1 - 1'w`: the
  cash rows have the opposite signs of those of an asset (see the table).
- A multiplier of an absent constraint is NaN: the budget when there is cash, the target return when none was set
  (and for the objectives that take none), the volatility constraint unless the objective is MAX_RETURN, and the
  leverage constraint, which does not exist in this version. The one other NaN is the vertical `s = sigma_min` case
  of MAX_RETURN above: the budget multiplier (without cash) and the bound multipliers.
- Each multiplier is the derivative of the optimal `f` with respect to the right-hand side of its constraint. A `<=`
  constraint therefore has a multiplier `<= 0`, like `beta` at an upper bound, and a `>=` constraint one `>= 0`. For
  MAX_RETURN, `lambda_vol <= 0` is the derivative with respect to the cap `s`, and 0 when the cap is inactive.
- When the constraints are dependent (a target-return row parallel to the budget row, two equal bounds), the
  multipliers are not unique; one valid set is returned.
- The weights and all multipliers are NaN when the status provides no point: always for `PO_SOLVE_INFEASIBLE` and
  `PO_SOLVE_UNBOUNDED`, and for `PO_SOLVE_MAX_ITER` of every objective solved as a QP (every objective except
  `PO_OBJECTIVE_RISK_PARITY`), because the dual iterate of the QP solver is not primal feasible. `po_result_residual`
  is then NaN too; otherwise it is 0 for every objective but `PO_OBJECTIVE_RISK_PARITY`.

For MIN_VARIANCE without a box or cash, `lambda_budget = 2 / A`: it is the derivative of the optimal variance with respect to the
budget `b`: the optimal variance is `b^2 / A`, whose derivative at `b = 1` is 2 / A.
