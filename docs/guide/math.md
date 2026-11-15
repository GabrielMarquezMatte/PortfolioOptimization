# Mathematics

The equations PortOpt solves, as the code writes them. This version covers `MIN_VARIANCE`, `UTILITY`,
`IMPORTANCE`, `KELLY`, `TRACKING_ERROR` and `MAX_RETURN`. Later objectives add their own sections.

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

    Q(t):  minimize  0.5 w' Sigma w - t mu~'w    over the budget, box and cash rows

(`G = Sigma` with the factor `L`, `g = -t mu~`; only `g` depends on `t`). Its optimum `w(t)` is piecewise affine in
`t` and its volatility is nondecreasing, from the minimum volatility `sigma_min` at `t = 0`. The solver finds the
`t*` where the volatility reaches `s`:

1. At `t = 0` it solves the minimum-variance point. `s` below `sigma_min` is `PO_SOLVE_INFEASIBLE`. Feasibility is
   decided there, once: a `Q(t)` that fails at `t > 0` has broken down along the path and is `PO_SOLVE_MAX_ITER`,
   never `PO_SOLVE_INFEASIBLE`.
2. It doubles `t`, from `t_1 = sigma_min^2 / |mu~|max` (with `sigma_min = 0`, all cash, from
   `min(s^2, max_i Sigma_ii) / |mu~|max`, so that a large cap does not start far past the last breakpoint), until
   the volatility reaches `s`, or until `w(t)` stops changing (the inactive cap below).
3. It bisects until both ends of the bracket have the same active set, takes the exact root `t*` of
   `|L'w|^2 = s^2` on that affine segment and solves `Q(t*)`. When `Q(t*)` has the same active set (it lies on the
   segment) and `|vol(w(t*)) - s|` is at most `1e-6` of the bracket's volatility range (`kCrossingSanity`), the
   crossing is accepted; otherwise `Q(t*)` narrows the bracket like any other point. A bracket whose ends are
   adjacent doubles is accepted as it is.

The reported point is the point of the segment at the root: the weights and the multipliers are the interpolation
of the two ends of the bracket at the fraction of the root. The KKT conditions of `Q(t)` on one active set are linear
in `(w, u, t)`, and the interpolation keeps the nonnegative multipliers nonnegative and the inactive rows feasible,
so that point is an exact solution of `Q(t*)` up to the rounding of the ends, and `vol(w) = s` to rounding. The
rounding of `vol(w(t*))` of a single QP solution, about `1e-15 cond(Sigma) s`, does not enter. (At adjacent ends
across a breakpoint, the lower end is reported, with `lambda_vol = -vol / t`.)

At `t*`, `Sigma w - t* mu~ = N u`, where `N u` collects the QP multipliers `u` of the rows. Dividing by `t*`:

    -mu~ = N (u / t*) - (1 / t*) Sigma w = N (u / t*) + lambda_vol Sigma w / vol(w),    lambda_vol = -s / t*

since `vol(w) = s`. The multipliers of the budget and bound rows are the QP's divided by `t*`, and the volatility
multiplier is `lambda_vol = -s / t*`, the derivative of the optimal `f = -mu~'w` with respect to `s` (<= 0).

**When `w(t)` is constant.** Two points of the path on the same active set have the same weights when the equality
rows and the active rows fix `w`: a vertex, decided combinatorially (the solver keeps the active normals independent,
so their rank is their count, and a count of `n` is a vertex). Otherwise the weights are compared with the rounding
the solver carries: it starts from the unconstrained minimizer `t Sigma^-1 mu~`, so a weight is known to

    16 n eps (max(1, |w|max) + t |Sigma^-1 mu~|max)

(the weight scale floored at 1, plus the size of that minimizer, which grows with `t |mu~|`; `16` is
`kNoiseFactor` in `src/problems/parametric_path.cpp`). Where `w` is constant, `Sigma w - t mu~ = N u(t)` for every `t`, so
`u(t) = u(t_b) + (t - t_b) b` with `N b = -mu~`. The slope `b` is computed from that equation by Householder least
squares over the active normals, not by differencing two QP solutions (whose rounding grows with `t`); a slope below
`-16 n eps |mu~|max` is a decrease, and an active multiplier that decreases reaches zero at the vertex exit
`t_b + u_k(t_b) / |b_k|`. A combination that comes out negative by rounding (within those bounds) is set to 0.

Special cases:

- **Inactive cap.** When `w(t)` is constant between two consecutive points of the climb and no active multiplier
  decreases, `w` stays optimal for every larger `t` and the volatility never reaches `s`. The weights are the
  maximum-return point of the box and `lambda_vol = 0`; the other multipliers are the slopes `b`.
- **`s = sigma_min`** (to `1e-12` relative). Let `t_b = sup{t : w(t) = w(0)}`; `t` is halved from `t_1` until `w(t)` is
  constant with `w(0)`. If `t_b` is infinite (no decrease), `lambda_vol = 0`. If `t_b` is finite (a vertex exit),
  `lambda_vol = -s / t_b`. The efficient frontier is vertical at `sigma_min` (`t_b = 0`) only when the weights moved
  beyond their rounding at a larger `t` on the active set of `w(0)` and no exit comes before that `t`: then
  `lambda_vol = -inf` and the budget multiplier (without cash) and the bound multipliers are NaN, although those
  constraints exist. This is the one case where a multiplier is NaN for an existing constraint. When the weights and
  the multipliers disagree, the multipliers decide (`0` or `-s / t_b`), never `-inf`.

Without a box, `MAX_RETURN` lies on the efficient frontier. With cash it is the capital market line,
`w = s Sigma^-1 mu~ / sqrt(mu~' Sigma^-1 mu~)` over the risky assets, and `lambda_vol = -sqrt(mu~' Sigma^-1 mu~)`,
the negative of the slope of the line.

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
