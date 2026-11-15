# Mathematics

The equations PortOpt solves, as the code writes them. This version covers `MIN_VARIANCE`, the only objective
solved so far. Later objectives add their own sections.

## Notation

- N risky assets. `w` in R^N are the weights of the risky assets, `mu` in R^N the expected returns and `Sigma` the
  N x N covariance matrix, symmetric positive definite. `sigma_i = sqrt(Sigma_ii)`.
- **Cash (the risk-free asset).** When the data carries a risk-free rate `rf` (`po_data_set_risk_free`), there is a
  further asset N + 1 with return `rf` and zero variance and covariances. Its weight is `w_0 = 1 - 1'w`, and it is
  always the last entry of every weight vector. `Sigma` stays that of the risky assets only and is never
  regularized. The return of a portfolio is `rf + (mu - rf 1)'w` and its variance is `w' Sigma w`. Without cash the
  budget constraint `1'w = 1` holds.
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

## The QP rows

`qp::Problem` is `min 0.5 x'Gx + g'x` subject to `E x = e` and `C x >= c`. MIN_VARIANCE sets the variables to the
risky weights, `x = w`, with

    G = 2 Sigma,    g = 0

so that `0.5 w'Gw = w' Sigma w` exactly and the multipliers of the QP are the user's (see below). The rows come from
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

The solver uses `G x + g = E'y + C'z` with `z >= 0` (`src/qp/problem.hpp`). Writing `f(w) = w' Sigma w` (the
minimization form, without the factor 1/2; a maximized objective would enter negated), the stationarity condition
at the optimum, over the risky weights, is

    grad f(w) = lambda_budget 1 + lambda_return grad(return)(w) + sum_i beta_i e_i - beta_0 1

- `lambda_budget` and `lambda_return` are `po_result_multiplier` for `PO_CONSTRAINT_BUDGET` and
  `PO_CONSTRAINT_TARGET_RETURN`. `grad(return)` is `mu`, or `mu - rf 1` with cash.
- The term of an absent constraint is omitted: there is no budget term with cash, no return term without a target
  return and no `beta_0` without cash.
- `beta_k` is the lower-bound multiplier minus the upper-bound multiplier of weight k (`po_result_bound_multipliers`):
  `>= 0` at an active lower bound and `<= 0` at an active upper bound. It is 0, not NaN, for a weight without a
  finite bound.
- `beta_0` is the entry of cash, the last one. It enters with a minus sign because the cash weight is `1 - 1'w`: the
  cash rows have the opposite signs of those of an asset (see the table).
- A multiplier of an absent constraint is NaN: the budget when there is cash, the target return when none was set,
  and the volatility and leverage constraints, which do not exist in this version.
- Each multiplier is the derivative of the optimal `f` with respect to the right-hand side of its constraint. A `<=`
  constraint therefore has a multiplier `<= 0`, like `beta` at an upper bound, and a `>=` constraint one `>= 0`.
- When the constraints are dependent (a target-return row parallel to the budget row, two equal bounds), the
  multipliers are not unique; one valid set is returned.
- The weights and all multipliers are NaN when the status provides no point: always for `PO_SOLVE_INFEASIBLE` and
  `PO_SOLVE_UNBOUNDED`, and for `PO_SOLVE_MAX_ITER` of every objective solved as a QP (every objective except
  `PO_OBJECTIVE_RISK_PARITY`), because the dual iterate of the QP solver is not primal feasible. `po_result_residual`
  is then NaN too; otherwise it is 0 for every objective but `PO_OBJECTIVE_RISK_PARITY`.

Without a box or cash, `lambda_budget = 2 / A`: it is the derivative of the optimal variance with respect to the
budget `b`: the optimal variance is `b^2 / A`, whose derivative at `b = 1` is 2 / A.
