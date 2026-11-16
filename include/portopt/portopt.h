#ifndef PORTOPT_PORTOPT_H
#define PORTOPT_PORTOPT_H

#include <stddef.h>
#include <stdint.h>

#if defined(_WIN32)
#if defined(PO_BUILDING)
#define PO_API __declspec(dllexport)
#else
#define PO_API __declspec(dllimport)
#endif
#else
#define PO_API __attribute__((visibility("default")))
#endif

#ifdef __cplusplus
extern "C"
{
#endif

#define PO_ABI_VERSION 0

PO_API int32_t po_abi_version(void);

/* Library version, "major.minor.patch". The pointer is valid for the lifetime of the process. */
PO_API const char* po_version(void);

typedef int32_t po_status;
#define PO_OK 0
#define PO_INVALID_ARGUMENT 1
#define PO_INVALID_HANDLE 2
#define PO_UNSUPPORTED 3
#define PO_OUT_OF_MEMORY 4
#define PO_INTERNAL_ERROR 5

/* Message of the last failed call on the calling thread ("" after a success). The pointer stays valid
   until the next PortOpt call on that thread. */
PO_API const char* po_last_error(void);

/* Handles are opaque. A null handle is PO_INVALID_HANDLE; a null output pointer or a null required input
   array is PO_INVALID_ARGUMENT (the one exception: mu and cov of po_data_create may be NULL, see there).
   After a failed po_*_create or po_solve, *out is NULL (when out itself is not null).
   Every *_free accepts NULL and does nothing. Handles are independent: freeing one never invalidates another.
   Functions that take const handles only read them, so they are safe to call in parallel on the same handles;
   a handle that is being modified (a setter) must not be used by another thread at the same time. */

typedef struct po_spec po_spec;
typedef struct po_data po_data;
typedef struct po_result po_result;

/* Objectives. This version solves every objective: PO_OBJECTIVE_MIN_VARIANCE, PO_OBJECTIVE_MAX_RETURN,
   PO_OBJECTIVE_MAX_SHARPE, PO_OBJECTIVE_UTILITY, PO_OBJECTIVE_IMPORTANCE, PO_OBJECTIVE_KELLY,
   PO_OBJECTIVE_TRACKING_ERROR, PO_OBJECTIVE_MAX_DIVERSIFICATION and PO_OBJECTIVE_RISK_PARITY. Values only
   grow, by addition.

   Each solved objective is the minimization of f over the weights, subject to the budget (the weights sum to
   one, cash included), the weight bounds, the gross exposure limit when po_spec_set_leverage was called
   (sum_i |w_i| <= L over the risky weights, cash not counted) and, for MIN_VARIANCE and TRACKING_ERROR when
   po_spec_set_target_return was called, the target return (mu' w plus rate times the cash weight equals the
   target), with mu~ = mu - rate * 1 (rate = 0 without cash), Sigma the covariance and
   sigma = sqrt(diag Sigma):
     MIN_VARIANCE    f = w' Sigma w
     MAX_RETURN      f = -mu~' w, subject to sqrt(w' Sigma w) <= sigma (po_spec_set_max_volatility)
     MAX_SHARPE      f = -mu~' w / sqrt(w' Sigma w), the Sharpe ratio negated
     UTILITY         f = (delta / 2) w' Sigma w - mu~' w, delta from po_spec_set_risk_aversion (default 5)
     IMPORTANCE      f = (1 - lambda) w' Sigma w - lambda mu~' w, lambda from po_spec_set_importance,
                     0 <= lambda < 1 (default 0.5)
     KELLY           f = w' Sigma w / (2 (1 + rate)) - mu~' w, UTILITY with delta = 1 / (1 + rate); needs
                     1 + rate > 0
     TRACKING_ERROR  f = |R~ w + rate * 1 - b|^2 with R~ = R - rate * 1 1', R the returns (T x n) and b the
                     benchmark; without cash R~ = R and the rate term vanishes
     MAX_DIVERSIFICATION
                     f = -sigma' w / sqrt(w' Sigma w), the diversification ratio negated
     RISK_PARITY     f = sum_i (c_i - b_i)^2 over the risky weights, with c_i = w_i (Sigma w)_i / w' Sigma w
                     the risk contributions and b the risk budget (po_spec_set_risk_budget, default 1 / n); it
                     takes the budget, the bounds, the target return when po_spec_set_target_return was called
                     (it needs mu) and the volatility sqrt(w' Sigma w) = s when po_spec_set_target_volatility
                     was called, but no cash and no leverage (below)
   MAX_SHARPE and MAX_DIVERSIFICATION do not change when the weights are scaled, so they are solved on the
   cone of the feasible set and the result is rescaled to satisfy the budget. With cash, both fix the cash
   weight at 0 (the risky weights then sum to one, and bounds on the cash weight that exclude 0 make the
   problem infeasible); the budget term of the stationarity condition below is then carried by beta_0, and the
   budget multiplier is NaN. A MAX_SHARPE or MAX_DIVERSIFICATION problem with no feasible portfolio of
   positive numerator (for example MAX_SHARPE long only with the rate above every mean) is
   PO_SOLVE_INFEASIBLE; one whose supremum exists but is not attained (for example MAX_SHARPE with unbounded
   short selling and 1' Sigma^-1 mu~ <= 0) is PO_SOLVE_UNBOUNDED. Near that boundary, a problem whose
   1' Sigma^-1 mu~ is positive but within the rounding of the data, so that its exact weights are of the order
   of the inverse of that rounding, is PO_SOLVE_UNBOUNDED too: its sign is not decided by the data. This
   applies only when the feasible set is unbounded; with a leverage, or when no two distinct weights have one
   an infinite upper and the other an infinite lower bound, the result is never PO_SOLVE_UNBOUNDED.
   Which setters and data each objective takes (a setter outside its row makes po_solve return
   PO_INVALID_ARGUMENT, "<setter> does not apply to <OBJECTIVE>"; a missing input, "<OBJECTIVE> requires
   <item>"):

     objective            requires                       also accepts, besides the bounds   leverage
     MIN_VARIANCE         covariance (mu with a target)  target_return                      yes
     MAX_RETURN           mu, covariance, max_volatility -                                  yes
     MAX_SHARPE           mu, covariance                 -                                  yes
     UTILITY              mu, covariance                 risk_aversion                      yes
     IMPORTANCE           mu, covariance                 importance                         yes
     KELLY                mu, covariance, rate > -1      -                                  yes
     TRACKING_ERROR       returns (mu with a target)     target_return                      yes
     MAX_DIVERSIFICATION  covariance                     -                                  yes
     RISK_PARITY          covariance (mu with a target)  target_return, target_volatility,  no
                                                         risk_budget

   po_spec_set_leverage is accepted by every objective except PO_OBJECTIVE_RISK_PARITY ("leverage does not
   apply to RISK_PARITY"). PO_OBJECTIVE_RISK_PARITY has no cash either: po_data_set_risk_free makes po_solve
   return PO_INVALID_ARGUMENT, "risk_free does not apply to RISK_PARITY". po_spec_set_target_volatility
   applies only to PO_OBJECTIVE_RISK_PARITY ("target_volatility does not apply to <OBJECTIVE>" elsewhere), and
   a risk budget with k entries on data with n assets is "risk_budget has k entries; expected n".

   PO_OBJECTIVE_RISK_PARITY is solved by a local method. PO_SOLVE_OPTIMAL means a KKT point of f over the
   constraints, reached from the positive-orthant Spinu solution (f = 0), not the global minimum: f is not
   convex, and in general there are at most 2^(n-1) normalized solutions of c = b, one per pair of opposite
   orthants, none when the weight sum of the ray is 0. Without a target return or a target volatility, with
   the Spinu point inside the bounds, the result is that point, every multiplier 0 and the residual at the
   rounding floor. A KKT point is accepted when its KKT residual and its feasibility are within the rounding
   noise of the problem. With a target volatility and few or no active bounds at cond(Sigma) around 1e8 the
   weights are determined only to the rounding of the KKT residual amplified by the conditioning of the
   problem: the weights of an OPTIMAL result can be up to about 1e-4 relative from the exact KKT point
   (measured up to 3.7e-5), while its objective value is right to 1e-7 (1 + f). That covers the weights and
   the value, not the multipliers: with a target volatility the multipliers are measured right only to about
   5e-10 relative at cond(Sigma) = 1e4, 3e-7 at 1e6 and 4e-2 at 1e8, and can be wrong by O(1) above that
   (an OPTIMAL result at cond(Sigma) >= 1e10 is not validated). Such a solve is also slow: about 10 s at
   n = 50, cond(Sigma) = 1e8, no box and a target volatility, and at most about 0.15 s for n <= 20.
   The Spinu point is computed before the feasibility check, so a Spinu iteration that does not converge
   reports PO_SOLVE_MAX_ITER even for an infeasible problem. */
#define PO_OBJECTIVE_MIN_VARIANCE 0
#define PO_OBJECTIVE_MAX_RETURN 1
#define PO_OBJECTIVE_MAX_SHARPE 2
#define PO_OBJECTIVE_UTILITY 3
#define PO_OBJECTIVE_IMPORTANCE 4
#define PO_OBJECTIVE_KELLY 5
#define PO_OBJECTIVE_TRACKING_ERROR 6
#define PO_OBJECTIVE_MAX_DIVERSIFICATION 7
#define PO_OBJECTIVE_RISK_PARITY 8

/* Solve status of a result. An infeasible problem is a result with this status, not an error of po_solve.
   PO_SOLVE_INFEASIBLE is certified: incompatible linear constraints, a target return outside the achievable
   range or, for PO_OBJECTIVE_RISK_PARITY, a target volatility below the minimum volatility of the feasible
   set. PO_SOLVE_MAX_ITER means the solver did not end at a point it can vouch for, and it carries no point
   (see po_result below). For PO_OBJECTIVE_RISK_PARITY that includes a target volatility the iteration cannot
   reach (for example above every volatility of a bounded set), a target return without a minimizer (no
   bounds), and a solution whose gross exposure sum_i |w_i| is above about 4.5e3, where the budget row 1'w = 1
   is no longer resolved in double (except at s = sigma_min, where the minimum-variance portfolio is returned
   without this check). A bounded problem with the budget row only (n = 20, cond(Sigma) = 1e8) can also end
   here: a known limitation. */
#define PO_SOLVE_OPTIMAL 0
#define PO_SOLVE_INFEASIBLE 1
#define PO_SOLVE_UNBOUNDED 2 /* the supremum exists but is not attained; weights and multipliers are NaN */
#define PO_SOLVE_MAX_ITER 3

/* Named constraints, for po_result_multiplier. */
#define PO_CONSTRAINT_BUDGET 0
#define PO_CONSTRAINT_TARGET_RETURN 1
#define PO_CONSTRAINT_VOLATILITY 2
#define PO_CONSTRAINT_LEVERAGE 3

/* Specification. A new spec has no weight bounds (every weight in (-inf, +inf)) and no target return. */

/* PO_INVALID_ARGUMENT if objective is not one of PO_OBJECTIVE_*. */
PO_API po_status po_spec_create(int32_t objective, po_spec** out);

/* Bounds on every weight, cash included. lower may be -inf and upper +inf; NaN, lower > upper, lower = +inf
   and upper = -inf are PO_INVALID_ARGUMENT and leave the spec unchanged. Clears any per-asset bounds: the
   last call of po_spec_set_bounds or po_spec_set_asset_bounds is the one that applies. */
PO_API po_status po_spec_set_bounds(po_spec* s, double lower, double upper);

/* One bound per weight: n must be the weight count of the data given to po_solve (the asset count, plus one
   for cash, which comes last), checked by po_solve. Each pair follows the rules of po_spec_set_bounds, and
   nothing is changed unless every pair is valid. Replaces the scalar bounds; n = 0 is PO_INVALID_ARGUMENT. */
PO_API po_status po_spec_set_asset_bounds(po_spec* s, const double* lower, const double* upper, size_t n);

/* Required portfolio return. r must be finite. */
PO_API po_status po_spec_set_target_return(po_spec* s, double r);

/* Risk aversion delta of PO_OBJECTIVE_UTILITY. delta must be positive and finite. */
PO_API po_status po_spec_set_risk_aversion(po_spec* s, double delta);

/* Weight lambda of the return in PO_OBJECTIVE_IMPORTANCE. 0 <= lambda < 1 (NaN is rejected). */
PO_API po_status po_spec_set_importance(po_spec* s, double lambda);

/* Volatility cap sigma of PO_OBJECTIVE_MAX_RETURN. sigma must be positive and finite.
   An invalid value of any of the three setters above leaves the spec unchanged and is PO_INVALID_ARGUMENT. */
PO_API po_status po_spec_set_max_volatility(po_spec* s, double sigma);

/* Volatility target sigma of PO_OBJECTIVE_RISK_PARITY: sqrt(w' Sigma w) = sigma. sigma must be positive and
   finite; an invalid value leaves the spec unchanged and is PO_INVALID_ARGUMENT. */
PO_API po_status po_spec_set_target_volatility(po_spec* s, double sigma);

/* Risk budget b of PO_OBJECTIVE_RISK_PARITY: n entries, each positive and finite, with |sum b - 1| <= n *
   DBL_EPSILON (the rounding of a sum of n terms). n must be the asset count of the data given to po_solve
   (checked there). A null b, n = 0 or a value outside the domain is PO_INVALID_ARGUMENT and leaves the spec
   unchanged. Without a budget it is 1 / n. */
PO_API po_status po_spec_set_risk_budget(po_spec* s, const double* b, size_t n);

/* Gross exposure limit: sum_i |w_i| <= max_gross over the risky weights (the cash weight is not counted).
   max_gross must be positive and finite; zero or negative, NaN and +-inf are PO_INVALID_ARGUMENT and leave
   the spec unchanged. Without cash a limit below 1 is infeasible, because sum_i |w_i| >= |1' w| = 1; po_solve
   reports it as PO_SOLVE_INFEASIBLE. The same holds with cash for PO_OBJECTIVE_MAX_SHARPE and
   PO_OBJECTIVE_MAX_DIVERSIFICATION, whose cash weight is fixed at 0, so that the risky weights sum to one. */
PO_API po_status po_spec_set_leverage(po_spec* s, double max_gross);

PO_API void po_spec_free(po_spec* s);

/* Data. The inputs are copied. mu (n entries) and cov (n x n, row-major) may each be NULL when the objective
   does not need them; po_solve and the portfolio functions report a missing one as PO_INVALID_ARGUMENT.
   n = 0, a non-finite entry, an asymmetric covariance or one that is not positive definite is
   PO_INVALID_ARGUMENT, with the reason in po_last_error(). */
PO_API po_status po_data_create(size_t n, const double* mu, const double* cov, po_data** out);

/* Adds a risk-free asset (cash) as the last weight: the weight count becomes the asset count plus one.
   Weights always sum to one, cash included. rate must be finite. */
PO_API po_status po_data_set_risk_free(po_data* d, double rate);

/* Return series for PO_OBJECTIVE_TRACKING_ERROR. returns is t x n (row-major, n = the asset count, cash
   excluded) and benchmark has t entries; both are copied and replace any earlier series. NULL returns or
   benchmark, t = 0, a size that overflows or a non-finite value is PO_INVALID_ARGUMENT, with the reason in
   po_last_error(), and the earlier series is kept. */
PO_API po_status po_data_set_returns(po_data* d, const double* returns, const double* benchmark, size_t t);

PO_API po_status po_data_asset_count(const po_data* d, size_t* out);

/* The asset count, plus one when a risk-free rate was set. */
PO_API po_status po_data_weight_count(const po_data* d, size_t* out);

PO_API void po_data_free(po_data* d);

/* Solves the spec on the data; both are copied, so the result does not depend on them afterwards. Every
   objective is solved, so PO_UNSUPPORTED does not occur for a valid objective. It checks, in this order, and
   returns PO_INVALID_ARGUMENT at the first failure: a setter outside the objective's row, data the objective
   requires, per-asset bounds that do not match the weight count. The core
   does not range-check the values: the setters are the trust boundary. Infeasibility is reported in the
   result. Safe in parallel with the same read-only spec and data. */
PO_API po_status po_solve(const po_spec* s, const po_data* d, po_result** out);

/* Result, immutable. The weights, the bound multipliers, the multipliers of the named constraints and the
   residual are NaN when the status provides no point: always for PO_SOLVE_INFEASIBLE and PO_SOLVE_UNBOUNDED,
   and for PO_SOLVE_MAX_ITER of every objective (the dual iterate of the QP solver is not primal feasible, and
   the iterate of RISK_PARITY is not a point the solver can vouch for).

   Multiplier convention. Let f be the objective in minimization form, as written for each objective above:
   for PO_OBJECTIVE_MIN_VARIANCE f(w) = w' Sigma w, without the factor 1/2, and a maximized objective enters
   negated. At the optimum, over the risky weights w (the first po_data_asset_count entries),
       grad f(w) = lambda_budget * 1 + lambda_return * grad(return) + lambda_vol * Sigma w / vol(w)
                   + lambda_lev * g + sum_i beta_i e_i - beta_0 * 1
   where lambda_budget, lambda_return, lambda_vol and lambda_lev are po_result_multiplier for
   PO_CONSTRAINT_BUDGET, PO_CONSTRAINT_TARGET_RETURN, PO_CONSTRAINT_VOLATILITY and PO_CONSTRAINT_LEVERAGE,
   vol(w) = sqrt(w' Sigma w), so Sigma w / vol(w) is the gradient of the volatility, grad(return) is mu
   (mu - rate * 1 with cash), g is in the subdifferential of the 1-norm of w: g_i = sign(w_i) where w_i is not
   0 and any value in [-1, 1] where w_i = 0, beta is po_result_bound_multipliers
   (lower-bound multiplier minus upper-bound multiplier: >= 0 at an active lower bound, <= 0 at an active
   upper bound) and beta_0 is its cash entry, which enters with a minus sign because the cash weight is
   1 - 1'w. The term of an absent constraint is omitted: no budget term with cash, no return term without a
   target return, no volatility term unless the objective is PO_OBJECTIVE_MAX_RETURN or
   PO_OBJECTIVE_RISK_PARITY with a target volatility, no leverage term without a leverage, no
   beta_0 without cash. For PO_OBJECTIVE_MAX_SHARPE and PO_OBJECTIVE_MAX_DIVERSIFICATION with cash the cash
   weight is 0 and not a variable, so the budget term appears as beta_0 and lambda_budget is NaN.

   Each multiplier is the derivative of the optimal f with respect to the right-hand side of its constraint,
   so a '<=' constraint has a multiplier <= 0, like beta at an upper bound. When the constraints are dependent
   the multipliers are not unique; one valid set is returned.

   For PO_OBJECTIVE_MAX_RETURN, PO_CONSTRAINT_VOLATILITY is lambda_vol <= 0, the derivative of the optimal f
   with respect to the cap sigma; it is 0 when the cap is inactive.

   For PO_OBJECTIVE_RISK_PARITY, PO_CONSTRAINT_VOLATILITY is the derivative of the optimal f with respect to
   the target volatility s, of either sign, and NaN when po_spec_set_target_volatility was not called (its
   accuracy is in the PO_OBJECTIVE_RISK_PARITY entry above: below that of the weights). The budget
   multiplier is finite (no cash), the target-return multiplier is NaN without a target return, the leverage
   multiplier is always NaN, and a bound multiplier is 0 for a weight with no finite bound. The risk
   contributions are homogeneous of degree 0, so w' grad f(w) = 0 and
       lambda_budget = -lambda_return * r - lambda_vol * s - sum_i beta_i w_i
   (a term of an absent constraint is 0). When s equals the minimum volatility of the feasible set to 1e-12
   relative, the minimum-variance portfolio is the only feasible point: the result is PO_SOLVE_OPTIMAL there,
   with lambda_vol = -inf and every other multiplier NaN, because the constraint gradients are dependent and
   no KKT multipliers exist (the same convention as PO_OBJECTIVE_MAX_RETURN below).

   PO_CONSTRAINT_LEVERAGE is lambda_lev <= 0, the derivative of the optimal f with respect to the limit; it is
   0 when the leverage is set but inactive.

   A bound multiplier is 0 (not NaN) for a weight with no finite bound. Apart from the two exceptions below,
   NaN is only for an absent named constraint: the budget constraint does not exist when there is cash (nor
   for PO_OBJECTIVE_MAX_SHARPE and PO_OBJECTIVE_MAX_DIVERSIFICATION with cash, see above), the target-return
   constraint when no target was set (and for the objectives that take none), the volatility constraint unless
   the objective is PO_OBJECTIVE_MAX_RETURN or PO_OBJECTIVE_RISK_PARITY with a target volatility, and the
   leverage constraint when po_spec_set_leverage was not called (always, for PO_OBJECTIVE_RISK_PARITY).

   Two documented exceptions, where existing constraints have NaN multipliers: for PO_OBJECTIVE_MAX_RETURN,
   when the cap equals the minimum volatility of the feasible set and the efficient frontier is vertical
   there, lambda_vol = -inf, and the budget multiplier (when present), the bound multipliers and the leverage
   multiplier (when a leverage was set, active or not) are NaN; for PO_OBJECTIVE_RISK_PARITY at s equal to the
   minimum volatility, lambda_vol = -inf and the budget and bound multipliers are NaN (above). */
PO_API po_status po_result_status(const po_result* r, int32_t* out);

/* The weight count of the data the result was solved on. */
PO_API po_status po_result_weight_count(const po_result* r, size_t* out);

/* The vector accessors require n to be exactly the weight count, otherwise PO_INVALID_ARGUMENT. */
PO_API po_status po_result_weights(const po_result* r, double* out, size_t n);

/* Total iterations of the solver, counted in a method-specific way; for diagnostics only. */
PO_API po_status po_result_iterations(const po_result* r, int32_t* out);

/* constraint is one of PO_CONSTRAINT_*, otherwise PO_INVALID_ARGUMENT. */
PO_API po_status po_result_multiplier(const po_result* r, int32_t constraint, double* out);
PO_API po_status po_result_bound_multipliers(const po_result* r, double* out, size_t n);

/* For PO_OBJECTIVE_RISK_PARITY, sum_i (c_i - b_i)^2 at the reported point, with
   c_i = w_i (Sigma w)_i / w' Sigma w the risk contributions and b the risk budget; 0 for every other
   objective; NaN when the weights are NaN. */
PO_API po_status po_result_residual(const po_result* r, double* out);
PO_API void po_result_free(po_result* r);

/* Portfolio metrics for weights w (n must be the weight count of d; cash last). The return is
   sum_i w_i mu_i + w_cash * rate and needs mu; the volatility is sqrt(w' Sigma w) over the risky weights only
   and needs the covariance. A missing input is PO_INVALID_ARGUMENT. */
PO_API po_status po_portfolio_return(const po_data* d, const double* w, size_t n, double* out);
PO_API po_status po_portfolio_volatility(const po_data* d, const double* w, size_t n, double* out);

/* return - delta / 2 w' Sigma w, with the return as above. Needs mu and the covariance, and delta must be
   positive and finite; otherwise PO_INVALID_ARGUMENT. */
PO_API po_status po_utility(const po_data* d, const double* w, size_t n, double delta, double* out);

/* sqrt(sum_t (r_t'w + w_0 rf - b_t)^2 / T): the root mean square the TRACKING_ERROR objective minimizes, with
   r_t the returns of period t over the risky weights, w_0 the cash weight and rf the rate (the cash term
   vanishes without cash), b the benchmark. Needs the returns (po_data_set_returns); otherwise
   PO_INVALID_ARGUMENT. */
PO_API po_status po_tracking_error(const po_data* d, const double* w, size_t n, double* out);

/* (return - rf) / volatility over the risky weights, with the return and the volatility as above and rf the
   rate (rf = 0 without cash). Needs mu and the covariance, otherwise PO_INVALID_ARGUMENT; NaN when the
   volatility is 0. For weights that sum to one it is mu~' w / vol(w), the quantity PO_OBJECTIVE_MAX_SHARPE
   maximizes. */
PO_API po_status po_sharpe_ratio(const po_data* d, const double* w, size_t n, double* out);

/* sigma' w / volatility over the risky weights, with sigma = sqrt(diag Sigma): the quantity
   PO_OBJECTIVE_MAX_DIVERSIFICATION maximizes. Needs the covariance (not mu), otherwise PO_INVALID_ARGUMENT;
   NaN when the volatility is 0. */
PO_API po_status po_diversification_ratio(const po_data* d, const double* w, size_t n, double* out);

/* Risk contributions c_i = w_i (Sigma w)_i / w' Sigma w over the risky weights: out has n entries (n must be
   the weight count of d, otherwise PO_INVALID_ARGUMENT), the cash entry is 0 and the entries sum to 1. All
   NaN when the volatility is 0. Needs the covariance (not mu), otherwise PO_INVALID_ARGUMENT. */
PO_API po_status po_risk_contributions(const po_data* d, const double* w, size_t n, double* out);

#ifdef __cplusplus
}
#endif

#endif
