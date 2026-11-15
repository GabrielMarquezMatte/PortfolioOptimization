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

/* Handles are opaque. A null handle is PO_INVALID_HANDLE; a null output pointer or a null required input array is
   PO_INVALID_ARGUMENT (the one exception: mu and cov of po_data_create may be NULL, see there). After a failed
   po_*_create or po_solve, *out is NULL (when out itself is not null).
   Every *_free accepts NULL and does nothing. Handles are independent: freeing one never invalidates another.
   Functions that take const handles only read them, so they are safe to call in parallel on the same handles;
   a handle that is being modified (a setter) must not be used by another thread at the same time. */

typedef struct po_spec po_spec;
typedef struct po_data po_data;
typedef struct po_result po_result;

/* Objectives. This version solves PO_OBJECTIVE_MIN_VARIANCE, PO_OBJECTIVE_MAX_RETURN, PO_OBJECTIVE_UTILITY,
   PO_OBJECTIVE_IMPORTANCE, PO_OBJECTIVE_KELLY and PO_OBJECTIVE_TRACKING_ERROR; po_solve returns
   PO_UNSUPPORTED for PO_OBJECTIVE_MAX_SHARPE, PO_OBJECTIVE_MAX_DIVERSIFICATION and PO_OBJECTIVE_RISK_PARITY.
   Values only grow, by addition.

   Each solved objective is the minimization of f over the weights, subject to the budget (the weights sum to
   one, cash included), the weight bounds and, for MIN_VARIANCE and TRACKING_ERROR when
   po_spec_set_target_return was called, the target return (mu' w plus rate times the cash weight equals the
   target), with mu~ = mu - rate * 1 (rate = 0 without cash) and Sigma the covariance:
     MIN_VARIANCE    f = w' Sigma w
     MAX_RETURN      f = -mu~' w, subject to sqrt(w' Sigma w) <= sigma (po_spec_set_max_volatility)
     UTILITY         f = (delta / 2) w' Sigma w - mu~' w, delta from po_spec_set_risk_aversion (default 5)
     IMPORTANCE      f = (1 - lambda) w' Sigma w - lambda mu~' w, lambda from po_spec_set_importance,
                     0 <= lambda < 1 (default 0.5)
     KELLY           f = w' Sigma w / (2 (1 + rate)) - mu~' w, UTILITY with delta = 1 / (1 + rate); needs
                     1 + rate > 0
     TRACKING_ERROR  f = |R~ w + rate * 1 - b|^2 with R~ = R - rate * 1 1', R the returns (T x n) and b the
                     benchmark; without cash R~ = R and the rate term vanishes
   Which setters and data each objective takes (a setter outside its row makes po_solve return
   PO_INVALID_ARGUMENT, "<setter> does not apply to <OBJECTIVE>"; a missing input, "<OBJECTIVE> requires
   <item>"):

     objective       requires                       also accepts, besides the bounds
     MIN_VARIANCE    covariance (mu with a target)  target_return
     MAX_RETURN      mu, covariance, max_volatility -
     UTILITY         mu, covariance                 risk_aversion
     IMPORTANCE      mu, covariance                 importance
     KELLY           mu, covariance, rate > -1      -
     TRACKING_ERROR  returns (mu with a target)     target_return */
#define PO_OBJECTIVE_MIN_VARIANCE 0
#define PO_OBJECTIVE_MAX_RETURN 1
#define PO_OBJECTIVE_MAX_SHARPE 2
#define PO_OBJECTIVE_UTILITY 3
#define PO_OBJECTIVE_IMPORTANCE 4
#define PO_OBJECTIVE_KELLY 5
#define PO_OBJECTIVE_TRACKING_ERROR 6
#define PO_OBJECTIVE_MAX_DIVERSIFICATION 7
#define PO_OBJECTIVE_RISK_PARITY 8

/* Solve status of a result. An infeasible problem is a result with this status, not an error of po_solve. */
#define PO_SOLVE_OPTIMAL 0
#define PO_SOLVE_INFEASIBLE 1
#define PO_SOLVE_UNBOUNDED 2
#define PO_SOLVE_MAX_ITER 3

/* Named constraints, for po_result_multiplier. */
#define PO_CONSTRAINT_BUDGET 0
#define PO_CONSTRAINT_TARGET_RETURN 1
#define PO_CONSTRAINT_VOLATILITY 2
#define PO_CONSTRAINT_LEVERAGE 3

/* Specification. A new spec has no weight bounds (every weight in (-inf, +inf)) and no target return. */

/* PO_INVALID_ARGUMENT if objective is not one of PO_OBJECTIVE_*. */
PO_API po_status po_spec_create(int32_t objective, po_spec** out);

/* Bounds on every weight, cash included. lower may be -inf and upper +inf; NaN, lower > upper, lower = +inf and
   upper = -inf are PO_INVALID_ARGUMENT and leave the spec unchanged. Clears any per-asset bounds: the last call of
   po_spec_set_bounds or po_spec_set_asset_bounds is the one that applies. */
PO_API po_status po_spec_set_bounds(po_spec* s, double lower, double upper);

/* One bound per weight: n must be the weight count of the data given to po_solve (the asset count, plus one for
   cash, which comes last), checked by po_solve. Each pair follows the rules of po_spec_set_bounds, and nothing is
   changed unless every pair is valid. Replaces the scalar bounds; n = 0 is PO_INVALID_ARGUMENT. */
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

PO_API void po_spec_free(po_spec* s);

/* Data. The inputs are copied. mu (n entries) and cov (n x n, row-major) may each be NULL when the objective does
   not need them; po_solve and the portfolio functions report a missing one as PO_INVALID_ARGUMENT. n = 0, a
   non-finite entry, an asymmetric covariance or one that is not positive definite is PO_INVALID_ARGUMENT, with the
   reason in po_last_error(). */
PO_API po_status po_data_create(size_t n, const double* mu, const double* cov, po_data** out);

/* Adds a risk-free asset (cash) as the last weight: the weight count becomes the asset count plus one. Weights
   always sum to one, cash included. rate must be finite. */
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

/* Solves the spec on the data; both are copied, so the result does not depend on them afterwards. Returns
   PO_UNSUPPORTED for an objective this version does not solve (see the objectives above). Otherwise it
   checks, in this order, and returns PO_INVALID_ARGUMENT at the first failure: a setter outside the
   objective's row, data the objective requires, per-asset bounds that do not match the weight count. The core
   does not range-check the values: the setters are the trust boundary. Infeasibility is reported in the
   result. Safe in parallel with the same read-only spec and data. */
PO_API po_status po_solve(const po_spec* s, const po_data* d, po_result** out);

/* Result, immutable. The weights, the bound multipliers and the multipliers of the named constraints are NaN when
   the status provides no point: always for PO_SOLVE_INFEASIBLE and PO_SOLVE_UNBOUNDED, and for PO_SOLVE_MAX_ITER
   of every objective solved as a QP (every objective except PO_OBJECTIVE_RISK_PARITY), because the dual iterate of
   the QP solver is not primal feasible. What PO_SOLVE_MAX_ITER reports for PO_OBJECTIVE_RISK_PARITY is documented
   when that objective is added.

   Multiplier convention. Let f be the objective in minimization form, as written for each objective above:
   for PO_OBJECTIVE_MIN_VARIANCE f(w) = w' Sigma w, without the factor 1/2, and a maximized objective enters
   negated. At the optimum, over the risky weights w (the first po_data_asset_count entries),
       grad f(w) = lambda_budget * 1 + lambda_return * grad(return) + lambda_vol * Sigma w / vol(w)
                   + sum_i beta_i e_i - beta_0 * 1
   where lambda_budget, lambda_return and lambda_vol are po_result_multiplier for PO_CONSTRAINT_BUDGET,
   PO_CONSTRAINT_TARGET_RETURN and PO_CONSTRAINT_VOLATILITY, vol(w) = sqrt(w' Sigma w), so Sigma w / vol(w) is
   the gradient of the volatility, grad(return) is mu (mu - rate * 1 with cash), beta is
   po_result_bound_multipliers
   (lower-bound multiplier minus upper-bound multiplier: >= 0 at an active lower bound, <= 0 at an active upper
   bound) and beta_0 is its cash entry, which enters with a minus sign because the cash weight is 1 - 1'w. The term
   of an absent constraint is omitted: no budget term with cash, no return term without a target return, no
   volatility term unless the objective is PO_OBJECTIVE_MAX_RETURN, no beta_0 without cash.

   Each multiplier is the derivative of the optimal f with respect to the right-hand side of its constraint, so a
   '<=' constraint has a multiplier <= 0, like beta at an upper bound. When the constraints are dependent the
   multipliers are not unique; one valid set is returned.

   For PO_OBJECTIVE_MAX_RETURN, PO_CONSTRAINT_VOLATILITY is lambda_vol <= 0, the derivative of the optimal f
   with respect to the cap sigma; it is 0 when the cap is inactive.

   A bound multiplier is 0 (not NaN) for a weight with no finite bound. NaN is only for an absent named
   constraint: the budget constraint does not exist when there is cash, the target-return constraint when no
   target was set (and for the objectives that take none), the volatility constraint unless the objective is
   PO_OBJECTIVE_MAX_RETURN, and the leverage constraint is absent in this version.

   One documented exception, for PO_OBJECTIVE_MAX_RETURN: when the cap equals the minimum volatility of the
   feasible set and the efficient frontier is vertical there, lambda_vol = -inf, and the budget multiplier
   (without cash) and the bound multipliers are NaN even though those constraints exist. */
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

/* For PO_OBJECTIVE_RISK_PARITY, sum_i (c_i - b_i)^2 at the reported point, with c_i = w_i (Sigma w)_i / w' Sigma w
   the risk contributions and b the risk budget; 0 for every other objective; NaN when the weights are NaN. */
PO_API po_status po_result_residual(const po_result* r, double* out);
PO_API void po_result_free(po_result* r);

/* Portfolio metrics for weights w (n must be the weight count of d; cash last). The return is
   sum_i w_i mu_i + w_cash * rate and needs mu; the volatility is sqrt(w' Sigma w) over the risky weights only and
   needs the covariance. A missing input is PO_INVALID_ARGUMENT. */
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

#ifdef __cplusplus
}
#endif

#endif
