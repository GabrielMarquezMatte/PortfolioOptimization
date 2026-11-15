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

/* Handles are opaque. A null handle is PO_INVALID_HANDLE, a null output pointer or null input array is
   PO_INVALID_ARGUMENT. After a failed po_*_create or po_solve, *out is NULL (when out itself is not null).
   Every *_free accepts NULL and does nothing. Handles are independent: freeing one never invalidates another.
   Functions that take const handles only read them, so they are safe to call in parallel on the same handles;
   a handle that is being modified (a setter) must not be used by another thread at the same time. */

typedef struct po_spec po_spec;
typedef struct po_data po_data;
typedef struct po_result po_result;

/* Objectives. Only PO_OBJECTIVE_MIN_VARIANCE can be solved in this version: po_solve returns PO_UNSUPPORTED for
   the others. Values only grow, by addition. */
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

PO_API void po_spec_free(po_spec* s);

/* Data. The inputs are copied. mu (n entries) and cov (n x n, row-major) may each be NULL when the objective does
   not need them; po_solve and the portfolio functions report a missing one as PO_INVALID_ARGUMENT. n = 0, a
   non-finite entry, an asymmetric covariance or one that is not positive definite is PO_INVALID_ARGUMENT, with the
   reason in po_last_error(). */
PO_API po_status po_data_create(size_t n, const double* mu, const double* cov, po_data** out);

/* Adds a risk-free asset (cash) as the last weight: the weight count becomes the asset count plus one. Weights
   always sum to one, cash included. rate must be finite. */
PO_API po_status po_data_set_risk_free(po_data* d, double rate);

PO_API po_status po_data_asset_count(const po_data* d, size_t* out);

/* The asset count, plus one when a risk-free rate was set. */
PO_API po_status po_data_weight_count(const po_data* d, size_t* out);

PO_API void po_data_free(po_data* d);

/* Solves the spec on the data; both are copied, so the result does not depend on them afterwards. Returns
   PO_UNSUPPORTED for an objective other than PO_OBJECTIVE_MIN_VARIANCE, PO_INVALID_ARGUMENT when the data lacks
   what the spec needs (the covariance; mu for a target return) or the per-asset bounds do not match the weight
   count. Infeasibility is reported in the result. Safe in parallel with the same read-only spec and data. */
PO_API po_status po_solve(const po_spec* s, const po_data* d, po_result** out);

/* Result, immutable. For a status other than PO_SOLVE_OPTIMAL the weights, bound multipliers and the multipliers
   of the named constraints are NaN.

   Multiplier convention. Let f be the objective in minimization form: for PO_OBJECTIVE_MIN_VARIANCE
   f(w) = w' Sigma w, without the factor 1/2, and a maximized objective enters negated. At the optimum, over the
   risky weights w (the first po_data_asset_count entries),
       grad f(w) = lambda_budget * 1 + lambda_return * grad(return)(w) + sum_i beta_i e_i - beta_0 * 1
   where lambda_budget and lambda_return are po_result_multiplier for PO_CONSTRAINT_BUDGET and
   PO_CONSTRAINT_TARGET_RETURN, beta is po_result_bound_multipliers (lower-bound multiplier minus upper-bound
   multiplier: >= 0 at an active lower bound, <= 0 at an active upper bound) and beta_0 is its cash entry, which
   enters with a minus sign because the cash weight is 1 - 1'w.

   A bound multiplier is 0 (not NaN) for a weight with no finite bound. NaN is only for an absent named
   constraint: the budget constraint does not exist when there is cash, the target-return constraint when no
   target was set, and the volatility and leverage constraints are absent in this version. */
PO_API po_status po_result_status(const po_result* r, int32_t* out);

/* The weight count of the data the result was solved on. */
PO_API po_status po_result_weight_count(const po_result* r, size_t* out);

/* The vector accessors require n to be exactly the weight count, otherwise PO_INVALID_ARGUMENT. */
PO_API po_status po_result_weights(const po_result* r, double* out, size_t n);
PO_API po_status po_result_iterations(const po_result* r, int32_t* out);

/* constraint is one of PO_CONSTRAINT_*, otherwise PO_INVALID_ARGUMENT. */
PO_API po_status po_result_multiplier(const po_result* r, int32_t constraint, double* out);
PO_API po_status po_result_bound_multipliers(const po_result* r, double* out, size_t n);
PO_API po_status po_result_residual(const po_result* r, double* out);
PO_API void po_result_free(po_result* r);

/* Portfolio metrics for weights w (n must be the weight count of d; cash last). The return is
   sum_i w_i mu_i + w_cash * rate and needs mu; the volatility is sqrt(w' Sigma w) over the risky weights only and
   needs the covariance. A missing input is PO_INVALID_ARGUMENT. */
PO_API po_status po_portfolio_return(const po_data* d, const double* w, size_t n, double* out);
PO_API po_status po_portfolio_volatility(const po_data* d, const double* w, size_t n, double* out);

#ifdef __cplusplus
}
#endif

#endif
