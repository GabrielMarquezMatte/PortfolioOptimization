#pragma once

#include <span>
#include <vector>

#include "core/data.hpp"
#include "core/spec.hpp"
#include "riskparity/objective.hpp"

namespace po::riskparity
{

// min f over the equalities of spec, in this order and only those present: the budget 1'w = 1, the target
// return mu'w = r and the target volatility vol(w) = s, and over the bounds. budget is the risk budget b of
// f.
struct KktProblem
{
    const Spec& spec;
    const Data& data;
    std::span<const double> budget;
};

// A candidate for grad f(w) = sum_k lambda_k grad h_k(w) + beta, the convention of core/result.hpp.
struct KktPoint
{
    std::vector<double> weights;
    // One per equality, in the order of KktProblem.
    std::vector<double> equality_multipliers;
    // Lower minus upper. A positive entry names the lower bound of its weight as active, a negative one the
    // upper bound.
    std::vector<double> bound_multipliers;
};

// The KKT system at a point. Every residual comes with its noise: how far it moves when each weight is
// rounded once, |dw_j| <= eps |w_j|, and each term of it is rounded once.
struct KktEvaluation
{
    // f, grad f, J and the exact Hessian of f.
    Linearization linearization;
    // hess f - volatility multiplier hess vol (n x n, row-major).
    std::vector<double> lagrangian_hessian;
    // grad h_k, one row per equality (m x n, row-major): 1, mu, Sigma w / vol(w).
    std::vector<double> rows;
    // h_k(w): 1'w - 1, mu'w - r, vol(w) - s.
    std::vector<double> gaps;
    std::vector<double> gap_noise;
    // grad f - sum_k lambda_k grad h_k - beta.
    std::vector<double> stationarity;
    std::vector<double> stationarity_noise;
};

KktEvaluation evaluate_kkt(const KktProblem& problem, const KktPoint& point);

// Each the largest ratio of a residual to its noise, so that a value of order 1 is a point that no rounding
// of the weights distinguishes from a KKT point. bounds is the violation of a bound; complementarity, for a
// nonzero beta_j, the smaller of |beta_j| and the distance of w_j to the bound that the sign of beta_j names.
// A NaN anywhere in the point makes its ratio NaN.
struct KktResidual
{
    double stationarity;
    double equalities;
    double bounds;
    double complementarity;
};

KktResidual residual_of(const KktProblem& problem, const KktPoint& point, const KktEvaluation& evaluation);

// The largest of the four; NaN when any is.
double largest(const KktResidual& residual);

// A point is a KKT point when its residual, feasibility included, is within this many times its noise
// (docs/guide/validation.md).
inline constexpr double kKktNoiseFactor = 16.0;

// A residual within its noise says nothing once the noise itself is large: weights that run off along an
// unbounded direction of the polyhedron, where f tends to a limit and its gradient to 0, satisfy every row to
// their rounding. So a KKT point must also resolve the budget row: the rounding of 1'w, eps sum_i |w_i|, is
// at most the feasibility tolerance of the QP solver. factor replaces kKktNoiseFactor for a tighter stopping
// rule inside an iteration.
bool is_kkt_point(const KktProblem& problem, const KktPoint& point, double factor = kKktNoiseFactor);

} // namespace po::riskparity
