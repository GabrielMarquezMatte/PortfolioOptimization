#pragma once

#include <optional>
#include <span>
#include <vector>

#include "core/data.hpp"
#include "core/spec.hpp"

namespace po::riskparity
{

// The augmented-Lagrangian term of the volatility row h3 = vol(w) - target:
// L3 = -multiplier h3 + penalty / 2 h3^2.
struct VolatilityTerm
{
    double multiplier;
    double penalty;
    double target;
};

// L = f + L3 (L3 absent without a volatility term) over the polyhedron of spec: the budget row, the
// target_return row and the bounds. budget is the risk budget b of f.
struct DescentProblem
{
    const Spec& spec;
    const Data& data;
    std::span<const double> budget;
    std::optional<VolatilityTerm> volatility;
};

// proximal sets the first proximal weight, tau_0 = proximal tr(2 J'J) / n at the start, positive; tolerance
// bounds the relative KKT residual.
struct DescentOptions
{
    double proximal;
    double tolerance;
    int max_iterations;
};

// The first proximal weight of every descent of RISK_PARITY (docs/guide/validation.md).
inline constexpr double kDescentProximal = 0.05;

enum class DescentStop
{
    // The relative KKT residual is at most the tolerance.
    Converged,
    IterationCap,
    // No step of the segment decreases L enough: the rounding floor of L.
    NoDecrease,
    // Q could not be factored after the bounded growth of tau; weights is the last iterate, the multipliers
    // and the residual are NaN.
    NotPositiveDefinite,
    // The QP had no solution, which a feasible start rules out in exact arithmetic; multipliers and residual
    // are NaN.
    QpFailed
};

// The multipliers are those of the QP solved at weights. At a fixed point they satisfy
//   grad L(w) = budget 1 + target_return mu + sum_i beta_i e_i,
// beta lower minus upper as in core/result.hpp, so that with grad L = grad f + (penalty h3 - multiplier)
// grad h3 they are the user's multipliers with volatility multiplier (multiplier - penalty h3).
struct DescentResult
{
    std::vector<double> weights;
    // budget, then target_return, those present.
    std::vector<double> equality_multipliers;
    std::vector<double> bound_multipliers;
    // The relative KKT residual at weights: stationarity and complementarity.
    double stationarity;
    int iterations;
    DescentStop stop;
};

// Minimizes L from a feasible start by proximal Gauss-Newton steps: each step solves the convex QP
//   min 1/2 d'Qd + grad L'd over the polyhedron shifted to w, Q = 2 J'J + tau I (+ penalty grad h3 grad h3'),
// and moves along the segment to its solution by Armijo backtracking, so every iterate stays in the
// polyhedron and L never increases (up to the rounding of the iterate across the rows). tau follows
// Levenberg-Marquardt: divided by 3 after a full step, multiplied by 4 after a shortened one, never below
// 100 n eps max_i (2 J'J)_ii, and multiplied by 10 while Q does not factor. Without a volatility term
// everything depends on Sigma only through J, so the iterates do not change when Sigma is scaled. A local
// method: a fixed point is a stationary point of L.
DescentResult descend(const DescentProblem& problem, std::span<const double> start,
                      const DescentOptions& options);

} // namespace po::riskparity
