#pragma once

#include <cstdint>
#include <limits>
#include <vector>

namespace po
{

enum class SolveStatus : std::int32_t
{
    Optimal = 0,
    Infeasible = 1,
    Unbounded = 2,
    MaxIterations = 3
};

// For the objective f in minimization form (MIN_VARIANCE: f = w'Sigma w, without 1/2; a maximized objective
// enters negated), at the optimum over the risky weights w:
//   grad f(w) = budget * 1 + target_return * grad(return) + volatility * Sigma w / vol(w)
//               + sum_i beta_i e_i - beta_0 * 1,
// with vol(w) = sqrt(w'Sigma w) (the volatility term only for MAX_RETURN), beta = bound_multipliers (lower
// minus upper: >= 0 at an active lower bound, <= 0 at an active upper bound) and beta_0 its cash entry,
// because the cash weight is 1 - 1'w. A multiplier of an absent constraint is NaN. One exception, for
// MAX_RETURN: at a cap equal to the minimum volatility where the efficient frontier is vertical,
// volatility = -inf and the budget multiplier (without cash) and the bound multipliers are NaN even though
// those constraints exist.
struct Result
{
    SolveStatus status = SolveStatus::Infeasible;
    std::vector<double> weights;
    double budget_multiplier = std::numeric_limits<double>::quiet_NaN();
    double target_return_multiplier = std::numeric_limits<double>::quiet_NaN();
    double volatility_multiplier = std::numeric_limits<double>::quiet_NaN();
    double leverage_multiplier = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> bound_multipliers;
    double residual = 0.0;
    int iterations = 0;
};

} // namespace po
