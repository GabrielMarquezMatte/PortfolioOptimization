#pragma once

#include <expected>
#include <vector>

#include "core/data.hpp"
#include "core/failure.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"

namespace po::riskparity
{

// spec.risk_budget, or 1 / n when it is empty.
std::vector<double> risk_budget_of(const Spec& spec, const Data& data);

// min sum_i (c_i - b_i)^2 over the budget row, c_i the risk contributions, started from the positive-orthant
// solution of Spinu. Without target_return, target_volatility or an active bound the minimizer is that solution:
// optimal, with every multiplier 0 and the residual at the rounding floor. A local method: Optimal means a
// stationary point, not the global minimum.
std::expected<Result, Failure> solve_risk_parity(const Spec& spec, const Data& data);

} // namespace po::riskparity
