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

// min f(w) = sum_i (c_i - b_i)^2 over the budget row, the target return, the target volatility and the bounds
// of spec, c the risk contributions and b the risk budget. In this order:
// - The positive-orthant solution of Spinu (f = 0). Without a target and inside the box it is the answer,
//   with every multiplier 0.
// - The minimum-variance portfolio over the linear rows and the bounds. Without one the problem is
//   Infeasible, and so is a target volatility s below its volatility sigma_min. At s = sigma_min (to
//   problems::kMinimumVolatilityTolerance) that portfolio is the only feasible point: Optimal, with a
//   volatility multiplier of -infinity and every other multiplier NaN, because there the gradients of the
//   constraints are dependent and in general no multipliers exist.
// - From the point of the polyhedron nearest to the Spinu solution, the Gauss-Newton descent over the
//   polyhedron, inside an augmented Lagrangian when there is a target volatility, and then Newton on the KKT
//   system.
// Optimal is a point whose KKT residual, feasibility included, is within the noise of the problem. When the
// iteration does not end at one (a target return on an unbounded line where f has no minimizer, a target
// volatility out of reach) the status is MaxIterations, with NaN weights, multipliers and residual. A local
// method: Optimal means a stationary point reached from the Spinu solution, not the global minimum.
// iterations counts the Newton steps of Spinu, the steps of every descent and the Newton steps on the KKT
// system.
std::expected<Result, Failure> solve_risk_parity(const Spec& spec, const Data& data);

} // namespace po::riskparity
