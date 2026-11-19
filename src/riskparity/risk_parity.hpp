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
// - One asset without a target volatility: w = 1, the only feasible point, with every multiplier 0.
// - From the point of the polyhedron nearest to the Spinu solution, the Gauss-Newton descent over the
//   polyhedron, inside an augmented Lagrangian when there is a target volatility, and then Newton on the KKT
//   system.
// Optimal is a point whose KKT residual, feasibility included, is within the noise of the problem. When the
// iteration does not end at one (a target return on an unbounded line where f has no minimizer, a target
// volatility out of reach) the status is MaxIterations, with NaN weights, multipliers and residual. So is a
// point whose budget row 1'w = 1 is not resolved by its own rounding, eps (sum_i |w_i| + 1) > 1e-12: a gross
// exposure sum_i |w_i| above about 4.5e3 is MaxIterations even when it is stationary.
// With a target volatility and few or no active bounds at cond(Sigma) around 1e8 the weights are determined
// only to the rounding of the KKT residual amplified by the conditioning of the problem: an Optimal result
// can be up to about 1e-4 relative from the exact KKT point (measured up to 2.8e-5), while its objective
// value is right to 1e-7 (1 + f). The multipliers are those of one Newton solve at the returned weights,
// kept when its safeguard accepts it (least squares at the weights otherwise). With a target volatility at
// least 1e-3 relative above sigma_min they are measured, not derived, within 8.1e-13 relative at
// cond(Sigma) = 1e4, 2.4e-12 at 1e6 and 9.3e-11 at 1e8 (the largest measured, n = 5 to 20, no box; 1.6e-10 for
// the bound multipliers at 1e8), at the weights those solves returned; further off, the multipliers follow the
// weights. Nearer sigma_min the error grows roughly like 1 / delta, s = sigma_min (1 + delta): measured
// up to 3.2e-12 at delta = 1e-4 and 4.3e-7 at 1e-10 (cond(Sigma) = 1e4), 6.2e-7 at 1e-10 over
// cond(Sigma) <= 1e8. Above 1e8 they are not validated (measured up to 2.1e-7 off at 1e10 on 8 shapes, and 0.25 on
// a wider grid before the portable generators;
// docs/guide/validation.md). A local method: Optimal means a stationary point reached
// from the Spinu solution, not the global minimum. iterations counts the Newton steps of Spinu, the steps of
// every descent and the Newton steps on the KKT system.
std::expected<Result, Failure> solve_risk_parity(const Spec& spec, const Data& data);

} // namespace po::riskparity
