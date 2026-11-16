#pragma once

#include <span>

#include "riskparity/kkt.hpp"

namespace po::riskparity
{

// The interval [-bound, bound] of lambda s and the largest rho s^2 (docs/guide/validation.md).
inline constexpr double kAugmentedMultiplierBound = 1e8;
inline constexpr double kAugmentedPenaltyCeiling = 1e12;
// An iterate ends the iteration when it is a KKT point to this many times its noise: tighter than the gate
// kKktNoiseFactor of the answer, so that the polish starts with a margin (docs/guide/validation.md).
inline constexpr double kAugmentedStopFactor = 2.0;
// With rho s^2 at its ceiling, this many outer iterations in a row without the gap falling to half of the
// smallest one seen there end the iteration, not converged: s is out of reach (docs/guide/validation.md).
inline constexpr int kAugmentedStallIterations = 3;

// tolerance bounds the relative KKT residual of the last descent and the relative gap |vol(w) / s - 1|; an
// iterate that is a KKT point to kAugmentedStopFactor times the noise of the problem (is_kkt_point) ends the
// iteration too, whatever the tolerance.
struct AugmentedOptions
{
    double tolerance;
    int max_outer_iterations;
    int max_inner_iterations;
};

struct AugmentedResult
{
    // The last iterate, with the multipliers of its descent and the volatility multiplier
    // (lambda - rho h3) of the convention; NaN multipliers when that descent could not be solved.
    KktPoint point;
    // rho s^2 at the last outer iteration: the penalty of the relative gap vol(w) / s - 1.
    double penalty;
    int outer_iterations;
    // The iterations of every descent.
    int inner_iterations;
    bool converged;
};

// The augmented Lagrangian of Andreani, Birgin, Martinez and Schuverdt (2007), Algorithm 3.1, for the one row
// h3 = vol(w) - s, the target volatility of spec; the budget, the target return and the bounds stay in the
// polyhedron of every subproblem, which descend minimizes:
//   L(w) = f(w) - lambda h3 + rho / 2 h3^2,  then  lambda <- lambda - rho h3.
// rho is multiplied by 10 when |h3| does not fall to half of its previous value, lambda is kept in a bounded
// interval, and the tolerance of the subproblem is divided by 10 at every outer iteration down to tolerance.
// lambda s and rho s^2 are the numbers that are bounded and grown, so that the iterates do not depend on the
// scale of Sigma. Not converged when the outer iterations run out or the gap stalls with the penalty at its
// ceiling, as they do when s is out of reach.
AugmentedResult augmented_lagrangian(const KktProblem& problem, std::span<const double> start,
                                     const AugmentedOptions& options);

} // namespace po::riskparity
