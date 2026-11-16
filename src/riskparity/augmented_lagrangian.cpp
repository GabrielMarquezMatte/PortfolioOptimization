#include "riskparity/augmented_lagrangian.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include "riskparity/objective.hpp"
#include "riskparity/polyhedral_descent.hpp"

namespace po::riskparity
{
namespace
{

constexpr double kInitialPenalty = 10.0;
constexpr double kPenaltyGrowth = 10.0;
constexpr double kRequiredDecrease = 0.5;
constexpr double kFirstTolerance = 1e-2;
constexpr double kToleranceDecrease = 0.1;

// lambda s and rho s^2, the multiplier and the penalty of the relative gap vol(w) / s - 1, with the
// tolerance of the next subproblem and the size of the last gap.
struct Outer
{
    double multiplier = 0.0;
    double penalty = kInitialPenalty;
    double tolerance;
    double previous_gap = std::numeric_limits<double>::infinity();
    // With the penalty at its ceiling: the smallest gap seen there and the outer iterations since it halved.
    double ceiling_gap = std::numeric_limits<double>::infinity();
    int stalled = 0;
};

double relative_gap(const Data& data, std::span<const double> weights, double target)
{
    double variance = compensated_dot(weights, covariance_times(data, weights));
    return std::sqrt(variance) / target - 1.0;
}

bool solved(DescentStop stop)
{
    return stop != DescentStop::NotPositiveDefinite && stop != DescentStop::QpFailed;
}

// The subproblem of the outer iterate, from the last point. Its multipliers with (lambda - rho h3), here
// (multiplier - penalty gap) / s, are those of the convention.
DescentResult subproblem(const KktProblem& problem, const Outer& outer, std::span<const double> start,
                         int max_iterations)
{
    double target = *problem.spec.target_volatility;
    VolatilityTerm term{outer.multiplier / target, outer.penalty / (target * target), target};
    DescentProblem descent{problem.spec, problem.data, problem.budget, term};
    return descend(descent, start, DescentOptions{kDescentProximal, outer.tolerance, max_iterations});
}

void record(const DescentResult& descent, double volatility_multiplier, KktPoint& point)
{
    point.weights = descent.weights;
    point.equality_multipliers = descent.equality_multipliers;
    point.equality_multipliers.push_back(volatility_multiplier);
    point.bound_multipliers = descent.bound_multipliers;
}

// Andreani, Birgin, Martinez and Schuverdt (2007), Algorithm 3.1, steps 3 and 4: the multiplier estimate is
// projected onto its interval, the penalty grows when the gap did not fall to the required fraction of the
// previous one, and the next subproblem is solved to a tighter tolerance.
void update(Outer& outer, double estimate, double gap, double final_tolerance)
{
    if (std::fabs(gap) > kRequiredDecrease * outer.previous_gap)
    {
        outer.penalty = std::fmin(kPenaltyGrowth * outer.penalty, kAugmentedPenaltyCeiling);
    }
    outer.previous_gap = std::fabs(gap);
    outer.multiplier = std::clamp(estimate, -kAugmentedMultiplierBound, kAugmentedMultiplierBound);
    outer.tolerance = std::fmax(final_tolerance, kToleranceDecrease * outer.tolerance);
}

// With the penalty at its ceiling the gap can only fall by the multiplier; when it does not halve for
// kAugmentedStallIterations outer iterations in a row, s is out of reach of the polyhedron.
bool stalls(Outer& outer, double gap)
{
    if (outer.penalty < kAugmentedPenaltyCeiling)
    {
        return false;
    }
    if (std::fabs(gap) <= kRequiredDecrease * outer.ceiling_gap)
    {
        outer.ceiling_gap = std::fabs(gap);
        outer.stalled = 0;
        return false;
    }
    ++outer.stalled;
    return outer.stalled >= kAugmentedStallIterations;
}

} // namespace

AugmentedResult augmented_lagrangian(const KktProblem& problem, std::span<const double> start,
                                     const AugmentedOptions& options)
{
    double target = *problem.spec.target_volatility;
    Outer outer;
    outer.tolerance = std::fmax(options.tolerance, kFirstTolerance);
    KktPoint first{std::vector<double>(start.begin(), start.end()), {}, {}};
    AugmentedResult result{std::move(first), outer.penalty, 0, 0, false};
    while (result.outer_iterations < options.max_outer_iterations)
    {
        DescentResult descent =
            subproblem(problem, outer, result.point.weights, options.max_inner_iterations);
        ++result.outer_iterations;
        result.inner_iterations += descent.iterations;
        result.penalty = outer.penalty;
        double gap = relative_gap(problem.data, descent.weights, target);
        double estimate = outer.multiplier - outer.penalty * gap;
        record(descent, estimate / target, result.point);
        if (!solved(descent.stop))
        {
            break;
        }
        double largest = std::fmax(descent.stationarity, std::fabs(gap));
        if (largest <= options.tolerance || is_kkt_point(problem, result.point, kAugmentedStopFactor))
        {
            result.converged = true;
            break;
        }
        update(outer, estimate, gap, options.tolerance);
        if (stalls(outer, gap))
        {
            break;
        }
    }
    return result;
}

} // namespace po::riskparity
