#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <span>

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "kkt.hpp"
#include "solve/solve.hpp"
#include "user_kkt.hpp"

namespace po::testing
{

inline double portfolio_volatility(const po::Data& data, std::span<const double> weights)
{
    std::span<const double> risky = weights.first(data.asset_count());
    return std::sqrt(dot(risky, covariance_times(data, risky)));
}

// mu'w over the risky weights plus rf times the cash weight.
inline double portfolio_return(const po::Data& data, std::span<const double> weights)
{
    double total = dot(weights.first(data.asset_count()), data.mu());
    if (data.risk_free().has_value())
    {
        total += *data.risk_free() * weights.back();
    }
    return total;
}

// MIN_VARIANCE over the bounds and the leverage of spec, with a target return when target is given.
inline po::Result frontier_point(const po::Spec& spec, const po::Data& data, std::optional<double> target)
{
    po::Spec frontier;
    frontier.lower = spec.lower;
    frontier.upper = spec.upper;
    frontier.weight_lower = spec.weight_lower;
    frontier.weight_upper = spec.weight_upper;
    frontier.leverage = spec.leverage;
    frontier.target_return = target;
    return po::solve(frontier, data).value();
}

// The return of MIN_VARIANCE(target_return = target) when it is feasible with a volatility of at most cap. Above
// the minimum-variance return this holds on an interval of targets, since the frontier volatility increases with
// the target. The return of the portfolio, not the target, is kept: the solver accepts a target return row to
// its feasibility tolerance, which has a floor of 1e-12 and would otherwise move a small reference.
inline std::optional<double> return_within(const po::Spec& spec, const po::Data& data, double cap, double target)
{
    po::Result point = frontier_point(spec, data, target);
    if (point.status != po::SolveStatus::Optimal || portfolio_volatility(data, point.weights) > cap)
    {
        return std::nullopt;
    }
    return portfolio_return(data, point.weights);
}

// Targets known to be reached (low, whose portfolio return is best) and not reached (high).
struct TargetBracket
{
    double low;
    double high;
    double best;
};

inline constexpr int kReferenceSteps = 400;

inline void widen(const po::Spec& spec, const po::Data& data, double cap, TargetBracket& bracket)
{
    double step = bracket.high - bracket.low;
    for (int doubling = 0; doubling < kReferenceSteps; ++doubling)
    {
        std::optional<double> reached = return_within(spec, data, cap, bracket.high);
        if (!reached)
        {
            return;
        }
        bracket.low = bracket.high;
        bracket.best = *reached;
        step *= 2.0;
        bracket.high = bracket.low + step;
    }
}

inline void narrow(const po::Spec& spec, const po::Data& data, double cap, TargetBracket& bracket)
{
    constexpr double kRelativeWidth = 1e-13;
    for (int bisection = 0; bisection < kReferenceSteps; ++bisection)
    {
        double middle = 0.5 * (bracket.low + bracket.high);
        double width = kRelativeWidth * std::max(std::fabs(bracket.low), std::fabs(bracket.high));
        if (bracket.high - bracket.low <= width || middle == bracket.low || middle == bracket.high)
        {
            return;
        }
        std::optional<double> reached = return_within(spec, data, cap, middle);
        if (!reached)
        {
            bracket.high = middle;
            continue;
        }
        bracket.low = middle;
        bracket.best = *reached;
    }
}

// MAX_RETURN by another route: the largest target return r with MIN_VARIANCE(target_return = r) feasible and of
// volatility <= cap, by bisection on r to 1e-13 relative, as the return of that portfolio; NaN when no feasible
// portfolio has volatility <= cap. The bounds and the leverage are those of spec; the frontier volatility is still
// increasing above the minimum-variance return, since the leverage is a convex constraint.
inline double max_return_reference(const po::Spec& spec, const po::Data& data, double cap)
{
    po::Result minimum = frontier_point(spec, data, std::nullopt);
    if (minimum.status != po::SolveStatus::Optimal || portfolio_volatility(data, minimum.weights) > cap)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    double low = portfolio_return(data, minimum.weights);
    double step = std::max(std::fabs(low), largest_magnitude(data.mu()));
    TargetBracket bracket{low, low + step, low};
    widen(spec, data, cap, bracket);
    narrow(spec, data, cap, bracket);
    return bracket.best;
}

} // namespace po::testing
