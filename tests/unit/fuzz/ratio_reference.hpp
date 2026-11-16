#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "core/data.hpp"
#include "core/objective.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "fuzz/frontier_reference.hpp"
#include "kkt.hpp"
#include "user_kkt.hpp"

namespace po::testing
{

// sigma = sqrt(diag Sigma).
inline std::vector<double> asset_volatilities(const po::Data& data)
{
    std::size_t n = data.asset_count();
    std::vector<double> sigma(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        sigma[i] = std::sqrt(data.covariance()[i * n + i]);
    }
    return sigma;
}

// v of the ratio objectives: mu - rf 1 for MAX_SHARPE (rf = 0 without cash), sigma for MAX_DIVERSIFICATION.
inline std::vector<double> ratio_vector(const po::Data& data, po::Objective objective)
{
    if (objective == po::Objective::MaxSharpe)
    {
        return excess_mu(data);
    }
    return asset_volatilities(data);
}

// v'w / vol(w) over the risky weights.
inline double ratio_of(const po::Data& data, std::span<const double> v, std::span<const double> weights)
{
    return dot(weights.first(data.asset_count()), v) / portfolio_volatility(data, weights);
}

// grad f of f = -v'w / vol(w) over the risky weights, as its two terms -v / vol and (v'w) Sigma w / vol^3.
inline ObjectiveGradient ratio_gradient(const po::Data& data, std::span<const double> v,
                                        const po::Result& result)
{
    std::span<const double> risky = std::span<const double>(result.weights).first(data.asset_count());
    std::vector<double> curvature = covariance_times(data, risky);
    double volatility = std::sqrt(dot(risky, curvature));
    double cube = volatility * volatility * volatility;
    double ratio_return = dot(risky, v);
    std::vector<double> linear(v.begin(), v.end());
    for (std::size_t i = 0; i < linear.size(); ++i)
    {
        linear[i] = -linear[i] / volatility;
        curvature[i] *= ratio_return / cube;
    }
    return sum_of({linear, curvature});
}

// The problem as written: with cash the cash weight is fixed at 0, so the checks hold it to the bounds
// [0, 0], and the risky weights to those of spec.
inline po::Spec with_cash_fixed(po::Spec spec, const po::Data& data)
{
    if (!data.risk_free().has_value())
    {
        return spec;
    }
    std::vector<double> lower;
    std::vector<double> upper;
    for (std::size_t k = 0; k < data.asset_count(); ++k)
    {
        lower.push_back(lower_of(spec, k));
        upper.push_back(upper_of(spec, k));
    }
    lower.push_back(0.0);
    upper.push_back(0.0);
    spec.asset_lower = std::move(lower);
    spec.asset_upper = std::move(upper);
    return spec;
}

struct RatioReference
{
    po::SolveStatus status;
    // The maximum of v'w / vol(w); NaN unless status is Optimal.
    double ratio;
};

// With cash fixed at 0 the ratio problem is the one without cash on the risky assets, 1'w = 1, with mean v:
// the frontier is that of MIN_VARIANCE over this data, the risky bounds and the leverage of spec.
struct RatioFrontier
{
    po::Spec spec;
    po::Data data;
};

inline RatioFrontier risky_frontier(const po::Spec& spec, const po::Data& data, std::span<const double> v)
{
    std::size_t n = data.asset_count();
    RatioFrontier frontier{po::Spec{}, po::Data::create(n, v, data.covariance()).value()};
    frontier.spec.lower = spec.lower;
    frontier.spec.upper = spec.upper;
    if (!spec.asset_lower.empty())
    {
        frontier.spec.asset_lower.assign(spec.asset_lower.begin(), spec.asset_lower.begin() + n);
    }
    if (!spec.asset_upper.empty())
    {
        frontier.spec.asset_upper.assign(spec.asset_upper.begin(), spec.asset_upper.begin() + n);
    }
    frontier.spec.leverage = spec.leverage;
    return frontier;
}

// The best ratio of the frontier points solved so far.
struct RatioSearch
{
    const RatioFrontier& frontier;
    double best;
};

inline double point_ratio(const RatioFrontier& frontier, const po::Result& point)
{
    double volatility = portfolio_volatility(frontier.data, point.weights);
    return portfolio_return(frontier.data, point.weights) / volatility;
}

// The ratio of MIN_VARIANCE(target_return = target), -infinity when it is not Optimal.
inline double ratio_at(RatioSearch& search, double target)
{
    po::Result point = frontier_point(search.frontier.spec, search.frontier.data, target);
    if (point.status != po::SolveStatus::Optimal)
    {
        return -std::numeric_limits<double>::infinity();
    }
    double ratio = point_ratio(search.frontier, point);
    search.best = std::fmax(search.best, ratio);
    return ratio;
}

// d/dr (r / vol(r)) = (2 vol^2 - r lambda_ret) / (2 vol^3), with lambda_ret = d vol^2 / dr the target return
// multiplier. By the multiplier convention w'grad f = 2 vol^2 = lambda_budget + lambda_ret r + beta'w +
// lambda_lev ||w||_1, so the sign is that of lambda_budget + beta'w + lambda_lev ||w||_1: terms of the order
// of r, not a difference of terms of the order of r^2, whose rounding hides an increase of the ratio at a
// large r.
inline bool ratio_increases(const po::Result& point)
{
    double slope = point.budget_multiplier + dot(point.bound_multipliers, point.weights);
    if (std::isfinite(point.leverage_multiplier))
    {
        double gross = 0.0;
        for (double weight : point.weights)
        {
            gross += std::fabs(weight);
        }
        slope += point.leverage_multiplier * gross;
    }
    return slope > 0.0;
}

// The ratio is quasi-concave in r where it is positive (vol(r) is convex), and increasing where it is not, so
// the search keeps the side of the better point, to 1e-13 relative.
inline void golden_section(RatioSearch& search, double low, double high)
{
    constexpr double kRelativeWidth = 1e-13;
    constexpr double kShrink = std::numbers::phi - 1.0;
    constexpr int kSteps = 400;
    double inner_low = high - kShrink * (high - low);
    double inner_high = low + kShrink * (high - low);
    double ratio_low = ratio_at(search, inner_low);
    double ratio_high = ratio_at(search, inner_high);
    for (int step = 0; step < kSteps; ++step)
    {
        if (high - low <= kRelativeWidth * std::fmax(std::fabs(low), std::fabs(high)))
        {
            return;
        }
        if (ratio_low >= ratio_high)
        {
            high = inner_high;
            inner_high = inner_low;
            ratio_high = ratio_low;
            inner_low = high - kShrink * (high - low);
            ratio_low = ratio_at(search, inner_low);
            continue;
        }
        low = inner_low;
        inner_low = inner_high;
        ratio_low = ratio_high;
        inner_high = low + kShrink * (high - low);
        ratio_high = ratio_at(search, inner_high);
    }
}

inline RatioReference positive_ratio(double best)
{
    if (best > 0.0)
    {
        return RatioReference{po::SolveStatus::Optimal, best};
    }
    return RatioReference{po::SolveStatus::Infeasible, std::numeric_limits<double>::quiet_NaN()};
}

// Doublings of the target step beyond the minimum-variance return: a ratio still increasing at the last,
// about 1e9 times the scale of the returns, is taken to increase without bound (docs/guide/validation.md).
inline constexpr int kRatioDoublings = 30;

// The efficient branch from the minimum-variance return r_mv, by doubling r - r_mv until the ratio stops
// increasing (the maximum lies behind) or no portfolio reaches r (the maximum lies before the largest return,
// found by bisection); then golden-section search on the ratio. UNBOUNDED when the ratio still increases
// after kRatioDoublings.
inline RatioReference search_branch(RatioSearch& search, double low, double step)
{
    for (int doubling = 0; doubling < kRatioDoublings; ++doubling)
    {
        double high = low + step;
        po::Result point = frontier_point(search.frontier.spec, search.frontier.data, high);
        if (point.status != po::SolveStatus::Optimal)
        {
            TargetBracket bracket{low, high, low};
            double uncapped = std::numeric_limits<double>::infinity();
            narrow(search.frontier.spec, search.frontier.data, uncapped, bracket);
            ratio_at(search, bracket.low);
            golden_section(search, low, bracket.low);
            return positive_ratio(search.best);
        }
        search.best = std::fmax(search.best, point_ratio(search.frontier, point));
        if (!ratio_increases(point))
        {
            golden_section(search, low, high);
            return positive_ratio(search.best);
        }
        low = high;
        step *= 2.0;
    }
    return RatioReference{po::SolveStatus::Unbounded, std::numeric_limits<double>::quiet_NaN()};
}

// The ratio maximum by another route: golden-section search on the target return r of
// MIN_VARIANCE(target_return = r, same bounds and leverage) over the efficient branch, maximizing r / vol(r)
// on the risky assets with mean v (with cash, v = mu - rf 1 is the excess return, and the cash weight is
// fixed at 0, so 1'w = 1 over the risky weights). INFEASIBLE when no r gives a positive ratio, or the cash
// bounds exclude 0; UNBOUNDED when the ratio increases along the whole branch.
inline RatioReference max_ratio_reference(const po::Spec& spec, const po::Data& data,
                                          std::span<const double> v)
{
    std::size_t n = data.asset_count();
    bool cash_excluded = data.risk_free().has_value() && (lower_of(spec, n) > 0.0 || upper_of(spec, n) < 0.0);
    RatioFrontier frontier = risky_frontier(spec, data, v);
    po::Result minimum = frontier_point(frontier.spec, frontier.data, std::nullopt);
    if (cash_excluded || minimum.status != po::SolveStatus::Optimal)
    {
        return RatioReference{po::SolveStatus::Infeasible, std::numeric_limits<double>::quiet_NaN()};
    }
    RatioSearch search{frontier, point_ratio(frontier, minimum)};
    double low = portfolio_return(frontier.data, minimum.weights);
    return search_branch(search, low, std::max(std::fabs(low), largest_magnitude(v)));
}

} // namespace po::testing
