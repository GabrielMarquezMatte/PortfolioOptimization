#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "kkt.hpp"

namespace po::testing
{

inline double dot(std::span<const double> a, std::span<const double> b)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        sum += a[i] * b[i];
    }
    return sum;
}

// Sigma w over the risky weights; a cash entry past the first n weights is ignored.
inline std::vector<double> covariance_times(const po::Data& data, std::span<const double> weights)
{
    std::size_t n = data.asset_count();
    std::vector<double> product(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        product[i] = dot(data.covariance().subspan(i * n, n), weights.first(n));
    }
    return product;
}

// mu - rf 1, with rf = 0 without cash.
inline std::vector<double> excess_mu(const po::Data& data)
{
    double rate = data.risk_free().value_or(0.0);
    std::vector<double> excess(data.mu().begin(), data.mu().end());
    for (double& value : excess)
    {
        value -= rate;
    }
    return excess;
}

inline double lower_of(const po::Spec& spec, std::size_t weight)
{
    if (spec.asset_lower.empty())
    {
        return spec.lower;
    }
    return spec.asset_lower[weight];
}

inline double upper_of(const po::Spec& spec, std::size_t weight)
{
    if (spec.asset_upper.empty())
    {
        return spec.upper;
    }
    return spec.asset_upper[weight];
}

inline void expect_user_equalities(const po::Spec& spec, const po::Data& data, const po::Result& result,
                                   double tolerance)
{
    double total = 0.0;
    for (double weight : result.weights)
    {
        total += weight;
    }
    EXPECT_NEAR(total, 1.0, tolerance);
    EXPECT_EQ(std::isnan(result.budget_multiplier), data.risk_free().has_value());
    EXPECT_EQ(std::isnan(result.target_return_multiplier), !spec.target_return.has_value());
    EXPECT_EQ(std::isnan(result.volatility_multiplier), !spec.max_volatility.has_value());
    EXPECT_TRUE(std::isnan(result.leverage_multiplier));
    if (!spec.target_return.has_value())
    {
        return;
    }
    std::span<const double> risky = std::span<const double>(result.weights).first(data.asset_count());
    double portfolio_return = dot(risky, data.mu());
    if (data.risk_free().has_value())
    {
        portfolio_return += *data.risk_free() * result.weights.back();
    }
    EXPECT_NEAR(portfolio_return, *spec.target_return, tolerance);
}

// Complementarity of a bound row relative to max(|bound|, |w|max), with no floor of 1: a bound at a small scale
// is held to the same relative accuracy as one at scale 1. The scale is that of all the weights, not of the one
// weight, because the rounding of a weight at its bound comes from the whole solution.
inline double bound_complementarity(double beta, double slack, double bound, double weight_scale)
{
    double row_scale = std::fmax(std::fabs(bound), weight_scale);
    if (row_scale == 0.0)
    {
        return std::fabs(slack);
    }
    return complementarity_of_row(beta, slack, row_scale);
}

inline void expect_user_bounds(const po::Spec& spec, const po::Result& result, double tolerance)
{
    double weight_scale = largest_magnitude(result.weights);
    for (std::size_t k = 0; k < result.weights.size(); ++k)
    {
        SCOPED_TRACE(k);
        double weight = result.weights[k];
        double beta = result.bound_multipliers[k];
        double lower = lower_of(spec, k);
        double upper = upper_of(spec, k);
        EXPECT_GE(weight, lower - tolerance);
        EXPECT_LE(weight, upper + tolerance);
        if (beta > 0.0)
        {
            EXPECT_LE(bound_complementarity(beta, weight - lower, lower, weight_scale), tolerance);
        }
        if (beta < 0.0)
        {
            EXPECT_LE(bound_complementarity(beta, upper - weight, upper, weight_scale), tolerance);
        }
    }
}

// grad f over the risky weights as the sum of the terms the objective is made of. The stationarity residual is
// relative to the largest magnitude of a term, so that a gradient that cancels to rounding (an unconstrained
// optimum) is not compared with its own rounding.
struct ObjectiveGradient
{
    std::vector<double> value;
    double scale = 0.0;
};

inline ObjectiveGradient sum_of(const std::vector<std::vector<double>>& terms)
{
    ObjectiveGradient gradient{std::vector<double>(terms.front().size(), 0.0), 0.0};
    for (const std::vector<double>& term : terms)
    {
        for (std::size_t i = 0; i < term.size(); ++i)
        {
            gradient.value[i] += term[i];
        }
        gradient.scale = worst_of(gradient.scale, largest_magnitude(term));
    }
    return gradient;
}

// Sigma w / vol, the gradient of the volatility over the risky weights; empty when the volatility multiplier is
// NaN.
inline std::vector<double> volatility_gradient(const po::Data& data, const po::Result& result)
{
    std::vector<double> gradient;
    if (std::isnan(result.volatility_multiplier))
    {
        return gradient;
    }
    gradient = covariance_times(data, result.weights);
    double volatility = std::sqrt(dot(std::span<const double>(result.weights).first(data.asset_count()), gradient));
    for (double& entry : gradient)
    {
        entry /= volatility;
    }
    return gradient;
}

// The constraint terms of the convention, each skipped when its multiplier is NaN:
//   budget * 1 + target_return * (mu - rf 1) + volatility * Sigma w / vol + beta_i e_i - beta_0 * 1.
struct ConstraintTerms
{
    std::vector<double> total;
    // The largest magnitude of a single term, the scale of the stationarity residual.
    double scale = 0.0;
};

inline ConstraintTerms constraint_terms(const po::Data& data, const po::Result& result)
{
    std::size_t n = data.asset_count();
    std::span<const double> beta = std::span<const double>(result.bound_multipliers).first(n);
    ConstraintTerms terms{std::vector<double>(beta.begin(), beta.end()), largest_magnitude(beta)};
    std::vector<double> return_gradient = excess_mu(data);
    std::vector<double> curvature = volatility_gradient(data, result);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (data.risk_free().has_value())
        {
            terms.total[i] -= result.bound_multipliers[n];
            terms.scale = worst_of(terms.scale, std::fabs(result.bound_multipliers[n]));
        }
        if (!std::isnan(result.budget_multiplier))
        {
            terms.total[i] += result.budget_multiplier;
            terms.scale = worst_of(terms.scale, std::fabs(result.budget_multiplier));
        }
        if (!std::isnan(result.target_return_multiplier))
        {
            double term = result.target_return_multiplier * return_gradient[i];
            terms.total[i] += term;
            terms.scale = worst_of(terms.scale, std::fabs(term));
        }
        if (!curvature.empty())
        {
            double term = result.volatility_multiplier * curvature[i];
            terms.total[i] += term;
            terms.scale = worst_of(terms.scale, std::fabs(term));
        }
    }
    return terms;
}

// The residual relative to the largest magnitude of a term of the gradient or of the constraints, with no floor of
// 1: a problem at a small scale is held to the same relative accuracy as one at scale 1.
inline void expect_user_stationarity(const po::Data& data, const po::Result& result,
                                     const ObjectiveGradient& objective_gradient, double tolerance)
{
    std::size_t n = data.asset_count();
    ASSERT_EQ(objective_gradient.value.size(), n);
    ConstraintTerms terms = constraint_terms(data, result);
    std::vector<double> residual(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        residual[i] = objective_gradient.value[i] - terms.total[i];
    }
    double scale = worst_of(terms.scale, objective_gradient.scale);
    EXPECT_LE(relative_to(largest_magnitude(residual), scale), tolerance);
}

// With a volatility cap s and a finite multiplier: vol <= s (1 + tolerance), lambda_vol <= 0, and
// lambda_vol (vol - s) small relative to |lambda_vol| s.
inline void expect_user_volatility(const po::Spec& spec, const po::Data& data, const po::Result& result,
                                   double tolerance)
{
    if (!spec.max_volatility.has_value() || !std::isfinite(result.volatility_multiplier))
    {
        return;
    }
    double cap = *spec.max_volatility;
    std::span<const double> risky = std::span<const double>(result.weights).first(data.asset_count());
    double volatility = std::sqrt(dot(risky, covariance_times(data, risky)));
    EXPECT_LE(volatility, cap * (1.0 + tolerance));
    EXPECT_LE(result.volatility_multiplier, 0.0);
    EXPECT_LE(complementarity_of_row(result.volatility_multiplier, volatility - cap, cap), tolerance);
}

// Checks, independently of the QP, the equalities, the bounds with per-row complementarity, the volatility cap
// with its complementarity, and
//   grad f(w) = budget * 1 + target_return * (mu - rf 1) + volatility * Sigma w / vol + beta - beta_0 * 1,
// where objective_gradient is grad f over the risky weights, with the scale of its terms.
inline void expect_user_kkt(const po::Spec& spec, const po::Data& data, const po::Result& result,
                            const ObjectiveGradient& objective_gradient, double tolerance = 1e-10)
{
    ASSERT_EQ(result.status, po::SolveStatus::Optimal);
    ASSERT_EQ(result.weights.size(), data.weight_count());
    ASSERT_EQ(result.bound_multipliers.size(), data.weight_count());
    EXPECT_EQ(result.residual, 0.0);
    expect_user_equalities(spec, data, result, tolerance);
    expect_user_bounds(spec, result, tolerance);
    expect_user_volatility(spec, data, result, tolerance);
    expect_user_stationarity(data, result, objective_gradient, tolerance);
}

} // namespace po::testing
