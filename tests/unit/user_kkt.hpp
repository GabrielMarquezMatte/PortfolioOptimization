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

inline void expect_user_bounds(const po::Spec& spec, const po::Result& result, double tolerance)
{
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
            double row_scale = std::fmax(1.0, std::fmax(std::fabs(lower), std::fabs(weight)));
            EXPECT_LE(complementarity_of_row(beta, weight - lower, row_scale), tolerance);
        }
        if (beta < 0.0)
        {
            double row_scale = std::fmax(1.0, std::fmax(std::fabs(upper), std::fabs(weight)));
            EXPECT_LE(complementarity_of_row(beta, upper - weight, row_scale), tolerance);
        }
    }
}

// The sum of the constraint terms of the convention, each skipped when its multiplier is NaN:
//   budget * 1 + target_return * (mu - rf 1) + volatility * Sigma w / vol + beta_i e_i - beta_0 * 1.
inline std::vector<double> constraint_gradient(const po::Data& data, const po::Result& result)
{
    std::size_t n = data.asset_count();
    std::vector<double> total(result.bound_multipliers.begin(), result.bound_multipliers.begin() + n);
    std::vector<double> return_gradient = excess_mu(data);
    std::vector<double> curvature = covariance_times(data, result.weights);
    double volatility = std::sqrt(dot(std::span<const double>(result.weights).first(n), curvature));
    for (std::size_t i = 0; i < n; ++i)
    {
        if (data.risk_free().has_value())
        {
            total[i] -= result.bound_multipliers[n];
        }
        if (!std::isnan(result.budget_multiplier))
        {
            total[i] += result.budget_multiplier;
        }
        if (!std::isnan(result.target_return_multiplier))
        {
            total[i] += result.target_return_multiplier * return_gradient[i];
        }
        if (!std::isnan(result.volatility_multiplier))
        {
            total[i] += result.volatility_multiplier * curvature[i] / volatility;
        }
    }
    return total;
}

inline void expect_user_stationarity(const po::Data& data, const po::Result& result,
                                     std::span<const double> objective_gradient, double tolerance)
{
    std::size_t n = data.asset_count();
    ASSERT_EQ(objective_gradient.size(), n);
    std::vector<double> constraints = constraint_gradient(data, result);
    std::vector<double> residual(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        residual[i] = objective_gradient[i] - constraints[i];
    }
    double scale = std::fmax(1.0, largest_magnitude(objective_gradient));
    scale = std::fmax(scale, largest_magnitude(result.bound_multipliers));
    EXPECT_LE(largest_magnitude(residual) / scale, tolerance);
}

// Checks, independently of the QP, the equalities, the bounds with per-row complementarity, and
//   grad f(w) = budget * 1 + target_return * (mu - rf 1) + volatility * Sigma w / vol + beta - beta_0 * 1,
// where objective_gradient is grad f over the risky weights.
inline void expect_user_kkt(const po::Spec& spec, const po::Data& data, const po::Result& result,
                            std::span<const double> objective_gradient, double tolerance = 1e-10)
{
    ASSERT_EQ(result.status, po::SolveStatus::Optimal);
    ASSERT_EQ(result.weights.size(), data.weight_count());
    ASSERT_EQ(result.bound_multipliers.size(), data.weight_count());
    EXPECT_EQ(result.residual, 0.0);
    expect_user_equalities(spec, data, result, tolerance);
    expect_user_bounds(spec, result, tolerance);
    expect_user_stationarity(data, result, objective_gradient, tolerance);
}

} // namespace po::testing
