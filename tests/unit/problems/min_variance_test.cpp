#include "problems/min_variance.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>
#include <span>
#include <vector>

#include "kkt.hpp"
#include "linalg/cholesky.hpp"
#include "problems/solve.hpp"
#include "qp/goldfarb_idnani.hpp"
#include "random_covariance.hpp"

namespace
{

using po::Data;
using po::Result;
using po::SolveStatus;
using po::Spec;

constexpr double kInfinity = std::numeric_limits<double>::infinity();
constexpr double kRiskFree = 0.02;

Data random_data(std::size_t n, std::uint64_t seed, double condition = 100.0)
{
    std::vector<double> mu = po::testing::random_vector(n, seed + 500, 0.03, 0.15);
    return Data::create(n, mu, po::testing::random_covariance(n, seed, condition)).value();
}

std::vector<double> covariance_solve(const Data& data, std::vector<double> rhs)
{
    std::size_t n = data.asset_count();
    std::vector<double> factor = po::linalg::cholesky(data.covariance(), n).value();
    po::linalg::solve_lower(factor, n, rhs);
    po::linalg::solve_lower_transposed(factor, n, rhs);
    return rhs;
}

double dot(std::span<const double> a, std::span<const double> b)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        sum += a[i] * b[i];
    }
    return sum;
}

std::vector<double> covariance_times(const Data& data, std::span<const double> weights)
{
    std::size_t n = data.asset_count();
    std::vector<double> product(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        product[i] = dot(data.covariance().subspan(i * n, n), weights.first(n));
    }
    return product;
}

double variance(const Data& data, const Result& result)
{
    std::span<const double> risky = std::span<const double>(result.weights).first(data.asset_count());
    return dot(risky, covariance_times(data, risky));
}

std::vector<double> excess_mu(const Data& data)
{
    double rate = data.risk_free().value_or(0.0);
    std::vector<double> excess(data.mu().begin(), data.mu().end());
    for (double& value : excess)
    {
        value -= rate;
    }
    return excess;
}

double lower_of(const Spec& spec, std::size_t weight)
{
    if (spec.asset_lower.empty())
    {
        return spec.lower;
    }
    return spec.asset_lower[weight];
}

double upper_of(const Spec& spec, std::size_t weight)
{
    if (spec.asset_upper.empty())
    {
        return spec.upper;
    }
    return spec.asset_upper[weight];
}

void expect_equalities(const Spec& spec, const Data& data, const Result& result, double tolerance)
{
    double total = 0.0;
    for (double weight : result.weights)
    {
        total += weight;
    }
    EXPECT_NEAR(total, 1.0, tolerance);
    EXPECT_EQ(std::isnan(result.budget_multiplier), data.risk_free().has_value());
    EXPECT_EQ(std::isnan(result.target_return_multiplier), !spec.target_return.has_value());
    EXPECT_TRUE(std::isnan(result.volatility_multiplier));
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

void expect_bounds(const Spec& spec, const Result& result, double tolerance)
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
            EXPECT_LE(po::testing::complementarity_of_row(beta, weight - lower, row_scale), tolerance);
        }
        if (beta < 0.0)
        {
            double row_scale = std::fmax(1.0, std::fmax(std::fabs(upper), std::fabs(weight)));
            EXPECT_LE(po::testing::complementarity_of_row(beta, upper - weight, row_scale), tolerance);
        }
    }
}

// The convention under test: grad(w'Sigma w) = budget * 1 + target * grad(return) + beta_i e_i - beta_0 * 1.
void expect_stationarity(const Spec& spec, const Data& data, const Result& result, double tolerance)
{
    std::size_t n = data.asset_count();
    std::vector<double> gradient = covariance_times(data, result.weights);
    std::vector<double> return_gradient = excess_mu(data);
    std::vector<double> residual(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        gradient[i] *= 2.0;
        residual[i] = gradient[i] - result.bound_multipliers[i];
        if (!data.risk_free().has_value())
        {
            residual[i] -= result.budget_multiplier;
        }
        if (data.risk_free().has_value())
        {
            residual[i] += result.bound_multipliers[n];
        }
        if (spec.target_return.has_value())
        {
            residual[i] -= result.target_return_multiplier * return_gradient[i];
        }
    }
    double scale = std::fmax(1.0, po::testing::largest_magnitude(gradient));
    scale = std::fmax(scale, po::testing::largest_magnitude(result.bound_multipliers));
    EXPECT_LE(po::testing::largest_magnitude(residual) / scale, tolerance);
}

void expect_user_kkt(const Spec& spec, const Data& data, const Result& result, double tolerance = 1e-10)
{
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    ASSERT_EQ(result.weights.size(), data.weight_count());
    ASSERT_EQ(result.bound_multipliers.size(), data.weight_count());
    EXPECT_EQ(result.residual, 0.0);
    expect_equalities(spec, data, result, tolerance);
    expect_bounds(spec, result, tolerance);
    expect_stationarity(spec, data, result, tolerance);
}

Result solve_checked(const Spec& spec, const Data& data, double tolerance = 1e-10)
{
    po::problems::Translation translation = po::problems::min_variance_problem(spec, data);
    po::testing::expect_kkt(translation.problem, po::qp::solve(translation.problem), tolerance);
    Result result = po::solve(spec, data).value();
    expect_user_kkt(spec, data, result, tolerance);
    return result;
}

void expect_weights_near(const std::vector<double>& actual, const std::vector<double>& expected, double relative)
{
    ASSERT_GE(actual.size(), expected.size());
    double tolerance = relative * po::testing::largest_magnitude(expected);
    for (std::size_t i = 0; i < expected.size(); ++i)
    {
        EXPECT_NEAR(actual[i], expected[i], tolerance);
    }
}

void expect_no_point(const Result& result, std::size_t weight_count)
{
    EXPECT_EQ(result.status, SolveStatus::Infeasible);
    ASSERT_EQ(result.weights.size(), weight_count);
    ASSERT_EQ(result.bound_multipliers.size(), weight_count);
    for (std::size_t k = 0; k < weight_count; ++k)
    {
        EXPECT_TRUE(std::isnan(result.weights[k]));
        EXPECT_TRUE(std::isnan(result.bound_multipliers[k]));
    }
    EXPECT_TRUE(std::isnan(result.budget_multiplier));
}

Spec long_only()
{
    Spec spec;
    spec.lower = 0.0;
    spec.upper = 1.0;
    return spec;
}

TEST(MinVariance, UnconstrainedMatchesClosedForm)
{
    constexpr std::size_t n = 6;
    Data data = random_data(n, 3);
    std::vector<double> direction = covariance_solve(data, std::vector<double>(n, 1.0));
    double a = dot(std::vector<double>(n, 1.0), direction);
    std::vector<double> expected(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        expected[i] = direction[i] / a;
    }
    Result result = solve_checked(Spec{}, data);
    expect_weights_near(result.weights, expected, 1e-12);
    EXPECT_NEAR(result.budget_multiplier, 2.0 / a, 1e-12 * 2.0 / a);
}

TEST(MinVariance, FrontierIsParabola)
{
    constexpr std::size_t n = 8;
    Data data = random_data(n, 11);
    std::vector<double> ones(n, 1.0);
    std::vector<double> inverse_ones = covariance_solve(data, ones);
    std::vector<double> inverse_mu = covariance_solve(data, std::vector<double>(data.mu().begin(), data.mu().end()));
    double a = dot(ones, inverse_ones);
    double b = dot(ones, inverse_mu);
    double c = dot(data.mu(), inverse_mu);
    for (double target : {0.0, 0.05, 0.08, 0.12, 0.2})
    {
        SCOPED_TRACE(target);
        Spec spec;
        spec.target_return = target;
        Result result = solve_checked(spec, data);
        double expected = (a * target * target - 2.0 * b * target + c) / (a * c - b * b);
        EXPECT_NEAR(variance(data, result), expected, 1e-10 * expected);
    }
}

TEST(MinVariance, LongOnlyBoundsKkt)
{
    for (std::uint64_t seed = 1; seed <= 20; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = random_data(20, seed);
        Result result = solve_checked(long_only(), data);
        for (double weight : result.weights)
        {
            EXPECT_GE(weight, -1e-14);
        }
    }
}

TEST(MinVariance, AssetBounds)
{
    constexpr std::size_t n = 6;
    Data data = random_data(n, 5);
    Result free = po::solve(Spec{}, data).value();
    auto largest = std::max_element(free.weights.begin(), free.weights.end());
    std::size_t capped = static_cast<std::size_t>(largest - free.weights.begin());
    Spec spec;
    spec.asset_lower.assign(n, -kInfinity);
    spec.asset_upper.assign(n, kInfinity);
    spec.asset_upper[capped] = 0.5 * *largest;
    Result result = solve_checked(spec, data);
    EXPECT_NEAR(result.weights[capped], spec.asset_upper[capped], 1e-14);
    EXPECT_LT(result.bound_multipliers[capped], 0.0);
}

TEST(MinVariance, InfeasibleBounds)
{
    Data data = random_data(3, 2);
    Spec spec;
    spec.upper = 0.3;
    expect_no_point(po::solve(spec, data).value(), 3);
}

TEST(MinVariance, CashWithoutTargetIsAllCash)
{
    constexpr std::size_t n = 5;
    Data data = random_data(n, 4);
    data.set_risk_free(kRiskFree);
    Result result = solve_checked(Spec{}, data);
    ASSERT_EQ(result.weights.size(), n + 1);
    for (std::size_t i = 0; i < n; ++i)
    {
        EXPECT_NEAR(result.weights[i], 0.0, 1e-15);
    }
    EXPECT_NEAR(result.weights[n], 1.0, 1e-15);
    EXPECT_TRUE(std::isnan(result.budget_multiplier));
}

TEST(MinVariance, CashFrontierIsCapitalMarketLine)
{
    constexpr std::size_t n = 7;
    Data data = random_data(n, 6);
    data.set_risk_free(kRiskFree);
    std::vector<double> excess = excess_mu(data);
    double sharpe_squared = dot(excess, covariance_solve(data, excess));
    for (double target : {0.0, 0.05, 0.1, 0.3})
    {
        SCOPED_TRACE(target);
        Spec spec;
        spec.target_return = target;
        Result result = solve_checked(spec, data);
        double expected = std::fabs(target - kRiskFree) / std::sqrt(sharpe_squared);
        EXPECT_NEAR(std::sqrt(variance(data, result)), expected, 1e-10 * expected);
    }
}

TEST(MinVariance, CashBoundsApply)
{
    constexpr std::size_t n = 6;
    Data data = random_data(n, 8);
    data.set_risk_free(kRiskFree);
    double best = *std::max_element(data.mu().begin(), data.mu().end());
    Spec spec;
    spec.target_return = kRiskFree + 0.9 * (best - kRiskFree);
    Result levered = po::solve(spec, data).value();
    ASSERT_LT(levered.weights[n], 0.0);
    spec.lower = 0.0;
    spec.upper = 1.0;
    Result result = solve_checked(spec, data);
    EXPECT_GE(result.weights[n], -1e-14);
    EXPECT_LE(result.weights[n], 1.0 + 1e-14);
}

TEST(MinVariance, AssetBoundsSizeIncludesCash)
{
    constexpr std::size_t n = 4;
    Data data = random_data(n, 9);
    data.set_risk_free(kRiskFree);
    Spec spec;
    spec.asset_lower.assign(n, 0.0);
    std::expected<Result, po::Failure> rejected = po::solve(spec, data);
    ASSERT_FALSE(rejected.has_value());
    EXPECT_EQ(rejected.error().kind, po::FailureKind::InvalidArgument);
    EXPECT_EQ(rejected.error().message, "asset_lower has 4 entries; expected 5");
    spec.asset_lower.assign(n + 1, 0.0);
    solve_checked(spec, data);
}

TEST(MinVariance, IllConditionedCovariance)
{
    Data data = random_data(10, 7, 1e10);
    solve_checked(long_only(), data, 1e-8);
}

TEST(MinVariance, TargetReturnEqualToCommonMean)
{
    constexpr std::size_t n = 5;
    Data data = Data::create(n, std::vector<double>(n, 0.1), po::testing::random_covariance(n, 10)).value();
    Spec spec;
    spec.target_return = 0.1;
    Result result = solve_checked(spec, data);
    Result free = po::solve(Spec{}, data).value();
    expect_weights_near(result.weights, free.weights, 1e-12);
}

TEST(MinVariance, FixedAsset)
{
    constexpr std::size_t n = 6;
    Data data = random_data(n, 12);
    Spec spec;
    spec.asset_lower.assign(n, -kInfinity);
    spec.asset_upper.assign(n, kInfinity);
    spec.asset_lower[0] = 0.25;
    spec.asset_upper[0] = 0.25;
    Result result = solve_checked(spec, data);
    EXPECT_NEAR(result.weights[0], 0.25, 1e-14);
}

TEST(MinVariance, SingleAsset)
{
    Data data = Data::create(1, std::vector<double>{0.1}, std::vector<double>{0.04}).value();
    Result result = solve_checked(Spec{}, data);
    EXPECT_NEAR(result.weights[0], 1.0, 1e-15);
    Spec spec;
    spec.target_return = 0.2;
    expect_no_point(po::solve(spec, data).value(), 1);
}

} // namespace
