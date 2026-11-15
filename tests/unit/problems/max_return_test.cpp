#include "problems/max_return.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <vector>

#include "kkt.hpp"
#include "linalg/cholesky.hpp"
#include "problems/solve.hpp"
#include "random_covariance.hpp"
#include "user_kkt.hpp"

namespace
{

using po::Data;
using po::Objective;
using po::Result;
using po::SolveStatus;
using po::Spec;
using po::testing::covariance_times;
using po::testing::dot;
using po::testing::excess_mu;

constexpr double kRiskFree = 0.02;

Data random_data(std::size_t n, std::uint64_t seed)
{
    std::vector<double> mu = po::testing::random_vector(n, seed + 500, 0.03, 0.15);
    return Data::create(n, mu, po::testing::random_covariance(n, seed)).value();
}

std::vector<double> covariance_solve(const Data& data, std::vector<double> rhs)
{
    std::size_t n = data.asset_count();
    std::vector<double> factor = po::linalg::cholesky(data.covariance(), n).value();
    po::linalg::solve_lower(factor, n, rhs);
    po::linalg::solve_lower_transposed(factor, n, rhs);
    return rhs;
}

Spec max_return(Spec bounds, double cap)
{
    bounds.objective = Objective::MaxReturn;
    bounds.max_volatility = cap;
    return bounds;
}

Spec long_only()
{
    Spec spec;
    spec.lower = 0.0;
    spec.upper = 1.0;
    return spec;
}

double volatility(const Data& data, std::span<const double> weights)
{
    std::span<const double> risky = weights.first(data.asset_count());
    return std::sqrt(dot(risky, covariance_times(data, risky)));
}

double minimum_volatility(const Spec& bounds, const Data& data)
{
    return volatility(data, po::solve(bounds, data).value().weights);
}

double excess_return(const Data& data, const Result& result)
{
    std::span<const double> risky = std::span<const double>(result.weights).first(data.asset_count());
    return dot(risky, excess_mu(data));
}

// f = -(mu - rf 1)'w, written here independently of the builder.
Result solve_checked(const Spec& spec, const Data& data, double tolerance = 1e-10)
{
    Result result = po::solve(spec, data).value();
    std::vector<double> gradient = excess_mu(data);
    for (double& entry : gradient)
    {
        entry = -entry;
    }
    po::testing::expect_user_kkt(spec, data, result, gradient, tolerance);
    return result;
}

void expect_near_relative(std::span<const double> actual, std::span<const double> expected, double relative)
{
    ASSERT_GE(actual.size(), expected.size());
    double tolerance = relative * po::testing::largest_magnitude(expected);
    for (std::size_t i = 0; i < expected.size(); ++i)
    {
        EXPECT_NEAR(actual[i], expected[i], tolerance);
    }
}

double largest_asset_volatility(const Data& data)
{
    std::size_t n = data.asset_count();
    double largest = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        largest = std::max(largest, std::sqrt(data.covariance()[i * n + i]));
    }
    return largest;
}

TEST(MaxReturn, UnconstrainedIsUpperFrontier)
{
    constexpr std::size_t n = 8;
    Data data = random_data(n, 11);
    std::vector<double> ones(n, 1.0);
    std::vector<double> inverse_ones = covariance_solve(data, ones);
    std::vector<double> inverse_mu =
        covariance_solve(data, std::vector<double>(data.mu().begin(), data.mu().end()));
    double a = dot(ones, inverse_ones);
    double b = dot(ones, inverse_mu);
    double c = dot(data.mu(), inverse_mu);
    double cap = 1.5 / std::sqrt(a);
    Result result = solve_checked(max_return(Spec{}, cap), data);
    double best = dot(result.weights, data.mu());
    double frontier_variance = (a * best * best - 2.0 * b * best + c) / (a * c - b * b);
    EXPECT_NEAR(frontier_variance, cap * cap, 1e-10 * cap * cap);
    EXPECT_GE(best, b / a);
    Spec target;
    target.target_return = best;
    expect_near_relative(result.weights, po::solve(target, data).value().weights, 1e-9);
    EXPECT_LT(result.volatility_multiplier, 0.0);
}

TEST(MaxReturn, CashIsCapitalMarketLine)
{
    constexpr std::size_t n = 7;
    constexpr double kCap = 0.15;
    Data data = random_data(n, 6);
    data.set_risk_free(kRiskFree);
    std::vector<double> excess = excess_mu(data);
    std::vector<double> direction = covariance_solve(data, excess);
    double scale = kCap / std::sqrt(dot(excess, direction));
    for (double& entry : direction)
    {
        entry *= scale;
    }
    Result result = solve_checked(max_return(Spec{}, kCap), data);
    expect_near_relative(std::span<const double>(result.weights).first(n), direction, 1e-10);
    EXPECT_TRUE(std::isnan(result.budget_multiplier));
}

TEST(MaxReturn, BelowMinimumIsInfeasible)
{
    Data data = random_data(6, 2);
    double minimum = minimum_volatility(long_only(), data);
    Result result = po::solve(max_return(long_only(), 0.9 * minimum), data).value();
    EXPECT_EQ(result.status, SolveStatus::Infeasible);
    ASSERT_EQ(result.weights.size(), 6u);
    for (double weight : result.weights)
    {
        EXPECT_TRUE(std::isnan(weight));
    }
    EXPECT_TRUE(std::isnan(result.volatility_multiplier));
    EXPECT_GT(result.iterations, 0);
}

// Three assets whose minimum variance with w >= 0 is all in asset 0. The lower-bound multipliers at that
// vertex are z_1 = covariance_01 - 0.01 - 0.03 t and z_2 = 0.015 - 0.15 t, so w(t) leaves it at
// t_b = min((covariance_01 - 0.01) / 0.03, 0.1).
Data vertex_data(double covariance_01)
{
    const std::vector<double> mu = {0.05, 0.08, 0.2};
    const std::vector<double> covariance = {
        0.01,          covariance_01, 0.025,
        covariance_01, 0.04,          0.03,
        0.025,         0.03,          0.25,
    };
    return Data::create(3, mu, covariance).value();
}

Spec nonnegative()
{
    Spec spec;
    spec.lower = 0.0;
    return spec;
}

TEST(MaxReturn, AtMinimumIsMinVariance)
{
    Data data = random_data(6, 2);
    Result minimum = po::solve(long_only(), data).value();
    Result result = po::solve(max_return(long_only(), volatility(data, minimum.weights)), data).value();
    EXPECT_EQ(result.status, SolveStatus::Optimal);
    expect_near_relative(result.weights, minimum.weights, 1e-12);
    EXPECT_EQ(result.volatility_multiplier, -std::numeric_limits<double>::infinity());
    EXPECT_TRUE(std::isnan(result.budget_multiplier));
    for (double beta : result.bound_multipliers)
    {
        EXPECT_TRUE(std::isnan(beta));
    }

    Data equal = Data::create(6, std::vector<double>(6, 0.1), po::testing::random_covariance(6, 2)).value();
    Result free = po::solve(Spec{}, equal).value();
    Result at_minimum = solve_checked(max_return(Spec{}, volatility(equal, free.weights)), equal);
    expect_near_relative(at_minimum.weights, free.weights, 1e-12);
    EXPECT_EQ(at_minimum.volatility_multiplier, 0.0);
}

TEST(MaxReturn, AtMinimumOnVertexIsRightDerivative)
{
    struct Case
    {
        double covariance_01;
        double exit;
    };
    for (Case vertex : {Case{0.012, 0.002 / 0.03}, Case{0.0105, 0.0005 / 0.03}})
    {
        SCOPED_TRACE(vertex.covariance_01);
        Data data = vertex_data(vertex.covariance_01);
        Result minimum = po::solve(nonnegative(), data).value();
        ASSERT_NEAR(minimum.weights[0], 1.0, 1e-14);
        double cap = volatility(data, minimum.weights);
        Result result = solve_checked(max_return(nonnegative(), cap), data);
        expect_near_relative(result.weights, minimum.weights, 1e-14);
        EXPECT_NEAR(result.volatility_multiplier, -cap / vertex.exit, 1e-12 * cap / vertex.exit);
        double step = 1e-6 * cap;
        double optimum = -excess_return(data, result);
        double far = -excess_return(data, solve_checked(max_return(nonnegative(), cap + step), data));
        double near = -excess_return(data, solve_checked(max_return(nonnegative(), cap + 0.5 * step), data));
        double derivative = (4.0 * near - far - 3.0 * optimum) / step;
        EXPECT_NEAR(result.volatility_multiplier, derivative, 1e-5 * std::fabs(derivative));
    }
}

TEST(MaxReturn, InactiveVolatilityLongOnly)
{
    constexpr std::size_t n = 8;
    Data data = random_data(n, 4);
    Spec bounds;
    bounds.lower = 0.0;
    Result result = solve_checked(max_return(bounds, 10.0 * largest_asset_volatility(data)), data);
    std::size_t best = static_cast<std::size_t>(std::max_element(data.mu().begin(), data.mu().end()) -
                                                data.mu().begin());
    for (std::size_t i = 0; i < n; ++i)
    {
        double expected = 0.0;
        if (i == best)
        {
            expected = 1.0;
        }
        EXPECT_NEAR(result.weights[i], expected, 1e-14);
    }
    EXPECT_EQ(result.volatility_multiplier, 0.0);
}

// The lower-bound multipliers at the minimum-variance vertex shrink as t grows: w(t) stays at that vertex
// over the first doubling and then leaves it.
TEST(MaxReturn, MinimumAtVertexLeavesIt)
{
    Data data = vertex_data(0.012);
    Result result = solve_checked(max_return(nonnegative(), 0.2), data);
    EXPECT_LT(result.volatility_multiplier, 0.0);
    EXPECT_NEAR(volatility(data, result.weights), 0.2, 1e-12 * 0.2);
}

TEST(MaxReturn, MultiplierIsDerivativeOfOptimum)
{
    constexpr std::size_t n = 10;
    Data data = random_data(n, 12);
    std::size_t best = static_cast<std::size_t>(std::max_element(data.mu().begin(), data.mu().end()) -
                                                data.mu().begin());
    double minimum = minimum_volatility(long_only(), data);
    double cap = 0.5 * (minimum + std::sqrt(data.covariance()[best * n + best]));
    double step = 1e-6 * cap;
    Result result = solve_checked(max_return(long_only(), cap), data);
    double above = -excess_return(data, solve_checked(max_return(long_only(), cap + step), data));
    double below = -excess_return(data, solve_checked(max_return(long_only(), cap - step), data));
    double derivative = (above - below) / (2.0 * step);
    EXPECT_LT(result.volatility_multiplier, 0.0);
    EXPECT_NEAR(result.volatility_multiplier, derivative, 1e-5 * std::fabs(derivative));
}

TEST(MaxReturn, DailyScale)
{
    constexpr std::size_t n = 9;
    std::vector<double> covariance = po::testing::random_covariance(n, 13);
    for (double& entry : covariance)
    {
        entry *= 1e-4;
    }
    std::vector<double> mu = po::testing::random_vector(n, 513, 0.03, 0.15);
    for (double& entry : mu)
    {
        entry *= 1e-2;
    }
    Data data = Data::create(n, mu, covariance).value();
    double cap = 1.2 * minimum_volatility(Spec{}, data);
    Result result = solve_checked(max_return(Spec{}, cap), data);
    EXPECT_NEAR(volatility(data, result.weights), cap, 1e-12 * cap);
    EXPECT_LT(result.volatility_multiplier, 0.0);
}

TEST(MaxReturn, EqualMeansIsMinVariance)
{
    constexpr std::size_t n = 7;
    Data data = Data::create(n, std::vector<double>(n, 0.1), po::testing::random_covariance(n, 14)).value();
    Result minimum = po::solve(Spec{}, data).value();
    Result result = solve_checked(max_return(Spec{}, 2.0 * volatility(data, minimum.weights)), data);
    expect_near_relative(result.weights, minimum.weights, 1e-12);
    EXPECT_EQ(result.volatility_multiplier, 0.0);
}

// Weights at zero in the minimum variance have bound multipliers whose slope in t is zero.
TEST(MaxReturn, EqualMeansLongOnly)
{
    int with_zero_weights = 0;
    for (std::uint64_t seed = 1; seed <= 40; ++seed)
    {
        SCOPED_TRACE(seed);
        std::size_t n = 4 + static_cast<std::size_t>(seed % 15);
        std::vector<double> covariance = po::testing::random_covariance(n, seed);
        Data data = Data::create(n, std::vector<double>(n, 0.1), covariance).value();
        Result minimum = po::solve(long_only(), data).value();
        if (std::ranges::none_of(minimum.weights, [](double weight) { return std::fabs(weight) <= 1e-14; }))
        {
            continue;
        }
        ++with_zero_weights;
        Result result = solve_checked(max_return(long_only(), 2.0 * volatility(data, minimum.weights)), data);
        expect_near_relative(result.weights, minimum.weights, 1e-12);
        EXPECT_EQ(result.volatility_multiplier, 0.0);
    }
    EXPECT_GE(with_zero_weights, 30);
}

TEST(MaxReturn, RandomLongOnly)
{
    for (std::uint64_t seed = 1; seed <= 30; ++seed)
    {
        SCOPED_TRACE(seed);
        std::size_t n = 3 + static_cast<std::size_t>(seed % 23);
        Data data = random_data(n, seed);
        double minimum = minimum_volatility(long_only(), data);
        double fraction = po::testing::random_vector(1, seed + 900, 0.0, 1.0)[0];
        double cap = minimum + fraction * (largest_asset_volatility(data) - minimum);
        Result result = solve_checked(max_return(long_only(), cap), data);
        double achieved = volatility(data, result.weights);
        if (result.volatility_multiplier == 0.0)
        {
            EXPECT_LT(achieved, cap);
            continue;
        }
        EXPECT_LT(result.volatility_multiplier, 0.0);
        EXPECT_NEAR(achieved, cap, 1e-12 * cap);
    }
}

} // namespace
