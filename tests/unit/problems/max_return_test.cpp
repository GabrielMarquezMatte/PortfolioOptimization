#include "problems/max_return.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <span>
#include <vector>

#include "kkt.hpp"
#include "linalg/cholesky.hpp"
#include "problems/parametric_path.hpp"
#include "qp/goldfarb_idnani.hpp"
#include "random_covariance.hpp"
#include "solve/solve.hpp"
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

Data random_data(std::size_t n, std::uint64_t seed, double condition = 100.0)
{
    std::vector<double> mu = po::testing::random_vector(n, seed + 500, 0.03, 0.15);
    return Data::create(n, mu, po::testing::random_covariance(n, seed, condition)).value();
}

// vol(w) = s of a solved point, relative.
constexpr double kVolatilityTolerance = 1e-12;

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
    po::testing::expect_user_kkt(spec, data, result, po::testing::sum_of({gradient}), tolerance);
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

// The QP solutions of Q(t) have a stationarity residual that grows with cond(Sigma), to 2.2e-10 at 1e7 on
// these problems (measured on qp::solve alone): above 1e6 the KKT conditions are held to 1e-8, the tolerance
// of an ill-conditioned covariance.
double sweep_tolerance(double condition)
{
    if (condition > 1e6)
    {
        return 1e-8;
    }
    return 1e-10;
}

// The QP rounding of vol(w(t*)) grows with cond(Sigma) (about 1e-15 cond s): the crossing is accepted on the
// segment it lies on, and the reported point is the exact root of that segment.
TEST(MaxReturn, ConditioningSweep)
{
    for (double condition : {1e2, 1e4, 1e5, 1e6, 1e7})
    {
        SCOPED_TRACE(condition);
        for (std::uint64_t seed = 1; seed <= 8; ++seed)
        {
            SCOPED_TRACE(seed);
            std::size_t n = 5 + static_cast<std::size_t>((7 * seed) % 36);
            Data data = random_data(n, seed, condition);
            double minimum = minimum_volatility(long_only(), data);
            for (double fraction : {0.1, 0.3, 0.6, 0.9})
            {
                SCOPED_TRACE(fraction);
                double cap = minimum + fraction * (largest_asset_volatility(data) - minimum);
                Result result = solve_checked(max_return(long_only(), cap), data, sweep_tolerance(condition));
                if (result.volatility_multiplier == 0.0)
                {
                    EXPECT_LT(volatility(data, result.weights), cap);
                    continue;
                }
                EXPECT_NEAR(volatility(data, result.weights), cap, kVolatilityTolerance * cap);
            }
        }
    }
}

// Cash fixed at 0.1 by its bounds, the risky weights long-only, cond(Sigma) = 1e8 and means 0.1 (1 + 1e-6 u_i): Q(t)
// far along the path starts at t Sigma^-1 mu~, far from the final w, and the rounding of the steps from there leaves
// the active bound of the cash off by more than the tolerance at w. Its other bound, violated by that alone, made a
// Q(t) infeasible, and the path broke down at t > 0.
TEST(MaxReturn, FixedCashIllConditioned)
{
    for (double fraction : {0.3, 0.9})
    {
        for (std::uint64_t seed = 1; seed <= 40; ++seed)
        {
            SCOPED_TRACE(testing::Message() << "fraction " << fraction << " seed " << seed);
            std::size_t n = 3 + seed % 6;
            std::vector<double> mu = po::testing::random_vector(n, seed + 500, -1.0, 1.0);
            for (double& value : mu)
            {
                value = 0.1 * (1.0 + 1e-6 * value);
            }
            Data data = Data::create(n, mu, po::testing::random_covariance(n, seed, 1e8)).value();
            data.set_risk_free(kRiskFree);
            Spec bounds;
            bounds.asset_lower.assign(n + 1, 0.0);
            bounds.asset_upper.assign(n + 1, 1.0);
            bounds.asset_lower[n] = 0.1;
            bounds.asset_upper[n] = 0.1;
            double cap = (1.0 + fraction) * minimum_volatility(bounds, data);
            Result result = solve_checked(max_return(bounds, cap), data, sweep_tolerance(1e8));
            ASSERT_EQ(result.status, SolveStatus::Optimal);
            EXPECT_NEAR(result.weights[n], 0.1, 1e-12);
        }
    }
}

// Without cash the path solves Q(t) for mu less its midrange (1'w = 1 makes it the same w(t)) and maps the
// budget multiplier back: adding k to every mean leaves the weights and the volatility multiplier and lowers
// the budget multiplier by k, at a crossing and with an inactive cap.
TEST(MaxReturn, BudgetMultiplierFollowsCommonMean)
{
    constexpr std::size_t n = 8;
    constexpr double kCommon = 1.0;
    Data data = random_data(n, 21);
    std::vector<double> lifted(data.mu().begin(), data.mu().end());
    for (double& value : lifted)
    {
        value += kCommon;
    }
    Data shifted = Data::create(n, lifted, std::vector<double>(data.covariance().begin(), data.covariance().end()))
                       .value();
    double minimum = minimum_volatility(long_only(), data);
    for (double cap : {0.5 * (minimum + largest_asset_volatility(data)), 10.0 * largest_asset_volatility(data)})
    {
        SCOPED_TRACE(cap);
        Result base = solve_checked(max_return(long_only(), cap), data);
        Result result = solve_checked(max_return(long_only(), cap), shifted);
        expect_near_relative(result.weights, base.weights, 1e-12);
        EXPECT_NEAR(result.volatility_multiplier, base.volatility_multiplier,
                    1e-12 * std::fabs(base.volatility_multiplier));
        EXPECT_NEAR(result.budget_multiplier, base.budget_multiplier - kCommon, 1e-12 * kCommon);
    }
}

// Every excess return is negative: all cash is the maximum return at any cap, and the cap is inactive.
TEST(MaxReturn, RiskFreeAboveAllMeans)
{
    for (std::uint64_t seed = 1; seed <= 6; ++seed)
    {
        SCOPED_TRACE(seed);
        std::size_t n = 2 + static_cast<std::size_t>((5 * seed) % 28);
        Data data = random_data(n, seed);
        data.set_risk_free(0.2);
        for (double cap : {1e-4, 0.01, 0.1, 0.5})
        {
            SCOPED_TRACE(cap);
            Result result = solve_checked(max_return(long_only(), cap), data);
            for (std::size_t i = 0; i < n; ++i)
            {
                EXPECT_NEAR(result.weights[i], 0.0, 1e-14);
            }
            EXPECT_NEAR(result.weights[n], 1.0, 1e-14);
            EXPECT_EQ(result.volatility_multiplier, 0.0);
        }
    }
}

Spec capped_long_only(double upper)
{
    Spec spec;
    spec.lower = 0.0;
    spec.upper = upper;
    return spec;
}

// With equal means the budget takes up all of t (mu - rf 1): w(t) is the minimum variance for every t, with
// interior weights, weights at zero and weights at the upper bound, and the cap is inactive at every s, even
// at s = sigma_min.
TEST(MaxReturn, EqualMeansTightUpperBounds)
{
    constexpr std::size_t n = 50;
    constexpr double kUpper = 0.05;
    for (std::uint64_t seed = 1; seed <= 4; ++seed)
    {
        SCOPED_TRACE(seed);
        std::vector<double> covariance = po::testing::random_covariance(n, seed, 1e4);
        Data data = Data::create(n, std::vector<double>(n, 0.1), covariance).value();
        Result minimum = po::solve(capped_long_only(kUpper), data).value();
        ASSERT_EQ(minimum.status, SolveStatus::Optimal);
        auto at_upper = [](double weight) { return weight > kUpper - 1e-14; };
        ASSERT_TRUE(std::ranges::any_of(minimum.weights, at_upper));
        double sigma_min = volatility(data, minimum.weights);
        for (double factor : {1.0, 1.01, 2.0, 1e4})
        {
            SCOPED_TRACE(factor);
            Result result = solve_checked(max_return(capped_long_only(kUpper), factor * sigma_min), data);
            expect_near_relative(result.weights, minimum.weights, 1e-10);
            EXPECT_EQ(result.volatility_multiplier, 0.0);
        }
    }
}

// Equal means with cash: w(t) grows along the capital market line until the cash weight reaches 0 and then
// stays at the long-only minimum variance of the risky assets, which is not a vertex.
TEST(MaxReturn, EqualMeansWithCash)
{
    for (std::size_t n : {std::size_t{6}, std::size_t{20}, std::size_t{40}})
    {
        SCOPED_TRACE(n);
        std::vector<double> covariance = po::testing::random_covariance(n, n, 1e4);
        Data risky = Data::create(n, std::vector<double>(n, 0.1), covariance).value();
        Result minimum = po::solve(long_only(), risky).value();
        double sigma_min = volatility(risky, minimum.weights);
        Data data = risky;
        data.set_risk_free(kRiskFree);
        for (double factor : {1.0, 2.0, 1e4})
        {
            SCOPED_TRACE(factor);
            Result result = solve_checked(max_return(long_only(), factor * sigma_min), data);
            expect_near_relative(result.weights, minimum.weights, 1e-10);
            EXPECT_NEAR(result.weights[n], 0.0, 1e-12);
        }
        Result below = solve_checked(max_return(long_only(), 0.5 * sigma_min), data);
        EXPECT_LT(below.volatility_multiplier, 0.0);
        EXPECT_NEAR(volatility(data, below.weights), 0.5 * sigma_min, kVolatilityTolerance * sigma_min);
    }
}

Spec leveraged(double leverage)
{
    Spec spec;
    spec.leverage = leverage;
    return spec;
}

std::size_t best_asset(const Data& data)
{
    return static_cast<std::size_t>(std::ranges::max_element(data.mu()) - data.mu().begin());
}

std::size_t worst_asset(const Data& data)
{
    return static_cast<std::size_t>(std::ranges::min_element(data.mu()) - data.mu().begin());
}

// max mu'w over 1'w = 1 and sum_i |w_i| <= L: (L + 1) / 2 in the best asset, (L - 1) / 2 short in the worst.
std::vector<double> leverage_vertex(const Data& data, double leverage)
{
    std::vector<double> weights(data.asset_count(), 0.0);
    weights[best_asset(data)] = 0.5 * (leverage + 1.0);
    weights[worst_asset(data)] = -0.5 * (leverage - 1.0);
    return weights;
}

// A cap above the volatility of the vertex of the leverage: w(t) reaches the vertex and stays there, the cap is
// inactive, and -mu = budget 1 + lambda_lev g with g = 1 at the best asset and -1 at the worst gives
// lambda_lev = (mu_worst - mu_best) / 2 and budget = -(mu_best + mu_worst) / 2.
TEST(MaxReturn, LeverageActiveAlongPath)
{
    constexpr std::size_t n = 8;
    constexpr double kLeverage = 1.5;
    for (std::uint64_t seed = 31; seed <= 34; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = random_data(n, seed);
        std::vector<double> vertex = leverage_vertex(data, kLeverage);
        ASSERT_LT(volatility(data, vertex), 2.0 * largest_asset_volatility(data));
        Result result = solve_checked(max_return(leveraged(kLeverage), 2.0 * largest_asset_volatility(data)), data);
        expect_near_relative(result.weights, vertex, 1e-12);
        EXPECT_EQ(result.volatility_multiplier, 0.0);
        double best = data.mu()[best_asset(data)];
        double worst = data.mu()[worst_asset(data)];
        EXPECT_LT(result.leverage_multiplier, 0.0);
        EXPECT_NEAR(result.leverage_multiplier, 0.5 * (worst - best), 1e-12 * best);
        EXPECT_NEAR(result.budget_multiplier, -0.5 * (best + worst), 1e-12 * best);
    }
}

double max_return_value(const Data& data, double cap, double leverage)
{
    return -excess_return(data, solve_checked(max_return(leveraged(leverage), cap), data));
}

// A cap between the minimum volatility under the leverage and the volatility of the vertex of the leverage, where
// both bind: each multiplier is the derivative of the optimum of f = -mu'w in its right-hand side.
TEST(MaxReturn, LeverageAndVolatilityBothActive)
{
    constexpr std::size_t n = 8;
    constexpr double kLeverage = 1.5;
    constexpr double kLeverageStep = 1e-4 * kLeverage;
    for (std::uint64_t seed = 31; seed <= 34; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = random_data(n, seed);
        double minimum = minimum_volatility(leveraged(kLeverage), data);
        double cap = minimum + 0.8 * (volatility(data, leverage_vertex(data, kLeverage)) - minimum);
        Result result = solve_checked(max_return(leveraged(kLeverage), cap), data);
        EXPECT_LT(result.volatility_multiplier, 0.0);
        EXPECT_LT(result.leverage_multiplier, 0.0);
        EXPECT_NEAR(volatility(data, result.weights), cap, kVolatilityTolerance * cap);
        double step = 1e-6 * cap;
        double cap_derivative = (max_return_value(data, cap + step, kLeverage) -
                                 max_return_value(data, cap - step, kLeverage)) /
                                (2.0 * step);
        EXPECT_NEAR(result.volatility_multiplier, cap_derivative, 1e-5 * std::fabs(cap_derivative));
        double leverage_derivative = (max_return_value(data, cap, kLeverage + kLeverageStep) -
                                      max_return_value(data, cap, kLeverage - kLeverageStep)) /
                                     (2.0 * kLeverageStep);
        EXPECT_NEAR(result.leverage_multiplier, leverage_derivative, 1e-6 * std::fabs(leverage_derivative));
    }
}

// Means 0.1 (1 + 1e-6 u_i), u_i uniform in [-1, 1], with cash at kRiskFree: once a face that holds the whole
// portfolio binds (the leverage, or the cash bound of long-only), only the spread of 1e-7 moves w(t), the cap is
// reached at t near 1e8 and the path solves Q(t) shifted along that face (docs/guide/validation.md).
Data nearly_equal_with_cash(std::size_t n, std::uint64_t seed)
{
    std::vector<double> mu = po::testing::random_vector(n, seed + 500, -1.0, 1.0);
    for (double& value : mu)
    {
        value = 0.1 * (1.0 + 1e-6 * value);
    }
    Data data = Data::create(n, mu, po::testing::random_covariance(n, seed, 1e4)).value();
    data.set_risk_free(kRiskFree);
    return data;
}

// Under the leverage with cash, the inactive cap leaves w = L e_best, where -mu~ = lambda_lev g with g_best = 1:
// lambda_lev = -mu~_best, the face multiplier mapped back from the shifted Q(t). At a cap below the volatility of
// that vertex both bind, and expect_user_kkt checks lambda_lev with lambda_vol.
TEST(MaxReturn, FaceShiftLeverageMultiplier)
{
    constexpr std::size_t n = 6;
    constexpr double kLeverage = 1.3;
    for (std::uint64_t seed = 41; seed <= 44; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = nearly_equal_with_cash(n, seed);
        std::size_t best = best_asset(data);
        double vertex = kLeverage * std::sqrt(data.covariance()[best * n + best]);
        Result limit = solve_checked(max_return(leveraged(kLeverage), 2.0 * vertex), data);
        EXPECT_EQ(limit.volatility_multiplier, 0.0);
        EXPECT_NEAR(limit.weights[best], kLeverage, 1e-12);
        EXPECT_NEAR(limit.leverage_multiplier, -(data.mu()[best] - kRiskFree), 1e-12 * data.mu()[best]);
        Result crossing = solve_checked(max_return(leveraged(kLeverage), 0.9 * vertex), data);
        EXPECT_LT(crossing.volatility_multiplier, 0.0);
        EXPECT_LT(crossing.leverage_multiplier, 0.0);
        EXPECT_NEAR(volatility(data, crossing.weights), 0.9 * vertex, kVolatilityTolerance * vertex);
    }
}

// Every weight nonnegative with cash: the cash bound w_cash >= 0 is the face. At the inactive cap w = e_best and cash = 0, where
// -mu~_i = beta_i - beta_0 gives beta_0 = mu~_best and beta_i = mu~_best - mu~_i: the cash multiplier is the face
// multiplier mapped back. Just below the volatility of the best asset (the end of the path, reached fully invested)
// both the cap and the cash bound bind in all 8 seeds, where the cash weight is 0 to its rounding.
TEST(MaxReturn, FaceShiftCashMultiplier)
{
    constexpr std::size_t n = 6;
    constexpr double kCashRounding = 1e-12;
    int on_face = 0;
    for (std::uint64_t seed = 41; seed <= 48; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = nearly_equal_with_cash(n, seed);
        std::size_t best = best_asset(data);
        double vertex = std::sqrt(data.covariance()[best * n + best]);
        Result limit = solve_checked(max_return(nonnegative(), 2.0 * vertex), data);
        EXPECT_EQ(limit.volatility_multiplier, 0.0);
        EXPECT_NEAR(limit.weights[best], 1.0, 1e-12);
        double best_excess = data.mu()[best] - kRiskFree;
        EXPECT_NEAR(limit.bound_multipliers[n], best_excess, 1e-12 * best_excess);
        for (std::size_t i = 0; i < n; ++i)
        {
            EXPECT_NEAR(limit.bound_multipliers[i], data.mu()[best] - data.mu()[i], 1e-12 * best_excess);
        }
        double cap = 0.999 * vertex;
        Result crossing = solve_checked(max_return(nonnegative(), cap), data);
        EXPECT_LT(crossing.volatility_multiplier, 0.0);
        EXPECT_NEAR(volatility(data, crossing.weights), cap, kVolatilityTolerance * cap);
        if (crossing.weights[n] > kCashRounding)
        {
            continue;
        }
        ++on_face;
        EXPECT_GT(crossing.bound_multipliers[n], 0.0);
    }
    EXPECT_EQ(on_face, 8);
}

bool same_bits(const std::vector<double>& a, const std::vector<double>& b)
{
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

void expect_same_point(const po::problems::PathPoint& cached, const po::problems::PathPoint& fresh)
{
    const po::qp::Solution& a = cached.solution;
    const po::qp::Solution& b = fresh.solution;
    EXPECT_EQ(a.status, b.status);
    EXPECT_EQ(a.iterations, b.iterations);
    EXPECT_TRUE(same_bits(a.x, b.x));
    EXPECT_TRUE(same_bits(a.equality_multipliers, b.equality_multipliers));
    EXPECT_TRUE(same_bits(a.inequality_multipliers, b.inequality_multipliers));
    EXPECT_EQ(a.active_inequalities, b.active_inequalities);
    EXPECT_TRUE(same_bits(a.gross_exposure_signs, b.gross_exposure_signs));
    EXPECT_TRUE(same_bits(a.gross_exposure_multipliers, b.gross_exposure_multipliers));
    EXPECT_EQ(std::memcmp(&cached.volatility, &fresh.volatility, sizeof(double)), 0);
    EXPECT_TRUE(same_bits(cached.shift.excess, fresh.shift.excess));
}

// Every point of the path is solved with the factor of Data, held unchanged by the path's problem, so the
// basis of the first point is reused at every later one: the points equal those of a path whose cache is
// emptied before each of them.
TEST(MaxReturn, PathBasisCacheIsBitIdenticalToAFreshBasis)
{
    constexpr std::size_t n = 12;
    Data data = random_data(n, 41);
    Data with_cash = random_data(n, 41);
    with_cash.set_risk_free(kRiskFree);
    Spec box = long_only();
    box.upper = 0.3;
    struct Case
    {
        Spec bounds;
        const Data* data;
    };
    for (const Case& shape : {Case{Spec{}, &data}, Case{long_only(), &data}, Case{box, &data},
                              Case{leveraged(1.5), &data}, Case{long_only(), &with_cash}})
    {
        Spec spec = max_return(shape.bounds, 1.0);
        po::problems::Path cached = po::problems::make_path(spec, *shape.data);
        po::problems::Path fresh = po::problems::make_path(spec, *shape.data);
        std::vector<double> factor = cached.translation.problem.hessian_factor;
        for (double t = 0.0; t < 40.0; t = std::max(2.0 * t, 1e-4))
        {
            SCOPED_TRACE(testing::Message() << "leverage " << spec.leverage.has_value() << ", t = " << t);
            po::problems::PathPoint point = po::problems::solve_at(cached, t);
            fresh.basis = po::qp::BasisCache{};
            expect_same_point(point, po::problems::solve_at(fresh, t));
            EXPECT_TRUE(same_bits(cached.translation.problem.hessian_factor, factor));
            EXPECT_TRUE(same_bits(cached.basis.factor, factor));
        }
        EXPECT_EQ(cached.iterations, fresh.iterations);
    }
}

} // namespace
