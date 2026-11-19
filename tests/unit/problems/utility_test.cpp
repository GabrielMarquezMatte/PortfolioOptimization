#include "problems/mean_variance.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "kkt.hpp"
#include "linalg/cholesky.hpp"
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

// f = (quadratic / 2) w'Sigma w - linear (mu - rf 1)'w, written here independently of the builder.
struct Scales
{
    double quadratic;
    double linear;
};

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

// grad f = quadratic Sigma w - linear (mu - rf 1), as its two terms.
po::testing::ObjectiveGradient objective_gradient(const Data& data, const Result& result, Scales scales)
{
    std::vector<double> curvature = covariance_times(data, result.weights);
    std::vector<double> excess = excess_mu(data);
    for (std::size_t i = 0; i < curvature.size(); ++i)
    {
        curvature[i] *= scales.quadratic;
        excess[i] *= -scales.linear;
    }
    return po::testing::sum_of({curvature, excess});
}

Result solve_checked(const Spec& spec, const Data& data, Scales scales, double tolerance = 1e-10)
{
    po::problems::MeanVariance builder_scales = po::problems::mean_variance_of(spec, data);
    po::problems::Translation translation = po::problems::mean_variance_problem(spec, data, builder_scales);
    po::testing::expect_kkt(translation.problem, po::qp::solve(translation.problem), tolerance);
    Result result = po::solve(spec, data).value();
    po::testing::expect_user_kkt(spec, data, result, objective_gradient(data, result, scales), tolerance);
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

// The same value, or both NaN.
void expect_same(double actual, double expected)
{
    EXPECT_EQ(std::isnan(actual), std::isnan(expected));
    if (!std::isnan(expected))
    {
        EXPECT_EQ(actual, expected);
    }
}

void expect_identical(const Result& actual, const Result& expected)
{
    EXPECT_EQ(actual.status, expected.status);
    EXPECT_EQ(actual.weights, expected.weights);
    EXPECT_EQ(actual.bound_multipliers, expected.bound_multipliers);
    expect_same(actual.budget_multiplier, expected.budget_multiplier);
    expect_same(actual.target_return_multiplier, expected.target_return_multiplier);
    expect_same(actual.volatility_multiplier, expected.volatility_multiplier);
    expect_same(actual.leverage_multiplier, expected.leverage_multiplier);
    EXPECT_EQ(actual.residual, expected.residual);
    EXPECT_EQ(actual.iterations, expected.iterations);
}

Spec utility(double risk_aversion)
{
    Spec spec;
    spec.objective = Objective::Utility;
    spec.risk_aversion = risk_aversion;
    return spec;
}

Spec importance(double lambda)
{
    Spec spec;
    spec.objective = Objective::Importance;
    spec.importance = lambda;
    return spec;
}

Spec kelly()
{
    Spec spec;
    spec.objective = Objective::Kelly;
    return spec;
}

Spec long_only(Spec spec)
{
    spec.lower = 0.0;
    spec.upper = 1.0;
    return spec;
}

std::vector<double> scaled(std::vector<double> values, double factor)
{
    for (double& value : values)
    {
        value *= factor;
    }
    return values;
}

TEST(Utility, UnconstrainedMatchesClosedForm)
{
    constexpr std::size_t n = 6;
    constexpr double kDelta = 3.0;
    Data data = random_data(n, 3);
    std::vector<double> ones(n, 1.0);
    std::vector<double> inverse_ones = covariance_solve(data, ones);
    std::vector<double> inverse_mu = covariance_solve(data, std::vector<double>(data.mu().begin(), data.mu().end()));
    double budget_share = (1.0 - dot(ones, inverse_mu) / kDelta) / dot(ones, inverse_ones);
    std::vector<double> expected(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        expected[i] = inverse_mu[i] / kDelta + budget_share * inverse_ones[i];
    }
    Result result = solve_checked(utility(kDelta), data, Scales{kDelta, 1.0});
    expect_near_relative(result.weights, expected, 1e-12);
}

TEST(Utility, CashMatchesClosedForm)
{
    constexpr std::size_t n = 7;
    constexpr double kDelta = 4.0;
    Data data = random_data(n, 6);
    data.set_risk_free(kRiskFree);
    std::vector<double> expected = scaled(covariance_solve(data, excess_mu(data)), 1.0 / kDelta);
    double cash = 1.0;
    for (double weight : expected)
    {
        cash -= weight;
    }
    expected.push_back(cash);
    Result result = solve_checked(utility(kDelta), data, Scales{kDelta, 1.0});
    expect_near_relative(result.weights, expected, 1e-12);
    EXPECT_TRUE(std::isnan(result.budget_multiplier));
}

TEST(Utility, DefaultRiskAversionIsFive)
{
    Data data = random_data(8, 4);
    Spec spec;
    spec.objective = Objective::Utility;
    Result result = solve_checked(spec, data, Scales{5.0, 1.0});
    expect_identical(result, po::solve(utility(5.0), data).value());
}

TEST(Utility, LongOnlyKkt)
{
    for (std::uint64_t seed = 1; seed <= 20; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = random_data(20, seed);
        Result result = solve_checked(long_only(utility(5.0)), data, Scales{5.0, 1.0});
        for (double weight : result.weights)
        {
            EXPECT_GE(weight, -1e-14);
        }
    }
}

TEST(Utility, LargeRiskAversionApproachesMinVariance)
{
    constexpr double kDelta = 1e8;
    Data data = random_data(8, 9);
    Result result = solve_checked(utility(kDelta), data, Scales{kDelta, 1.0});
    Result minimum = po::solve(Spec{}, data).value();
    expect_near_relative(result.weights, minimum.weights, 1e-6);
}

Spec with_asset_bounds(Spec spec, std::size_t count, double fixed)
{
    spec.asset_lower.assign(count, 0.0);
    spec.asset_upper.assign(count, 1.0);
    spec.asset_lower[0] = fixed;
    spec.asset_upper[0] = fixed;
    return spec;
}

// A small risk aversion puts the unconstrained start near 1e4 times the final weights, and the rounding of the steps
// from it leaves the active bound of the fixed weight off by more than the tolerance there: its other bound, violated
// by that alone, made the problem infeasible.
TEST(Utility, FixedWeightAtSmallRiskAversion)
{
    constexpr std::size_t n = 6;
    constexpr double kSmallDelta = 1e-5;
    for (std::uint64_t seed = 1; seed <= 40; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = random_data(n, seed);
        Result result = solve_checked(with_asset_bounds(utility(kSmallDelta), n, 0.2), data,
                                      Scales{kSmallDelta, 1.0});
        ASSERT_EQ(result.status, SolveStatus::Optimal);
        EXPECT_NEAR(result.weights[0], 0.2, 1e-12);
    }
}

// Cash last, fixed by its bounds, with the risky weights long-only.
TEST(Utility, FixedCashAtSmallRiskAversion)
{
    constexpr std::size_t n = 5;
    constexpr double kSmallDelta = 1e-4;
    for (std::uint64_t seed = 1; seed <= 40; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = random_data(n, seed);
        data.set_risk_free(kRiskFree);
        Spec spec = utility(kSmallDelta);
        spec.asset_lower.assign(n + 1, 0.0);
        spec.asset_upper.assign(n + 1, 1.0);
        spec.asset_lower[n] = 0.15;
        spec.asset_upper[n] = 0.15;
        Result result = solve_checked(spec, data, Scales{kSmallDelta, 1.0});
        ASSERT_EQ(result.status, SolveStatus::Optimal);
        EXPECT_NEAR(result.weights[n], 0.15, 1e-12);
    }
}

// The leverage 1 without bounds and without cash: the feasible set is the simplex, and the row of the gross exposure
// family at s = 1 is the negation of the budget, violated by the budget's residual alone. A small risk aversion puts
// the unconstrained start far from the final weights.
TEST(Utility, UnitLeverageIsTheSimplex)
{
    for (double risk_aversion : {1e-4, 1e-5})
    {
        for (std::uint64_t seed = 1; seed <= 40; ++seed)
        {
            SCOPED_TRACE(testing::Message() << "risk aversion " << risk_aversion << " seed " << seed);
            std::size_t n = 2 + static_cast<std::size_t>(seed % 7);
            Data data = random_data(n, seed);
            Spec spec = utility(risk_aversion);
            spec.leverage = 1.0;
            Result result = solve_checked(spec, data, Scales{risk_aversion, 1.0});
            ASSERT_EQ(result.status, SolveStatus::Optimal);
            for (double weight : result.weights)
            {
                EXPECT_GE(weight, -1e-12);
            }
        }
    }
}

TEST(Importance, ZeroIsMinVariance)
{
    Data data = random_data(10, 2);
    Result result = solve_checked(long_only(importance(0.0)), data, Scales{2.0, 0.0});
    expect_identical(result, po::solve(long_only(Spec{}), data).value());
}

TEST(Importance, MatchesUtility)
{
    constexpr double kLambda = 0.3;
    Data data = random_data(12, 5);
    Spec spec = long_only(importance(kLambda));
    Result result = solve_checked(spec, data, Scales{2.0 * (1.0 - kLambda), kLambda});
    double delta = 2.0 * (1.0 - kLambda) / kLambda;
    Result reference = solve_checked(long_only(utility(delta)), data, Scales{delta, 1.0});
    expect_near_relative(result.weights, reference.weights, 1e-12);
    std::vector<double> expected_multipliers = scaled(reference.bound_multipliers, kLambda);
    expected_multipliers.push_back(kLambda * reference.budget_multiplier);
    std::vector<double> multipliers = result.bound_multipliers;
    multipliers.push_back(result.budget_multiplier);
    expect_near_relative(multipliers, expected_multipliers, 1e-10);
}

TEST(Importance, NearOne)
{
    constexpr double kLambda = 0.999999;
    Data data = random_data(10, 7);
    Result result = solve_checked(long_only(importance(kLambda)), data, Scales{2.0 * (1.0 - kLambda), kLambda});
    EXPECT_EQ(result.status, SolveStatus::Optimal);
}

TEST(Kelly, WithoutCashIsUtilityOne)
{
    Data data = random_data(10, 8);
    Result unconstrained = solve_checked(kelly(), data, Scales{1.0, 1.0});
    expect_near_relative(unconstrained.weights, po::solve(utility(1.0), data).value().weights, 1e-12);
    Result bounded = solve_checked(long_only(kelly()), data, Scales{1.0, 1.0});
    expect_near_relative(bounded.weights, po::solve(long_only(utility(1.0)), data).value().weights, 1e-12);
}

TEST(Kelly, WithCashIsUtilityOfRate)
{
    constexpr std::size_t n = 9;
    constexpr double kDelta = 1.0 / (1.0 + kRiskFree);
    Data data = random_data(n, 10);
    data.set_risk_free(kRiskFree);
    Result unconstrained = solve_checked(kelly(), data, Scales{kDelta, 1.0});
    expect_near_relative(unconstrained.weights, po::solve(utility(kDelta), data).value().weights, 1e-12);
    std::vector<double> expected = scaled(covariance_solve(data, excess_mu(data)), 1.0 + kRiskFree);
    expect_near_relative(unconstrained.weights, expected, 1e-12);
    Result bounded = solve_checked(long_only(kelly()), data, Scales{kDelta, 1.0});
    expect_near_relative(bounded.weights, po::solve(long_only(utility(kDelta)), data).value().weights, 1e-12);
}

} // namespace
