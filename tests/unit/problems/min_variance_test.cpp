#include "problems/mean_variance.hpp"

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
#include "qp/goldfarb_idnani.hpp"
#include "random_covariance.hpp"
#include "solve/solve.hpp"
#include "user_kkt.hpp"

namespace
{

using po::Data;
using po::Result;
using po::SolveStatus;
using po::Spec;
using po::testing::covariance_times;
using po::testing::dot;
using po::testing::excess_mu;

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

Result solve_checked(const Spec& spec, const Data& data, double tolerance = 1e-10)
{
    po::problems::MeanVariance scales = po::problems::mean_variance_of(spec, data);
    po::problems::Translation translation = po::problems::mean_variance_problem(spec, data, scales);
    po::testing::expect_kkt(translation.problem, po::qp::solve(translation.problem), tolerance);
    Result result = po::solve(spec, data).value();
    std::vector<double> gradient = covariance_times(data, result.weights);
    for (double& entry : gradient)
    {
        entry *= 2.0;
    }
    po::testing::expect_user_kkt(spec, data, result, po::testing::sum_of({gradient}), tolerance);
    return result;
}

double variance(const Data& data, const Result& result)
{
    std::span<const double> risky = std::span<const double>(result.weights).first(data.asset_count());
    return dot(risky, covariance_times(data, risky));
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
    EXPECT_TRUE(std::isnan(result.residual));
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

TEST(MinVariance, WeightBounds)
{
    constexpr std::size_t n = 6;
    Data data = random_data(n, 5);
    Result free = po::solve(Spec{}, data).value();
    auto largest = std::max_element(free.weights.begin(), free.weights.end());
    std::size_t capped = static_cast<std::size_t>(largest - free.weights.begin());
    Spec spec;
    spec.weight_lower.assign(n, -kInfinity);
    spec.weight_upper.assign(n, kInfinity);
    spec.weight_upper[capped] = 0.5 * *largest;
    Result result = solve_checked(spec, data);
    EXPECT_NEAR(result.weights[capped], spec.weight_upper[capped], 1e-14);
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

TEST(MinVariance, WeightBoundsSizeIncludesCash)
{
    constexpr std::size_t n = 4;
    Data data = random_data(n, 9);
    data.set_risk_free(kRiskFree);
    Spec spec;
    spec.weight_lower.assign(n, 0.0);
    std::expected<Result, po::Failure> rejected = po::solve(spec, data);
    ASSERT_FALSE(rejected.has_value());
    EXPECT_EQ(rejected.error().kind, po::FailureKind::InvalidArgument);
    EXPECT_EQ(rejected.error().message, "weight_lower has 4 entries; expected 5");
    spec.weight_lower.assign(n + 1, 0.0);
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

// The target is the largest return under 0 <= w <= 0.4, taken by filling the assets in order of mu (0.4, 0.4, 0.2
// at assets 5, 4, 0): the feasible set is that one vertex, where seven bounds, the budget and the target are tight,
// n + 1 rows. With mu_0 close to mu_4 the target is nearly a combination of the others: the last tight row
// depends on the active ones with a combination of order 1e6, and was declared infeasible on its rounding. The
// vertex is known only to that conditioning: the weights to 1e-10.
TEST(MinVariance, TargetReturnAtTheLargestAttainable)
{
    constexpr std::size_t n = 8;
    std::vector<double> mu = {0.1338745880192623,   -0.092463009515773423, 0.09538240615686816,
                              -0.010237515481265905, 0.13387509740870435,  0.13889560735303444,
                              0.024383871528689163,  0.12144568928812466};
    std::vector<double> covariance = {
        0.068304314954539092,   -0.0083312779102571507, -0.0021625493490939236, -0.0022370484942518896,
        -0.044207324519999223,  0.020080979100564653,   -0.0061652606531945164, 0.006342088736113939,
        -0.0083312779102571507, 0.013384098457535448,   0.014900733366591536,   -0.0096509372607783417,
        0.0065282611309635486,  -0.0085126147941023689, -0.0086500007671099291, 0.0020409389274801707,
        -0.0021625493490939236, 0.014900733366591536,   0.057729546284346137,   -0.018241095074620409,
        0.020424541245805238,   0.03830597102222022,    -0.020414464344774803,  -0.0097943855050318385,
        -0.0022370484942518896, -0.0096509372607783417, -0.018241095074620409,  0.026335529695318938,
        0.0095124629883779153,  -0.0078047571124178275, 0.029239597989153415,   -0.0073372615957962781,
        -0.044207324519999223,  0.0065282611309635486,  0.020424541245805238,   0.0095124629883779153,
        0.065912591044277685,   -0.0045332984878143356, 0.022902681024665139,   -0.0094537305955773639,
        0.020080979100564653,   -0.0085126147941023689, 0.03830597102222022,    -0.0078047571124178275,
        -0.0045332984878143356, 0.17511697880471597,    -0.063040406212934566,  -0.026012298670492321,
        -0.0061652606531945164, -0.0086500007671099291, -0.020414464344774803,  0.029239597989153415,
        0.022902681024665139,   -0.063040406212934566,  0.12062312245280171,    -0.012387930889907512,
        0.006342088736113939,   0.0020409389274801707,  -0.0097943855050318385, -0.0073372615957962781,
        -0.0094537305955773639, -0.026012298670492321,  -0.012387930889907512,  0.018745380207707916};
    Data data = Data::create(n, mu, covariance).value();
    Spec spec;
    spec.lower = 0.0;
    spec.upper = 0.4;
    spec.target_return = 0.13588319950854799;
    Result result = po::solve(spec, data).value();
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    expect_weights_near(result.weights, {0.2, 0.0, 0.0, 0.0, 0.4, 0.4, 0.0, 0.0}, 1e-9);
    solve_checked(spec, data);
}

TEST(MinVariance, FixedAsset)
{
    constexpr std::size_t n = 6;
    Data data = random_data(n, 12);
    Spec spec;
    spec.weight_lower.assign(n, -kInfinity);
    spec.weight_upper.assign(n, kInfinity);
    spec.weight_lower[0] = 0.25;
    spec.weight_upper[0] = 0.25;
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
