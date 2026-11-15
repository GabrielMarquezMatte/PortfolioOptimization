#include "problems/tracking_error.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <vector>

#include "kkt.hpp"
#include "linalg/cholesky.hpp"
#include "problems/solve.hpp"
#include "qp/goldfarb_idnani.hpp"
#include "random_covariance.hpp"
#include "user_kkt.hpp"

namespace
{

using po::Data;
using po::Failure;
using po::FailureKind;
using po::Objective;
using po::Result;
using po::SolveStatus;
using po::Spec;
using po::testing::dot;

constexpr double kRiskFree = 0.001;
constexpr const char* kRankMessage = "returns minus the risk-free rate do not have full column rank";

struct Series
{
    std::vector<double> returns;
    std::vector<double> benchmark;
    std::size_t periods;
};

Spec tracking_error()
{
    Spec spec;
    spec.objective = Objective::TrackingError;
    return spec;
}

Spec long_only(Spec spec)
{
    spec.lower = 0.0;
    return spec;
}

std::vector<double> random_returns(std::size_t periods, std::size_t n, std::uint64_t seed)
{
    return po::testing::random_vector(periods * n, seed, -0.05, 0.05);
}

std::vector<double> weights_summing_to_one(std::size_t n, std::uint64_t seed)
{
    std::vector<double> weights = po::testing::random_vector(n, seed, -0.2, 1.0);
    double total = 0.0;
    for (double weight : weights)
    {
        total += weight;
    }
    for (double& weight : weights)
    {
        weight /= total;
    }
    return weights;
}

// R w for the risky weights, row by row, written independently of the builder.
std::vector<double> product(const std::vector<double>& returns, std::size_t periods,
                             std::span<const double> weights)
{
    std::size_t n = weights.size();
    std::vector<double> values(periods);
    for (std::size_t t = 0; t < periods; ++t)
    {
        values[t] = dot(std::span<const double>(returns).subspan(t * n, n), weights);
    }
    return values;
}

Data data_with(const Series& series, std::size_t n, bool cash)
{
    Data data = Data::create(n, {}, {}).value();
    EXPECT_TRUE(data.set_returns(series.returns, series.benchmark, series.periods).has_value());
    if (cash)
    {
        data.set_risk_free(kRiskFree);
    }
    return data;
}

// grad f = 2 R~'(R~ w + rf 1 - b), with R~ = R - rf 1 1' and rf = 0 without cash.
std::vector<double> objective_gradient(const Data& data, std::span<const double> weights)
{
    std::size_t n = data.asset_count();
    std::size_t periods = data.observation_count();
    double rate = data.risk_free().value_or(0.0);
    std::vector<double> gradient(n, 0.0);
    for (std::size_t t = 0; t < periods; ++t)
    {
        double residual = rate - data.benchmark()[t];
        for (std::size_t i = 0; i < n; ++i)
        {
            residual += (data.returns()[t * n + i] - rate) * weights[i];
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            gradient[i] += 2.0 * (data.returns()[t * n + i] - rate) * residual;
        }
    }
    return gradient;
}

double tracking_error_value(const Data& data, std::span<const double> weights)
{
    std::size_t n = data.asset_count();
    double rate = data.risk_free().value_or(0.0);
    double total = 0.0;
    for (std::size_t t = 0; t < data.observation_count(); ++t)
    {
        double residual = rate - data.benchmark()[t];
        for (std::size_t i = 0; i < n; ++i)
        {
            residual += (data.returns()[t * n + i] - rate) * weights[i];
        }
        total += residual * residual;
    }
    return total;
}

Result solve_checked(const Spec& spec, const Data& data, double tolerance = 1e-10)
{
    po::problems::Translation translation = po::problems::tracking_error_problem(spec, data).value();
    po::testing::expect_kkt(translation.problem, po::qp::solve(translation.problem), tolerance);
    Result result = po::solve(spec, data).value();
    po::testing::expect_user_kkt(spec, data, result, objective_gradient(data, result.weights), tolerance);
    return result;
}

void expect_rank_failure(const Spec& spec, const Data& data)
{
    std::expected<po::problems::Translation, Failure> translation =
        po::problems::tracking_error_problem(spec, data);
    ASSERT_FALSE(translation.has_value());
    EXPECT_EQ(translation.error().kind, FailureKind::InvalidArgument);
    EXPECT_EQ(translation.error().message, kRankMessage);
    std::expected<Result, Failure> result = po::solve(spec, data);
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error().kind, FailureKind::InvalidArgument);
    EXPECT_EQ(result.error().message, kRankMessage);
}

std::vector<double> solve_with_normal_matrix(const std::vector<double>& returns, std::size_t periods,
                                             std::size_t n, std::vector<double> rhs)
{
    std::vector<double> normal(n * n, 0.0);
    for (std::size_t t = 0; t < periods; ++t)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            for (std::size_t j = 0; j < n; ++j)
            {
                normal[i * n + j] += returns[t * n + i] * returns[t * n + j];
            }
        }
    }
    std::vector<double> factor = po::linalg::cholesky(normal, n).value();
    po::linalg::solve_lower(factor, n, rhs);
    po::linalg::solve_lower_transposed(factor, n, rhs);
    return rhs;
}

std::vector<double> transposed_product(const std::vector<double>& returns, std::size_t periods, std::size_t n,
                                       const std::vector<double>& values)
{
    std::vector<double> result(n, 0.0);
    for (std::size_t t = 0; t < periods; ++t)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            result[i] += returns[t * n + i] * values[t];
        }
    }
    return result;
}

TEST(TrackingError, RecoversExactReplication)
{
    constexpr std::size_t periods = 60;
    constexpr std::size_t n = 5;
    std::vector<double> returns = random_returns(periods, n, 1);
    std::vector<double> expected = weights_summing_to_one(n, 2);
    Series series{returns, product(returns, periods, expected), periods};
    Data data = data_with(series, n, false);
    Result result = solve_checked(tracking_error(), data);
    ASSERT_EQ(result.weights.size(), n);
    for (std::size_t i = 0; i < n; ++i)
    {
        EXPECT_NEAR(result.weights[i], expected[i], 1e-10);
    }
    EXPECT_NEAR(tracking_error_value(data, result.weights), 0.0, 1e-24);
}

TEST(TrackingError, MatchesConstrainedLeastSquares)
{
    constexpr std::size_t periods = 60;
    constexpr std::size_t n = 6;
    std::vector<double> returns = random_returns(periods, n, 3);
    std::vector<double> benchmark = po::testing::random_vector(periods, 4, -0.05, 0.05);
    Data data = data_with(Series{returns, benchmark, periods}, n, false);
    std::vector<double> ones(n, 1.0);
    std::vector<double> inverse_ones = solve_with_normal_matrix(returns, periods, n, ones);
    std::vector<double> inverse_moment =
        solve_with_normal_matrix(returns, periods, n, transposed_product(returns, periods, n, benchmark));
    double budget_share = (1.0 - dot(ones, inverse_moment)) / dot(ones, inverse_ones);
    Result result = solve_checked(tracking_error(), data);
    for (std::size_t i = 0; i < n; ++i)
    {
        EXPECT_NEAR(result.weights[i], inverse_moment[i] + budget_share * inverse_ones[i], 1e-10);
    }
}

TEST(TrackingError, CashReplication)
{
    constexpr std::size_t periods = 60;
    constexpr std::size_t n = 5;
    std::vector<double> returns = random_returns(periods, n, 5);
    std::vector<double> expected = po::testing::random_vector(n, 6, -0.3, 0.6);
    std::vector<double> excess_returns = returns;
    for (double& value : excess_returns)
    {
        value -= kRiskFree;
    }
    std::vector<double> benchmark = product(excess_returns, periods, expected);
    for (double& value : benchmark)
    {
        value += kRiskFree;
    }
    Data data = data_with(Series{returns, benchmark, periods}, n, true);
    Result result = solve_checked(tracking_error(), data);
    ASSERT_EQ(result.weights.size(), n + 1);
    double cash = 1.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        EXPECT_NEAR(result.weights[i], expected[i], 1e-10);
        cash -= expected[i];
    }
    EXPECT_NEAR(result.weights[n], cash, 1e-10);
    EXPECT_NEAR(tracking_error_value(data, result.weights), 0.0, 1e-24);
}

TEST(TrackingError, CashWithZeroRateIsNotSingular)
{
    constexpr std::size_t periods = 40;
    constexpr std::size_t n = 4;
    Series series{random_returns(periods, n, 7), po::testing::random_vector(periods, 8, -0.05, 0.05),
                  periods};
    Data data = data_with(series, n, true);
    data.set_risk_free(0.0);
    Result result = solve_checked(tracking_error(), data);
    EXPECT_EQ(result.status, SolveStatus::Optimal);
    ASSERT_EQ(result.weights.size(), n + 1);
}

TEST(TrackingError, TargetReturnAndLongOnlyKkt)
{
    constexpr std::size_t periods = 50;
    constexpr std::size_t n = 8;
    for (std::uint64_t seed = 1; seed <= 10; ++seed)
    {
        SCOPED_TRACE(seed);
        Series series{random_returns(periods, n, seed),
                      po::testing::random_vector(periods, seed + 100, -0.05, 0.05), periods};
        std::vector<double> mu = po::testing::random_vector(n, seed + 200, 0.03, 0.15);
        Data data = Data::create(n, mu, {}).value();
        ASSERT_TRUE(data.set_returns(series.returns, series.benchmark, periods).has_value());
        Spec spec = long_only(tracking_error());
        spec.target_return = 0.09;
        Result result = solve_checked(spec, data);
        EXPECT_FALSE(std::isnan(result.target_return_multiplier));
        for (double weight : result.weights)
        {
            EXPECT_GE(weight, -1e-14);
        }
    }
}

TEST(TrackingError, CashTargetReturnAndLongOnlyKkt)
{
    constexpr std::size_t periods = 50;
    constexpr std::size_t n = 8;
    for (std::uint64_t seed = 1; seed <= 10; ++seed)
    {
        SCOPED_TRACE(seed);
        Series series{random_returns(periods, n, seed + 300),
                      po::testing::random_vector(periods, seed + 400, -0.05, 0.05), periods};
        std::vector<double> mu = po::testing::random_vector(n, seed + 500, 0.03, 0.15);
        Data data = Data::create(n, mu, {}).value();
        ASSERT_TRUE(data.set_returns(series.returns, series.benchmark, periods).has_value());
        data.set_risk_free(kRiskFree);
        Spec spec = long_only(tracking_error());
        spec.target_return = 0.06;
        Result result = solve_checked(spec, data);
        EXPECT_GE(result.weights[n], -1e-14);
    }
}

TEST(TrackingError, RankDeficientIsInvalid)
{
    constexpr std::size_t n = 5;
    Series fewer{random_returns(3, n, 9), po::testing::random_vector(3, 10, -0.05, 0.05), 3};
    expect_rank_failure(tracking_error(), data_with(fewer, n, false));

    constexpr std::size_t periods = 30;
    constexpr std::size_t assets = 4;
    std::vector<double> duplicated = random_returns(periods, assets, 11);
    for (std::size_t t = 0; t < periods; ++t)
    {
        duplicated[t * assets + 3] = duplicated[t * assets + 1];
    }
    std::vector<double> benchmark = po::testing::random_vector(periods, 12, -0.05, 0.05);
    expect_rank_failure(tracking_error(), data_with(Series{duplicated, benchmark, periods}, assets, false));

    std::vector<double> constant = random_returns(periods, assets, 13);
    for (std::size_t t = 0; t < periods; ++t)
    {
        constant[t * assets + 2] = kRiskFree;
    }
    expect_rank_failure(tracking_error(), data_with(Series{constant, benchmark, periods}, assets, true));
    Data unrelated_to_cash = data_with(Series{constant, benchmark, periods}, assets, false);
    EXPECT_TRUE(po::solve(tracking_error(), unrelated_to_cash).has_value());
}

TEST(TrackingError, WithoutCovariance)
{
    constexpr std::size_t periods = 30;
    constexpr std::size_t n = 3;
    Data data = Data::create(n, {}, {}).value();
    ASSERT_FALSE(data.has_covariance());
    std::vector<double> returns = random_returns(periods, n, 14);
    std::vector<double> benchmark = po::testing::random_vector(periods, 15, -0.05, 0.05);
    ASSERT_TRUE(data.set_returns(returns, benchmark, periods).has_value());
    std::expected<Result, Failure> result = po::solve(tracking_error(), data);
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->status, SolveStatus::Optimal);
}

} // namespace
