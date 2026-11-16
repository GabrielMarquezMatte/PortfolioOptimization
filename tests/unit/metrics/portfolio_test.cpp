#include "metrics/portfolio.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

#include "linalg/cholesky.hpp"
#include "problems/solve.hpp"
#include "random_covariance.hpp"

namespace
{

po::Data two_assets()
{
    return po::Data::create(2, std::vector<double>{0.1, 0.2}, std::vector<double>{0.04, 0.01, 0.01, 0.09}).value();
}

TEST(Metrics, ReturnAndVolatility)
{
    po::Data data = two_assets();
    std::vector<double> weights = {0.5, 0.5};
    EXPECT_NEAR(po::metrics::portfolio_return(data, weights), 0.15, 1e-15);
    EXPECT_NEAR(po::metrics::portfolio_volatility(data, weights), std::sqrt(0.0375), 1e-15);
}

TEST(Metrics, CashEntersTheReturnOnly)
{
    po::Data data = two_assets();
    data.set_risk_free(0.02);
    std::vector<double> weights = {0.4, 0.4, 0.2};
    EXPECT_NEAR(po::metrics::portfolio_return(data, weights), 0.124, 1e-15);
    EXPECT_NEAR(po::metrics::portfolio_volatility(data, weights), std::sqrt(0.024), 1e-15);
}

TEST(Metrics, AllCashHasZeroVolatility)
{
    po::Data data = two_assets();
    data.set_risk_free(0.02);
    std::vector<double> weights = {0.0, 0.0, 1.0};
    EXPECT_EQ(po::metrics::portfolio_return(data, weights), 0.02);
    EXPECT_EQ(po::metrics::portfolio_volatility(data, weights), 0.0);
}

constexpr std::size_t kPeriods = 3;

po::Data two_assets_with_returns()
{
    po::Data data = two_assets();
    std::vector<double> returns = {0.01, 0.02, 0.03, -0.01, -0.02, 0.04};
    std::vector<double> benchmark = {0.005, 0.01, 0.0};
    EXPECT_TRUE(data.set_returns(returns, benchmark, kPeriods).has_value());
    return data;
}

TEST(Metrics, UtilityAndTrackingError)
{
    po::Data data = two_assets_with_returns();
    std::vector<double> weights = {0.5, 0.5};
    EXPECT_NEAR(po::metrics::utility(data, weights, 4.0), 0.15 - 2.0 * 0.0375, 1e-15);
    EXPECT_NEAR(po::metrics::tracking_error(data, weights), std::sqrt(2.0e-4 / 3.0), 1e-15);
}

TEST(Metrics, UtilityAndTrackingErrorWithCashIncludeTheCashReturn)
{
    po::Data data = two_assets_with_returns();
    data.set_risk_free(0.02);
    std::vector<double> weights = {0.4, 0.4, 0.2};
    EXPECT_NEAR(po::metrics::utility(data, weights, 2.0), 0.124 - 0.024, 1e-15);
    EXPECT_NEAR(po::metrics::tracking_error(data, weights), std::sqrt(2.69e-4 / 3.0), 1e-15);
}

double root_of_excess_objective(const po::Data& data, const std::vector<double>& weights)
{
    std::size_t n = data.asset_count();
    double rate = data.risk_free().value_or(0.0);
    double objective = 0.0;
    for (std::size_t t = 0; t < data.observation_count(); ++t)
    {
        double fitted = rate;
        for (std::size_t i = 0; i < n; ++i)
        {
            fitted += (data.returns()[t * n + i] - rate) * weights[i];
        }
        double gap = fitted - data.benchmark()[t];
        objective += gap * gap;
    }
    return std::sqrt(objective / static_cast<double>(data.observation_count()));
}

TEST(Metrics, TrackingErrorAtOptimumIsRootOfObjective)
{
    constexpr std::size_t kAssets = 3;
    constexpr std::size_t kObservations = 20;
    for (bool with_cash : {false, true})
    {
        SCOPED_TRACE(with_cash);
        po::Data data = po::Data::create(kAssets, {}, po::testing::random_covariance(kAssets, 5)).value();
        std::vector<double> returns = po::testing::random_vector(kObservations * kAssets, 51, -0.05, 0.05);
        std::vector<double> benchmark = po::testing::random_vector(kObservations, 52, -0.03, 0.03);
        ASSERT_TRUE(data.set_returns(returns, benchmark, kObservations).has_value());
        if (with_cash)
        {
            data.set_risk_free(0.004);
        }
        po::Spec spec;
        spec.objective = po::Objective::TrackingError;
        spec.lower = 0.0;
        po::Result result = po::solve(spec, data).value();
        ASSERT_EQ(result.status, po::SolveStatus::Optimal);
        double expected = root_of_excess_objective(data, result.weights);
        EXPECT_NEAR(po::metrics::tracking_error(data, result.weights), expected, 1e-15);
    }
}

TEST(Metrics, SharpeAndDiversification)
{
    po::Data data = two_assets();
    std::vector<double> weights = {0.5, 0.5};
    EXPECT_NEAR(po::metrics::sharpe_ratio(data, weights), 0.15 / std::sqrt(0.0375), 1e-14);
    EXPECT_NEAR(po::metrics::diversification_ratio(data, weights), 0.25 / std::sqrt(0.0375), 1e-14);

    data.set_risk_free(0.02);
    std::vector<double> with_cash = {0.4, 0.4, 0.2};
    EXPECT_NEAR(po::metrics::sharpe_ratio(data, with_cash), 0.104 / std::sqrt(0.024), 1e-14);
    EXPECT_NEAR(po::metrics::diversification_ratio(data, with_cash), 0.2 / std::sqrt(0.024), 1e-14);

    std::vector<double> all_cash = {0.0, 0.0, 1.0};
    EXPECT_TRUE(std::isnan(po::metrics::sharpe_ratio(data, all_cash)));
    EXPECT_TRUE(std::isnan(po::metrics::diversification_ratio(data, all_cash)));
}

po::Spec sharpe_spec()
{
    po::Spec spec;
    spec.objective = po::Objective::MaxSharpe;
    return spec;
}

// -f* at the unconstrained tangency: sqrt(mu~' Sigma^-1 mu~).
double tangency_ratio(const po::Data& data)
{
    std::size_t n = data.asset_count();
    double rate = data.risk_free().value_or(0.0);
    std::vector<double> excess(data.mu().begin(), data.mu().end());
    for (double& value : excess)
    {
        value -= rate;
    }
    std::vector<double> solved = excess;
    std::vector<double> factor = po::linalg::cholesky(data.covariance(), n).value();
    po::linalg::solve_lower(factor, n, solved);
    po::linalg::solve_lower_transposed(factor, n, solved);
    double total = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        total += excess[i] * solved[i];
    }
    return std::sqrt(total);
}

TEST(Metrics, SharpeAtOptimumIsMaximal)
{
    constexpr std::size_t kAssets = 5;
    for (bool with_cash : {false, true})
    {
        SCOPED_TRACE(with_cash);
        po::Data data = po::Data::create(kAssets, po::testing::random_vector(kAssets, 61, 0.03, 0.15),
                                         po::testing::random_covariance(kAssets, 6))
                            .value();
        if (with_cash)
        {
            data.set_risk_free(0.02);
        }
        po::Result tangency = po::solve(sharpe_spec(), data).value();
        ASSERT_EQ(tangency.status, po::SolveStatus::Optimal);
        double expected = tangency_ratio(data);
        EXPECT_NEAR(po::metrics::sharpe_ratio(data, tangency.weights), expected, 1e-10 * expected);

        po::Spec long_only = sharpe_spec();
        long_only.lower = 0.0;
        po::Result bounded = po::solve(long_only, data).value();
        ASSERT_EQ(bounded.status, po::SolveStatus::Optimal);
        double best = po::metrics::sharpe_ratio(data, bounded.weights);
        EXPECT_LE(best, expected * (1.0 + 1e-10));
        for (std::size_t asset = 0; asset < kAssets; ++asset)
        {
            std::vector<double> single(data.weight_count(), 0.0);
            single[asset] = 1.0;
            EXPECT_GE(best, po::metrics::sharpe_ratio(data, single) * (1.0 - 1e-10));
        }
    }
}

} // namespace
