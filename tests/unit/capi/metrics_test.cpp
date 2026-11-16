#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <portopt/portopt.h>
#include <vector>

#include "capi_scope.hpp"

namespace
{

using po::testing::DataPtr;

const std::vector<double> kMu = {0.1, 0.2};
const std::vector<double> kCovariance = {0.04, 0.01, 0.01, 0.09};

TEST(CapiMetrics, ReturnAndVolatility)
{
    DataPtr data = po::testing::make_data(2, kMu, kCovariance);
    const double weights[2] = {0.5, 0.5};
    double value = 0.0;
    ASSERT_EQ(po_portfolio_return(data.get(), weights, 2, &value), PO_OK);
    EXPECT_NEAR(value, 0.15, 1e-15);
    ASSERT_EQ(po_portfolio_volatility(data.get(), weights, 2, &value), PO_OK);
    EXPECT_NEAR(value, std::sqrt(0.0375), 1e-15);
}

TEST(CapiMetrics, CashEntersTheReturnOnly)
{
    DataPtr data = po::testing::make_data(2, kMu, kCovariance);
    ASSERT_EQ(po_data_set_risk_free(data.get(), 0.02), PO_OK);
    const double weights[3] = {0.4, 0.4, 0.2};
    double value = 0.0;
    ASSERT_EQ(po_portfolio_return(data.get(), weights, 3, &value), PO_OK);
    EXPECT_NEAR(value, 0.124, 1e-15);
    ASSERT_EQ(po_portfolio_volatility(data.get(), weights, 3, &value), PO_OK);
    EXPECT_NEAR(value, std::sqrt(0.024), 1e-15);
}

TEST(CapiMetrics, WrongWeightCountIsInvalid)
{
    DataPtr data = po::testing::make_data(2, kMu, kCovariance);
    const double weights[3] = {0.4, 0.4, 0.2};
    double value = 1.0;
    EXPECT_EQ(po_portfolio_return(data.get(), weights, 3, &value), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "n must be 2, got 3");
    EXPECT_EQ(po_portfolio_volatility(data.get(), weights, 1, &value), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "n must be 2, got 1");
    ASSERT_EQ(po_data_set_risk_free(data.get(), 0.02), PO_OK);
    EXPECT_EQ(po_portfolio_return(data.get(), weights, 2, &value), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "n must be 3, got 2");
    EXPECT_EQ(value, 1.0);
}

TEST(CapiMetrics, MissingMuIsInvalid)
{
    DataPtr data = po::testing::make_data(2, {}, kCovariance);
    const double weights[2] = {0.5, 0.5};
    double value = 1.0;
    EXPECT_EQ(po_portfolio_return(data.get(), weights, 2, &value), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "the return needs mu");
    EXPECT_EQ(po_portfolio_volatility(data.get(), weights, 2, &value), PO_OK);
}

TEST(CapiMetrics, MissingCovarianceIsInvalid)
{
    DataPtr data = po::testing::make_data(2, kMu, {});
    const double weights[2] = {0.5, 0.5};
    double value = 1.0;
    EXPECT_EQ(po_portfolio_volatility(data.get(), weights, 2, &value), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "the volatility needs the covariance");
    EXPECT_EQ(po_portfolio_return(data.get(), weights, 2, &value), PO_OK);
}

DataPtr data_with_returns()
{
    DataPtr data = po::testing::make_data(2, kMu, kCovariance);
    po::testing::set_returns_ok(data.get(), {0.01, 0.02, 0.03, -0.01, -0.02, 0.04}, {0.005, 0.01, 0.0});
    return data;
}

TEST(CapiMetrics, UtilityAndTrackingError)
{
    DataPtr data = data_with_returns();
    const double weights[2] = {0.5, 0.5};
    double value = 0.0;
    ASSERT_EQ(po_utility(data.get(), weights, 2, 4.0, &value), PO_OK);
    EXPECT_NEAR(value, 0.15 - 2.0 * 0.0375, 1e-15);
    ASSERT_EQ(po_tracking_error(data.get(), weights, 2, &value), PO_OK);
    EXPECT_NEAR(value, std::sqrt(2.0e-4 / 3.0), 1e-15);

    ASSERT_EQ(po_data_set_risk_free(data.get(), 0.02), PO_OK);
    const double with_cash[3] = {0.4, 0.4, 0.2};
    ASSERT_EQ(po_utility(data.get(), with_cash, 3, 2.0, &value), PO_OK);
    EXPECT_NEAR(value, 0.124 - 0.024, 1e-15);
    ASSERT_EQ(po_tracking_error(data.get(), with_cash, 3, &value), PO_OK);
    EXPECT_NEAR(value, std::sqrt(2.69e-4 / 3.0), 1e-15);
}

TEST(CapiMetrics, TrackingErrorWithoutReturnsIsInvalid)
{
    DataPtr data = po::testing::make_data(2, kMu, kCovariance);
    const double weights[2] = {0.5, 0.5};
    double value = 1.0;
    EXPECT_EQ(po_tracking_error(data.get(), weights, 2, &value), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "the tracking error needs returns");
    EXPECT_EQ(value, 1.0);
}

TEST(CapiMetrics, UtilityNeedsMuCovarianceAndPositiveDelta)
{
    const double weights[2] = {0.5, 0.5};
    double value = 1.0;
    DataPtr without_mu = po::testing::make_data(2, {}, kCovariance);
    EXPECT_EQ(po_utility(without_mu.get(), weights, 2, 2.0, &value), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "the utility needs mu");
    DataPtr without_covariance = po::testing::make_data(2, kMu, {});
    EXPECT_EQ(po_utility(without_covariance.get(), weights, 2, 2.0, &value), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "the utility needs the covariance");

    DataPtr data = po::testing::make_data(2, kMu, kCovariance);
    for (double delta : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                         std::numeric_limits<double>::infinity()})
    {
        SCOPED_TRACE(delta);
        EXPECT_EQ(po_utility(data.get(), weights, 2, delta, &value), PO_INVALID_ARGUMENT);
        EXPECT_STREQ(po_last_error(), "delta must be positive and finite");
    }
    EXPECT_EQ(value, 1.0);
}

TEST(CapiMetrics, NewMetricsCheckWeightsAndNulls)
{
    DataPtr data = data_with_returns();
    const double weights[2] = {0.5, 0.5};
    double value = 1.0;
    EXPECT_EQ(po_utility(nullptr, weights, 2, 2.0, &value), PO_INVALID_HANDLE);
    EXPECT_EQ(po_tracking_error(nullptr, weights, 2, &value), PO_INVALID_HANDLE);
    EXPECT_EQ(po_utility(data.get(), nullptr, 2, 2.0, &value), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_tracking_error(data.get(), nullptr, 2, &value), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_utility(data.get(), weights, 2, 2.0, nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_tracking_error(data.get(), weights, 2, nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_utility(data.get(), weights, 3, 2.0, &value), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "n must be 2, got 3");
    EXPECT_EQ(po_tracking_error(data.get(), weights, 1, &value), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "n must be 2, got 1");
    EXPECT_EQ(value, 1.0);
}

TEST(CapiMetrics, SharpeAndDiversification)
{
    DataPtr data = po::testing::make_data(2, kMu, kCovariance);
    const double weights[2] = {0.5, 0.5};
    double value = 0.0;
    ASSERT_EQ(po_sharpe_ratio(data.get(), weights, 2, &value), PO_OK);
    EXPECT_NEAR(value, 0.15 / std::sqrt(0.0375), 1e-14);
    ASSERT_EQ(po_diversification_ratio(data.get(), weights, 2, &value), PO_OK);
    EXPECT_NEAR(value, 0.25 / std::sqrt(0.0375), 1e-14);

    ASSERT_EQ(po_data_set_risk_free(data.get(), 0.02), PO_OK);
    const double with_cash[3] = {0.4, 0.4, 0.2};
    ASSERT_EQ(po_sharpe_ratio(data.get(), with_cash, 3, &value), PO_OK);
    EXPECT_NEAR(value, 0.104 / std::sqrt(0.024), 1e-14);
    ASSERT_EQ(po_diversification_ratio(data.get(), with_cash, 3, &value), PO_OK);
    EXPECT_NEAR(value, 0.2 / std::sqrt(0.024), 1e-14);

    const double all_cash[3] = {0.0, 0.0, 1.0};
    value = 1.0;
    ASSERT_EQ(po_sharpe_ratio(data.get(), all_cash, 3, &value), PO_OK);
    EXPECT_TRUE(std::isnan(value));
    value = 1.0;
    ASSERT_EQ(po_diversification_ratio(data.get(), all_cash, 3, &value), PO_OK);
    EXPECT_TRUE(std::isnan(value));
}

TEST(CapiMetrics, RatiosNeedTheirInputs)
{
    const double weights[2] = {0.5, 0.5};
    double value = 1.0;
    DataPtr without_mu = po::testing::make_data(2, {}, kCovariance);
    EXPECT_EQ(po_sharpe_ratio(without_mu.get(), weights, 2, &value), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "the Sharpe ratio needs mu");
    EXPECT_EQ(po_diversification_ratio(without_mu.get(), weights, 2, &value), PO_OK);

    DataPtr without_covariance = po::testing::make_data(2, kMu, {});
    EXPECT_EQ(po_sharpe_ratio(without_covariance.get(), weights, 2, &value), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "the Sharpe ratio needs the covariance");
    value = 1.0;
    EXPECT_EQ(po_diversification_ratio(without_covariance.get(), weights, 2, &value), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "the diversification ratio needs the covariance");
    EXPECT_EQ(value, 1.0);
}

TEST(CapiMetrics, RatiosCheckWeightsAndNulls)
{
    DataPtr data = po::testing::make_data(2, kMu, kCovariance);
    const double weights[2] = {0.5, 0.5};
    double value = 1.0;
    EXPECT_EQ(po_sharpe_ratio(nullptr, weights, 2, &value), PO_INVALID_HANDLE);
    EXPECT_EQ(po_diversification_ratio(nullptr, weights, 2, &value), PO_INVALID_HANDLE);
    EXPECT_EQ(po_sharpe_ratio(data.get(), nullptr, 2, &value), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_diversification_ratio(data.get(), nullptr, 2, &value), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_sharpe_ratio(data.get(), weights, 2, nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_diversification_ratio(data.get(), weights, 2, nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_sharpe_ratio(data.get(), weights, 3, &value), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "n must be 2, got 3");
    EXPECT_EQ(po_diversification_ratio(data.get(), weights, 1, &value), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "n must be 2, got 1");
    EXPECT_EQ(value, 1.0);
}

TEST(CapiMetrics, NullArguments)
{
    DataPtr data = po::testing::make_data(2, kMu, kCovariance);
    const double weights[2] = {0.5, 0.5};
    double value = 0.0;
    EXPECT_EQ(po_portfolio_return(nullptr, weights, 2, &value), PO_INVALID_HANDLE);
    EXPECT_EQ(po_portfolio_volatility(nullptr, weights, 2, &value), PO_INVALID_HANDLE);
    EXPECT_EQ(po_portfolio_return(data.get(), nullptr, 2, &value), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_portfolio_volatility(data.get(), nullptr, 2, &value), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_portfolio_return(data.get(), weights, 2, nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_portfolio_volatility(data.get(), weights, 2, nullptr), PO_INVALID_ARGUMENT);
}

} // namespace
