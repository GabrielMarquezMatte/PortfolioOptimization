#include <gtest/gtest.h>

#include <cmath>
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
