#include "metrics/portfolio.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

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

} // namespace
