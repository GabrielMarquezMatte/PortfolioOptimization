#include "estimation/bayes_stein.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "estimation/inputs.hpp"

namespace
{

using po::Data;
using po::estimation::bayes_stein;
using po::estimation::check_inputs;

const std::vector<double> kCovariance = {0.0004, 0.0001, 0.0001, 0.0009};

Data make_data(const std::vector<double>& mu, const std::vector<double>& covariance)
{
    auto data = Data::create(2, mu, covariance);
    EXPECT_TRUE(data.has_value());
    return std::move(*data);
}

std::vector<double> scaled(const std::vector<double>& x, double factor)
{
    std::vector<double> out = x;
    for (double& value : out)
    {
        value *= factor;
    }
    return out;
}

void expect_check_failure(const Data& data, bool needs_mu, const char* message)
{
    auto failure = check_inputs(data, needs_mu);
    ASSERT_TRUE(failure.has_value());
    EXPECT_EQ(failure->kind, po::FailureKind::InvalidArgument);
    EXPECT_EQ(failure->message, message);
}

TEST(BayesStein, ConstantMeanIsItsOwnTarget)
{
    std::vector<double> mu_hat = {0.001, 0.001};
    auto result = bayes_stein(make_data(mu_hat, kCovariance), 60);
    ASSERT_TRUE(result.has_value());
    EXPECT_NEAR(result->grand_mean, 0.001, 1e-18);
    EXPECT_NEAR(result->weight, 1.0, 1e-12);
    for (std::size_t i = 0; i < 2; ++i)
    {
        EXPECT_NEAR(result->mu[i], mu_hat[i], 1e-18);
    }
}

TEST(BayesStein, MatchesTheClosedForm)
{
    std::vector<double> mu_hat = {0.001, 0.003};
    std::size_t t = 60;
    double n = 2.0;
    double k = (60.0 - n - 2.0) / (60.0 - 1.0);
    // Sigma^-1 = k S^-1, with S^-1 from the explicit 2 x 2 inverse.
    double det = 0.0004 * 0.0009 - 0.0001 * 0.0001;
    double inv[2][2] = {{k * 0.0009 / det, -k * 0.0001 / det}, {-k * 0.0001 / det, k * 0.0004 / det}};
    double inv_mu0 = inv[0][0] * mu_hat[0] + inv[0][1] * mu_hat[1];
    double inv_mu1 = inv[1][0] * mu_hat[0] + inv[1][1] * mu_hat[1];
    double inv_one = inv[0][0] + inv[0][1] + inv[1][0] + inv[1][1];
    double mu_0 = (inv_mu0 + inv_mu1) / inv_one;
    double d0 = mu_hat[0] - mu_0;
    double d1 = mu_hat[1] - mu_0;
    double q = d0 * (inv[0][0] * d0 + inv[0][1] * d1) + d1 * (inv[1][0] * d0 + inv[1][1] * d1);
    double w = (n + 2.0) / ((n + 2.0) + static_cast<double>(t) * q);
    auto result = bayes_stein(make_data(mu_hat, kCovariance), t);
    ASSERT_TRUE(result.has_value());
    EXPECT_NEAR(result->grand_mean, mu_0, 1e-15);
    EXPECT_NEAR(result->weight, w, 1e-12);
    EXPECT_NEAR(result->mu[0], (1.0 - w) * mu_hat[0] + w * mu_0, 1e-15);
    EXPECT_NEAR(result->mu[1], (1.0 - w) * mu_hat[1] + w * mu_0, 1e-15);
}

TEST(BayesStein, WeightIgnoresTheUnit)
{
    std::vector<double> mu_hat = {0.001, 0.003};
    auto base = bayes_stein(make_data(mu_hat, kCovariance), 60);
    auto rescaled = bayes_stein(make_data(scaled(mu_hat, 100.0), scaled(kCovariance, 1e4)), 60);
    ASSERT_TRUE(base.has_value());
    ASSERT_TRUE(rescaled.has_value());
    EXPECT_LE(std::abs(rescaled->weight - base->weight), 1e-12 * base->weight);
}

TEST(BayesStein, WeightDependsOnTheFrequency)
{
    std::vector<double> mu_hat = {0.001, 0.003};
    auto base = bayes_stein(make_data(mu_hat, kCovariance), 60);
    auto annual = bayes_stein(make_data(scaled(mu_hat, 252.0), scaled(kCovariance, 252.0)), 60);
    ASSERT_TRUE(base.has_value());
    ASSERT_TRUE(annual.has_value());
    EXPECT_GT(std::abs(annual->weight - base->weight), 1e-3 * base->weight);
}

TEST(BayesStein, TooFewObservationsIsRejected)
{
    auto result = bayes_stein(make_data({0.001, 0.003}, kCovariance), 4);
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error().kind, po::FailureKind::InvalidArgument);
    EXPECT_EQ(result.error().message, "bayes-stein needs t > n + 2; got t = 4 for n = 2");
}

TEST(EstimationInputs, RiskFreeDataIsRejected)
{
    Data data = make_data({0.001, 0.003}, kCovariance);
    EXPECT_FALSE(check_inputs(data, true).has_value());
    data.set_risk_free(0.01);
    expect_check_failure(data, true, "data must not have a risk-free rate");
}

TEST(EstimationInputs, MissingMuIsRejected)
{
    expect_check_failure(make_data({}, kCovariance), true, "data has no mu");
    EXPECT_FALSE(check_inputs(make_data({}, kCovariance), false).has_value());
}

TEST(EstimationInputs, MissingCovarianceIsRejected)
{
    expect_check_failure(make_data({0.001, 0.003}, {}), false, "data has no covariance");
}

} // namespace
