#include "estimation/moments.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <limits>
#include <vector>

#include "random_covariance.hpp"

namespace
{

using po::estimation::sample_moments;

void expect_invalid(const std::vector<double>& returns, std::size_t t, std::size_t n, const char* message)
{
    auto moments = sample_moments(returns, t, n);
    ASSERT_FALSE(moments.has_value());
    EXPECT_EQ(moments.error().kind, po::FailureKind::InvalidArgument);
    EXPECT_EQ(moments.error().message, message);
}

TEST(Moments, SmallPanelByHand)
{
    std::vector<double> returns = {0.01, 0.02, 0.03, -0.01, -0.01, 0.02};
    auto moments = sample_moments(returns, 3, 2);
    ASSERT_TRUE(moments.has_value());
    std::vector<double> mu = {0.01, 0.01};
    std::vector<double> covariance = {0.0004, -0.0003, -0.0003, 0.0003};
    for (std::size_t i = 0; i < 2; ++i)
    {
        EXPECT_NEAR(moments->mu[i], mu[i], 1e-18);
    }
    for (std::size_t i = 0; i < 4; ++i)
    {
        EXPECT_NEAR(moments->covariance[i], covariance[i], 1e-18);
    }
}

TEST(Moments, CovarianceIsExactlySymmetric)
{
    std::size_t n = 4;
    std::vector<double> returns = po::testing::random_vector(50 * n, 7, -0.05, 0.05);
    auto moments = sample_moments(returns, 50, n);
    ASSERT_TRUE(moments.has_value());
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            EXPECT_EQ(moments->covariance[i * n + j], moments->covariance[j * n + i]);
        }
    }
}

TEST(Moments, OneObservationIsRejected)
{
    expect_invalid({0.01, 0.02}, 1, 2, "returns need at least two observations");
}

TEST(Moments, NonFiniteIsRejected)
{
    expect_invalid({0.01, 0.02, std::numeric_limits<double>::quiet_NaN(), 0.01}, 2, 2,
                   "returns is not finite at (1, 0)");
}

TEST(Moments, ZeroAssetsIsRejected)
{
    expect_invalid({}, 3, 0, "asset count must be positive");
}

} // namespace
