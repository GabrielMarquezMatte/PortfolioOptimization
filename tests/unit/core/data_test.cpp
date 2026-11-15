#include "core/data.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <expected>
#include <limits>
#include <vector>

#include "random_covariance.hpp"

namespace
{

using po::Data;
using po::Failure;
using po::FailureKind;

void expect_invalid(const std::expected<Data, Failure>& data, const char* message)
{
    ASSERT_FALSE(data.has_value());
    EXPECT_EQ(data.error().kind, FailureKind::InvalidArgument);
    EXPECT_EQ(data.error().message, message);
}

TEST(Data, RejectsNonFinite)
{
    std::vector<double> mu = {0.1, std::numeric_limits<double>::quiet_NaN(), 0.2};
    std::vector<double> covariance = po::testing::random_covariance(3, 1);
    expect_invalid(Data::create(3, mu, covariance), "mu is not finite at 1");
    mu[1] = 0.05;
    covariance[1 * 3 + 2] = std::numeric_limits<double>::infinity();
    expect_invalid(Data::create(3, mu, covariance), "covariance is not finite at (1, 2)");
}

TEST(Data, RejectsAsymmetric)
{
    std::vector<double> covariance = po::testing::random_covariance(3, 2);
    covariance[1 * 3 + 2] += 1e-6;
    expect_invalid(Data::create(3, {}, covariance), "covariance is not symmetric at (1, 2)");
}

TEST(Data, RejectsIndefinite)
{
    std::vector<double> covariance = {1.0, 2.0, 2.0, 1.0};
    expect_invalid(Data::create(2, {}, covariance), "covariance is not positive definite");
}

TEST(Data, RejectsMismatchedSizes)
{
    std::vector<double> covariance = po::testing::random_covariance(3, 3);
    expect_invalid(Data::create(3, std::vector<double>{0.1, 0.2}, covariance), "mu has 2 entries; expected 3");
    expect_invalid(Data::create(2, {}, covariance), "covariance has 9 entries; expected 4");
    expect_invalid(Data::create(0, {}, {}), "asset count must be positive");
}

TEST(Data, AcceptsAbsentMuAndCovariance)
{
    Data data = Data::create(4, {}, {}).value();
    EXPECT_EQ(data.asset_count(), 4u);
    EXPECT_EQ(data.weight_count(), 4u);
    EXPECT_FALSE(data.has_mu());
    EXPECT_FALSE(data.has_covariance());
    EXPECT_TRUE(data.mu().empty());
    EXPECT_TRUE(data.covariance().empty());
    EXPECT_FALSE(data.risk_free().has_value());
    data.set_risk_free(0.02);
    EXPECT_EQ(data.risk_free(), 0.02);
    EXPECT_EQ(data.weight_count(), 5u);
}

TEST(Data, SymmetrizesWithinTolerance)
{
    constexpr std::size_t n = 3;
    std::vector<double> covariance = po::testing::random_covariance(n, 4);
    double lower = covariance[2 * n + 0];
    double upper = lower + 1e-12 * covariance[0];
    covariance[0 * n + 2] = upper;
    std::vector<double> mu = {0.1, 0.2, 0.3};
    Data data = Data::create(n, mu, covariance).value();
    EXPECT_TRUE(data.has_mu());
    EXPECT_TRUE(data.has_covariance());
    EXPECT_EQ(std::vector<double>(data.mu().begin(), data.mu().end()), mu);
    EXPECT_EQ(data.covariance()[0 * n + 2], 0.5 * (lower + upper));
    EXPECT_EQ(data.covariance()[2 * n + 0], 0.5 * (lower + upper));
    EXPECT_EQ(data.covariance()[1 * n + 1], covariance[1 * n + 1]);
}

} // namespace
