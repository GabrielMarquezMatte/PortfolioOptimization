#include "core/data.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>
#include <random>
#include <span>
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

void expect_invalid_returns(const std::expected<void, Failure>& result, const char* message)
{
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error().kind, FailureKind::InvalidArgument);
    EXPECT_EQ(result.error().message, message);
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

std::vector<double> sample_covariance(std::size_t n, std::size_t periods, std::uint64_t seed)
{
    std::mt19937_64 engine(seed);
    std::vector<double> returns(periods * n);
    for (double& value : returns)
    {
        value = 0.01 * po::testing::standard_normal(engine);
    }
    std::vector<double> mean(n, 0.0);
    for (std::size_t t = 0; t < periods; ++t)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            mean[i] += returns[t * n + i] / static_cast<double>(periods);
        }
    }
    std::vector<double> covariance(n * n, 0.0);
    for (std::size_t t = 0; t < periods; ++t)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            for (std::size_t j = 0; j < n; ++j)
            {
                covariance[i * n + j] += (returns[t * n + i] - mean[i]) * (returns[t * n + j] - mean[j]);
            }
        }
    }
    for (double& value : covariance)
    {
        value /= static_cast<double>(periods - 1);
    }
    return covariance;
}

TEST(Data, RejectsDuplicatedAsset)
{
    std::vector<double> covariance = {0.04, 0.04, 0.01, 0.04, 0.04, 0.01, 0.01, 0.01, 0.09};
    expect_invalid(Data::create(3, {}, covariance), "covariance is not positive definite");
}

// With fewer periods than assets the sample covariance is singular: its last pivots are rounding noise, of
// either sign.
TEST(Data, RejectsRankDeficientSampleCovariance)
{
    constexpr std::size_t n = 10;
    for (std::uint64_t seed = 1; seed <= 20; ++seed)
    {
        SCOPED_TRACE(seed);
        std::vector<double> covariance = sample_covariance(n, n - 1, seed);
        expect_invalid(Data::create(n, {}, covariance), "covariance is not positive definite");
    }
}

TEST(Data, StoresCovarianceFactor)
{
    constexpr std::size_t n = 4;
    std::vector<double> covariance = po::testing::random_covariance(n, 6);
    Data data = Data::create(n, {}, covariance).value();
    std::span<const double> factor = data.covariance_factor();
    ASSERT_EQ(factor.size(), n * n);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            double product = 0.0;
            for (std::size_t k = 0; k < n; ++k)
            {
                product += factor[i * n + k] * factor[j * n + k];
            }
            EXPECT_NEAR(product, covariance[i * n + j], 1e-14 * covariance[0]);
        }
    }
    EXPECT_TRUE(Data::create(n, {}, {}).value().covariance_factor().empty());
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

TEST(Data, SetReturnsCopiesAndValidates)
{
    constexpr std::size_t n = 3;
    constexpr std::size_t periods = 4;
    Data data = Data::create(n, {}, {}).value();
    EXPECT_FALSE(data.has_returns());
    EXPECT_EQ(data.observation_count(), 0u);
    EXPECT_TRUE(data.returns().empty());
    EXPECT_TRUE(data.benchmark().empty());
    std::vector<double> returns = po::testing::random_vector(periods * n, 11, -0.05, 0.05);
    std::vector<double> benchmark = po::testing::random_vector(periods, 12, -0.05, 0.05);
    ASSERT_TRUE(data.set_returns(returns, benchmark, periods).has_value());
    EXPECT_TRUE(data.has_returns());
    EXPECT_EQ(data.observation_count(), periods);
    EXPECT_EQ(std::vector<double>(data.returns().begin(), data.returns().end()), returns);
    EXPECT_EQ(std::vector<double>(data.benchmark().begin(), data.benchmark().end()), benchmark);
    returns[0] = 9.0;
    EXPECT_NE(data.returns()[0], 9.0);

    expect_invalid_returns(data.set_returns(returns, benchmark, 0), "returns need at least one observation");
    expect_invalid_returns(data.set_returns(returns, benchmark, periods - 1),
                           "returns has 12 entries; expected 9");
    std::span<const double> short_benchmark = std::span<const double>(benchmark).first(periods - 1);
    expect_invalid_returns(data.set_returns(returns, short_benchmark, periods),
                           "benchmark has 3 entries; expected 4");
    constexpr std::size_t huge = std::numeric_limits<std::size_t>::max() / 2;
    expect_invalid_returns(data.set_returns(returns, benchmark, huge), "returns size is too large");
    returns[5] = std::numeric_limits<double>::quiet_NaN();
    expect_invalid_returns(data.set_returns(returns, benchmark, periods), "returns is not finite at (1, 2)");
    returns[5] = 0.01;
    benchmark[3] = std::numeric_limits<double>::infinity();
    expect_invalid_returns(data.set_returns(returns, benchmark, periods), "benchmark is not finite at 3");
    EXPECT_EQ(data.observation_count(), periods);
    EXPECT_NE(data.returns()[0], 9.0);
}

} // namespace
