#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <portopt/portopt.h>
#include <vector>

#include "capi_scope.hpp"
#include "random_covariance.hpp"

namespace
{

using po::testing::DataPtr;
using po::testing::SpecPtr;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr std::size_t kAssets = 4;

std::vector<double> solved_weights(const po_data* data)
{
    SpecPtr spec = po::testing::make_spec();
    return po::testing::weights_of(po::testing::solve_ok(spec.get(), data).get());
}

void expect_create_fails(std::size_t n, const double* mu, const double* covariance, const char* message)
{
    DataPtr existing = po::testing::make_data(kAssets, {}, {});
    po_data* data = existing.get();
    EXPECT_EQ(po_data_create(n, mu, covariance, &data), PO_INVALID_ARGUMENT);
    EXPECT_EQ(data, nullptr);
    EXPECT_STREQ(po_last_error(), message);
}

TEST(CapiData, CopiesInputs)
{
    std::vector<double> mu = {0.1, 0.2, 0.3, 0.4};
    std::vector<double> covariance = po::testing::random_covariance(kAssets, 3);
    DataPtr data = po::testing::make_data(kAssets, mu, covariance);
    std::vector<double> before = solved_weights(data.get());
    std::fill(mu.begin(), mu.end(), kNaN);
    std::fill(covariance.begin(), covariance.end(), kNaN);
    EXPECT_EQ(solved_weights(data.get()), before);
}

TEST(CapiData, NullCovarianceAllowedAtCreate)
{
    std::vector<double> mu = {0.1, 0.2, 0.3, 0.4};
    DataPtr data = po::testing::make_data(kAssets, mu, {});
    ASSERT_NE(data, nullptr);
    SpecPtr spec = po::testing::make_spec();
    po_result* result = nullptr;
    EXPECT_EQ(po_solve(spec.get(), data.get(), &result), PO_INVALID_ARGUMENT);
    EXPECT_EQ(result, nullptr);
    EXPECT_STREQ(po_last_error(), "MIN_VARIANCE requires covariance");
}

TEST(CapiData, RejectsNaN)
{
    std::vector<double> mu = {0.1, kNaN, 0.3, 0.4};
    std::vector<double> covariance = po::testing::random_covariance(kAssets, 3);
    expect_create_fails(kAssets, mu.data(), nullptr, "mu is not finite at 1");
    mu[1] = 0.2;
    covariance[6] = kNaN;
    expect_create_fails(kAssets, mu.data(), covariance.data(), "covariance is not finite at (1, 2)");
}

TEST(CapiData, RejectsZeroAssets)
{
    expect_create_fails(0, nullptr, nullptr, "asset count must be positive");
}

TEST(CapiData, RejectsNotPositiveDefiniteCovariance)
{
    std::vector<double> covariance(kAssets * kAssets, 1.0);
    expect_create_fails(kAssets, nullptr, covariance.data(), "covariance is not positive definite");
}

TEST(CapiData, RejectsDuplicatedAsset)
{
    std::vector<double> covariance = {0.04, 0.04, 0.01, 0.04, 0.04, 0.01, 0.01, 0.01, 0.09};
    expect_create_fails(3, nullptr, covariance.data(), "covariance is not positive definite");
}

TEST(CapiData, RejectsCountWhoseSquareOverflows)
{
    const std::size_t huge = std::numeric_limits<std::size_t>::max();
    expect_create_fails(huge, nullptr, nullptr, "n is too large");
}

TEST(CapiData, WeightCountWithCash)
{
    DataPtr data = po::testing::make_data(kAssets, {}, {});
    std::size_t assets = 0;
    std::size_t weights = 0;
    ASSERT_EQ(po_data_weight_count(data.get(), &weights), PO_OK);
    EXPECT_EQ(weights, kAssets);
    ASSERT_EQ(po_data_set_risk_free(data.get(), 0.02), PO_OK);
    ASSERT_EQ(po_data_weight_count(data.get(), &weights), PO_OK);
    ASSERT_EQ(po_data_asset_count(data.get(), &assets), PO_OK);
    EXPECT_EQ(weights, kAssets + 1);
    EXPECT_EQ(assets, kAssets);
}

TEST(CapiData, RejectsNonFiniteRiskFreeRate)
{
    DataPtr data = po::testing::make_data(kAssets, {}, {});
    EXPECT_EQ(po_data_set_risk_free(data.get(), kNaN), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_data_set_risk_free(data.get(), std::numeric_limits<double>::infinity()), PO_INVALID_ARGUMENT);
    std::size_t weights = 0;
    ASSERT_EQ(po_data_weight_count(data.get(), &weights), PO_OK);
    EXPECT_EQ(weights, kAssets);
}

} // namespace
