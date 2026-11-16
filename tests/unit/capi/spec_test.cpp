#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <portopt/portopt.h>
#include <vector>

#include "capi_scope.hpp"
#include "random_covariance.hpp"

namespace
{

using po::testing::DataPtr;
using po::testing::SpecPtr;

constexpr double kInfinity = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr std::size_t kAssets = 3;

DataPtr three_asset_data()
{
    return po::testing::make_data(kAssets, {}, po::testing::random_covariance(kAssets, 7));
}

DataPtr three_asset_data_with_mu()
{
    return po::testing::make_data(kAssets, po::testing::random_vector(kAssets, 3, 0.03, 0.15),
                                  po::testing::random_covariance(kAssets, 7));
}

TEST(CapiSpec, RejectsUnknownObjective)
{
    for (int32_t objective : {-1, 9})
    {
        SCOPED_TRACE(objective);
        SpecPtr existing = po::testing::make_spec();
        po_spec* spec = existing.get();
        EXPECT_EQ(po_spec_create(objective, &spec), PO_INVALID_ARGUMENT);
        EXPECT_EQ(spec, nullptr);
        EXPECT_STREQ(po_last_error(), "objective is not one of PO_OBJECTIVE_*");
    }
}

TEST(CapiSpec, AcceptsEveryObjective)
{
    for (int32_t objective = PO_OBJECTIVE_MIN_VARIANCE; objective <= PO_OBJECTIVE_RISK_PARITY; ++objective)
    {
        SCOPED_TRACE(objective);
        SpecPtr spec = po::testing::make_spec(objective);
        EXPECT_NE(spec, nullptr);
        EXPECT_STREQ(po_last_error(), "");
    }
}

TEST(CapiSpec, CreateRejectsNullOut)
{
    EXPECT_EQ(po_spec_create(PO_OBJECTIVE_MIN_VARIANCE, nullptr), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "out must not be null");
}

TEST(CapiSpec, InvalidBoundsLeaveSpecUnchanged)
{
    DataPtr data = three_asset_data();
    SpecPtr spec = po::testing::make_spec();
    ASSERT_EQ(po_spec_set_bounds(spec.get(), 0.0, 1.0), PO_OK);
    std::vector<double> expected = po::testing::weights_of(po::testing::solve_ok(spec.get(), data.get()).get());

    EXPECT_EQ(po_spec_set_bounds(spec.get(), kNaN, 1.0), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_spec_set_bounds(spec.get(), 0.0, kNaN), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_spec_set_bounds(spec.get(), 0.5, 0.2), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_spec_set_bounds(spec.get(), kInfinity, kInfinity), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_spec_set_bounds(spec.get(), -kInfinity, -kInfinity), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(),
                 "bounds must not be NaN and must satisfy lower <= upper, lower < +inf and upper > -inf");

    const double lower[kAssets] = {0.0, 0.0, 0.0};
    const double upper[kAssets] = {1.0, kNaN, 1.0};
    EXPECT_EQ(po_spec_set_asset_bounds(spec.get(), lower, upper, kAssets), PO_INVALID_ARGUMENT);
    const double crossed[kAssets] = {1.0, 1.0, -1.0};
    EXPECT_EQ(po_spec_set_asset_bounds(spec.get(), crossed, lower, kAssets), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_spec_set_target_return(spec.get(), kNaN), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_spec_set_target_return(spec.get(), kInfinity), PO_INVALID_ARGUMENT);

    EXPECT_EQ(po::testing::weights_of(po::testing::solve_ok(spec.get(), data.get()).get()), expected);
}

TEST(CapiSpec, AssetBoundsRejectNullArraysAndZeroCount)
{
    SpecPtr spec = po::testing::make_spec();
    const double values[kAssets] = {0.0, 0.0, 0.0};
    EXPECT_EQ(po_spec_set_asset_bounds(spec.get(), nullptr, values, kAssets), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_spec_set_asset_bounds(spec.get(), values, nullptr, kAssets), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_spec_set_asset_bounds(spec.get(), values, values, 0), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "n must be positive");
}

TEST(CapiSpec, LastBoundsCallWins)
{
    DataPtr data = three_asset_data();
    SpecPtr spec = po::testing::make_spec();
    const double lower[kAssets] = {-kInfinity, -kInfinity, -kInfinity};
    const double upper[kAssets] = {0.2, 0.2, 0.2};

    ASSERT_EQ(po_spec_set_bounds(spec.get(), 0.0, 1.0), PO_OK);
    ASSERT_EQ(po_spec_set_asset_bounds(spec.get(), lower, upper, kAssets), PO_OK);
    EXPECT_EQ(po::testing::status_of(po::testing::solve_ok(spec.get(), data.get()).get()), PO_SOLVE_INFEASIBLE);

    ASSERT_EQ(po_spec_set_bounds(spec.get(), 0.0, 1.0), PO_OK);
    EXPECT_EQ(po::testing::status_of(po::testing::solve_ok(spec.get(), data.get()).get()), PO_SOLVE_OPTIMAL);

    ASSERT_EQ(po_spec_set_bounds(spec.get(), 0.0, 0.2), PO_OK);
    EXPECT_EQ(po::testing::status_of(po::testing::solve_ok(spec.get(), data.get()).get()), PO_SOLVE_INFEASIBLE);

    const double wide[kAssets] = {1.0, 1.0, 1.0};
    ASSERT_EQ(po_spec_set_asset_bounds(spec.get(), lower, wide, kAssets), PO_OK);
    EXPECT_EQ(po::testing::status_of(po::testing::solve_ok(spec.get(), data.get()).get()), PO_SOLVE_OPTIMAL);
}

TEST(CapiSpec, WrongAssetBoundCountFailsAtSolve)
{
    DataPtr data = three_asset_data();
    SpecPtr spec = po::testing::make_spec();
    const double lower[2] = {0.0, 0.0};
    const double upper[2] = {1.0, 1.0};
    ASSERT_EQ(po_spec_set_asset_bounds(spec.get(), lower, upper, 2), PO_OK);
    po_result* result = nullptr;
    EXPECT_EQ(po_solve(spec.get(), data.get(), &result), PO_INVALID_ARGUMENT);
    EXPECT_EQ(result, nullptr);
    EXPECT_STREQ(po_last_error(), "asset_lower has 2 entries; expected 3");
}

using Setter = po_status (*)(po_spec*, double);

void expect_domain(int32_t objective, Setter setter, double valid, const std::vector<double>& invalid,
                   const char* message)
{
    DataPtr data = three_asset_data_with_mu();
    SpecPtr spec = po::testing::make_spec(objective);
    ASSERT_EQ(po_spec_set_bounds(spec.get(), 0.0, 1.0), PO_OK);
    ASSERT_EQ(setter(spec.get(), valid), PO_OK);
    std::vector<double> expected = po::testing::weights_of(po::testing::solve_ok(spec.get(), data.get()).get());
    for (double value : invalid)
    {
        SCOPED_TRACE(value);
        EXPECT_EQ(setter(spec.get(), value), PO_INVALID_ARGUMENT);
        EXPECT_STREQ(po_last_error(), message);
        EXPECT_EQ(po::testing::weights_of(po::testing::solve_ok(spec.get(), data.get()).get()), expected);
    }
}

TEST(CapiSpec, RiskAversionDomain)
{
    expect_domain(PO_OBJECTIVE_UTILITY, po_spec_set_risk_aversion, 2.0, {0.0, -1.0, kNaN, kInfinity},
                  "delta must be positive and finite");
}

TEST(CapiSpec, ImportanceDomain)
{
    expect_domain(PO_OBJECTIVE_IMPORTANCE, po_spec_set_importance, 0.3, {1.0, -0.1, kNaN, 2.0},
                  "lambda must satisfy 0 <= lambda < 1");
    SpecPtr spec = po::testing::make_spec(PO_OBJECTIVE_IMPORTANCE);
    EXPECT_EQ(po_spec_set_importance(spec.get(), 0.0), PO_OK);
}

TEST(CapiSpec, MaxVolatilityDomain)
{
    expect_domain(PO_OBJECTIVE_MAX_RETURN, po_spec_set_max_volatility, 5.0, {0.0, -0.1, kNaN, kInfinity},
                  "sigma must be positive and finite");
}

TEST(CapiSpec, LeverageDomain)
{
    expect_domain(PO_OBJECTIVE_MIN_VARIANCE, po_spec_set_leverage, 2.0, {0.0, -1.0, kNaN, kInfinity, -kInfinity},
                  "max_gross must be positive and finite");
    EXPECT_EQ(po_spec_set_leverage(nullptr, 2.0), PO_INVALID_HANDLE);
}

} // namespace
