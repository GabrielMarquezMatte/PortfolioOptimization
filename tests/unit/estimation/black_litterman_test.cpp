#include "estimation/black_litterman.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

namespace
{

using po::Data;
using po::estimation::equilibrium;
using po::estimation::omega_he_litterman;
using po::estimation::omega_idzorek;
using po::estimation::Views;

const std::vector<double> kCovariance = {0.04, 0.01, 0.01, 0.09};

Data make_data(std::size_t n, const std::vector<double>& covariance)
{
    auto data = Data::create(n, {}, covariance);
    EXPECT_TRUE(data.has_value());
    return std::move(*data);
}

Views make_views(std::size_t n, std::vector<double> p)
{
    auto views = Views::create(1, n, p, std::vector<double>{0.02});
    EXPECT_TRUE(views.has_value());
    return std::move(*views);
}

template <typename T>
void expect_failure(const std::expected<T, po::Failure>& result, const char* message)
{
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error().kind, po::FailureKind::InvalidArgument);
    EXPECT_EQ(result.error().message, message);
}

TEST(Equilibrium, ByHand)
{
    auto pi = equilibrium(make_data(2, kCovariance), std::vector<double>{0.6, 0.4}, 2.5);
    ASSERT_TRUE(pi.has_value());
    EXPECT_NEAR((*pi)[0], 0.07, 1e-16);
    EXPECT_NEAR((*pi)[1], 0.105, 1e-16);
}

TEST(Equilibrium, FloatingPointSumIsAccepted)
{
    std::vector<double> diagonal = {0.04, 0, 0, 0, 0.09, 0, 0, 0, 0.16};
    double third = 1.0 / 3.0;
    EXPECT_TRUE(equilibrium(make_data(3, diagonal), std::vector<double>{third, third, third}, 2.5).has_value());
}

TEST(Equilibrium, WeightsOffByOnePercentAreRejected)
{
    expect_failure(equilibrium(make_data(2, kCovariance), std::vector<double>{0.6, 0.39}, 2.5),
                   "w_mkt must sum to 1 within 1e-8");
}

TEST(Equilibrium, DeltaMustBePositive)
{
    Data data = make_data(2, kCovariance);
    std::vector<double> w = {0.6, 0.4};
    expect_failure(equilibrium(data, w, 0.0), "delta must be positive and finite");
    expect_failure(equilibrium(data, w, std::nan("")), "delta must be positive and finite");
}

TEST(Equilibrium, InputsAreChecked)
{
    Data data = make_data(2, kCovariance);
    expect_failure(equilibrium(data, std::vector<double>{1.0}, 2.5), "w_mkt has 1 entries; expected 2");
    expect_failure(equilibrium(data, std::vector<double>{std::nan(""), 0.4}, 2.5), "w_mkt is not finite at 0");
    auto no_covariance = Data::create(2, std::vector<double>{0.01, 0.02}, {});
    ASSERT_TRUE(no_covariance.has_value());
    expect_failure(equilibrium(*no_covariance, std::vector<double>{0.6, 0.4}, 2.5), "data has no covariance");
}

TEST(Omega, HeLitterman)
{
    auto omega = omega_he_litterman(make_data(2, kCovariance), make_views(2, {1.0, -1.0}), 0.05);
    ASSERT_TRUE(omega.has_value());
    ASSERT_EQ(omega->size(), 1u);
    EXPECT_NEAR((*omega)[0], 0.0055, 1e-18);
}

TEST(Omega, HeLittermanIsDiagonal)
{
    auto views = Views::create(2, 2, std::vector<double>{1, 0, 0, 1}, std::vector<double>{0.01, 0.02});
    ASSERT_TRUE(views.has_value());
    auto omega = omega_he_litterman(make_data(2, kCovariance), *views, 0.05);
    ASSERT_TRUE(omega.has_value());
    EXPECT_NEAR((*omega)[0], 0.05 * 0.04, 1e-18);
    EXPECT_EQ((*omega)[1], 0.0);
    EXPECT_EQ((*omega)[2], 0.0);
    EXPECT_NEAR((*omega)[3], 0.05 * 0.09, 1e-18);
}

TEST(Omega, Idzorek)
{
    Data data = make_data(2, kCovariance);
    Views views = make_views(2, {1.0, -1.0});
    auto half = omega_idzorek(data, views, 0.05, std::vector<double>{0.5});
    ASSERT_TRUE(half.has_value());
    EXPECT_NEAR((*half)[0], 0.0055, 1e-18);
    auto full = omega_idzorek(data, views, 0.05, std::vector<double>{1.0});
    ASSERT_TRUE(full.has_value());
    EXPECT_EQ((*full)[0], 0.0);
    auto quarter = omega_idzorek(data, views, 0.05, std::vector<double>{0.25});
    ASSERT_TRUE(quarter.has_value());
    EXPECT_NEAR((*quarter)[0], 0.0165, 1e-18);
}

TEST(Omega, IdzorekConfidenceIsChecked)
{
    Data data = make_data(2, kCovariance);
    Views views = make_views(2, {1.0, -1.0});
    expect_failure(omega_idzorek(data, views, 0.05, std::vector<double>{0.0}), "confidence must be in (0, 1] at 0");
    expect_failure(omega_idzorek(data, views, 0.05, std::vector<double>{1.5}), "confidence must be in (0, 1] at 0");
    expect_failure(omega_idzorek(data, views, 0.05, std::vector<double>{std::nan("")}),
                   "confidence must be in (0, 1] at 0");
    expect_failure(omega_idzorek(data, views, 0.05, std::vector<double>{0.5, 0.5}),
                   "confidence has 2 entries; expected 1");
}

TEST(Omega, AssetCountsMustAgree)
{
    expect_failure(omega_he_litterman(make_data(2, kCovariance), make_views(3, {1.0, -1.0, 0.0}), 0.05),
                   "views have 3 assets; data has 2");
}

TEST(Omega, TauMustBePositive)
{
    Data data = make_data(2, kCovariance);
    Views views = make_views(2, {1.0, -1.0});
    expect_failure(omega_he_litterman(data, views, 0.0), "tau must be positive and finite");
    expect_failure(omega_idzorek(data, views, std::nan(""), std::vector<double>{0.5}),
                   "tau must be positive and finite");
}

TEST(Omega, NeedsACovariance)
{
    Views views = make_views(2, {1.0, -1.0});
    auto no_covariance = Data::create(2, std::vector<double>{0.01, 0.02}, {});
    ASSERT_TRUE(no_covariance.has_value());
    expect_failure(omega_he_litterman(*no_covariance, views, 0.05), "data has no covariance");
}

} // namespace
