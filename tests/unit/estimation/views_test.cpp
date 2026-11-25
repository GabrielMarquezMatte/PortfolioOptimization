#include "estimation/views.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

namespace
{

using po::estimation::Views;

void expect_failure(const std::expected<Views, po::Failure>& result, const char* message)
{
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error().kind, po::FailureKind::InvalidArgument);
    EXPECT_EQ(result.error().message, message);
}

Views two_views()
{
    auto views = Views::create(2, 2, std::vector<double>{1, 0, 0, 1}, std::vector<double>{0.01, 0.02});
    EXPECT_TRUE(views.has_value());
    return std::move(*views);
}

TEST(Views, KeepsWhatItWasGiven)
{
    auto views = Views::create(1, 2, std::vector<double>{1.0, -1.0}, std::vector<double>{0.02});
    ASSERT_TRUE(views.has_value());
    EXPECT_EQ(views->view_count(), 1u);
    EXPECT_EQ(views->asset_count(), 2u);
    EXPECT_EQ(views->p()[1], -1.0);
    EXPECT_EQ(views->q()[0], 0.02);
    EXPECT_FALSE(views->has_omega());
}

TEST(Views, RedundantViewsAreRejected)
{
    expect_failure(Views::create(2, 2, std::vector<double>{1, -1, 2, -2}, std::vector<double>{0.01, 0.02}),
                   "views are redundant: P does not have full row rank");
}

TEST(Views, MoreViewsThanAssetsAreRejected)
{
    expect_failure(Views::create(3, 2, std::vector<double>(6, 1.0), std::vector<double>(3, 0.0)),
                   "k must not exceed n");
}

TEST(Views, EmptyViewsAreAccepted)
{
    auto views = Views::create(0, 2, {}, {});
    ASSERT_TRUE(views.has_value());
    EXPECT_EQ(views->view_count(), 0u);
    EXPECT_TRUE(views->p().empty());
}

TEST(Views, ShapeAndFinitenessAreChecked)
{
    expect_failure(Views::create(1, 0, {}, std::vector<double>{0.0}), "asset count must be positive");
    expect_failure(Views::create(1, 2, std::vector<double>{1.0}, std::vector<double>{0.0}),
                   "p has 1 entries; expected 2");
    expect_failure(Views::create(1, 2, std::vector<double>{1.0, -1.0}, std::vector<double>{}),
                   "q has 0 entries; expected 1");
    expect_failure(Views::create(1, 2, std::vector<double>{1.0, std::nan("")}, std::vector<double>{0.0}),
                   "p is not finite at (0, 1)");
    expect_failure(Views::create(1, 2, std::vector<double>{1.0, -1.0}, std::vector<double>{HUGE_VAL}),
                   "q is not finite at 0");
}

TEST(Views, OmegaIsCopiedAndKept)
{
    Views views = two_views();
    ASSERT_TRUE(views.set_omega(std::vector<double>{0.1, 0.01, 0.01, 0.2}).has_value());
    EXPECT_TRUE(views.has_omega());
    EXPECT_EQ(views.omega()[3], 0.2);
}

TEST(Views, OmegaMustBeSymmetricWithNonNegativeDiagonal)
{
    Views views = two_views();
    ASSERT_TRUE(views.set_omega(std::vector<double>{0.1, 0.0, 0.0, 0.2}).has_value());
    auto asymmetric = views.set_omega(std::vector<double>{0.1, 0.01, 0.02, 0.2});
    ASSERT_FALSE(asymmetric.has_value());
    EXPECT_EQ(asymmetric.error().message, "omega is not symmetric at (0, 1)");
    auto negative = views.set_omega(std::vector<double>{0.1, 0.0, 0.0, -0.2});
    ASSERT_FALSE(negative.has_value());
    EXPECT_EQ(negative.error().message, "omega has a negative diagonal entry at 1");
    EXPECT_EQ(views.omega()[0], 0.1);
    EXPECT_EQ(views.omega()[3], 0.2);
}

TEST(Views, OmegaSizeAndFinitenessAreChecked)
{
    Views views = two_views();
    auto size = views.set_omega(std::vector<double>{1.0, 2.0});
    ASSERT_FALSE(size.has_value());
    EXPECT_EQ(size.error().message, "omega has 2 entries; expected 4");
    auto nan = views.set_omega(std::vector<double>{0.1, 0.0, 0.0, std::nan("")});
    ASSERT_FALSE(nan.has_value());
    EXPECT_EQ(nan.error().message, "omega is not finite at (1, 1)");
    EXPECT_FALSE(views.has_omega());
}

} // namespace
