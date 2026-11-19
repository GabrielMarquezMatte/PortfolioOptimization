#include "problems/parametric_path.hpp"

#include <gtest/gtest.h>

#include <cstddef>
#include <optional>
#include <vector>

#include "core/data.hpp"
#include "core/spec.hpp"

namespace
{

using po::Data;
using po::Spec;
using po::problems::Crossing;
using po::problems::Path;
using po::problems::PathPoint;

// One generated row active, x = (0.6, 0, -0.4): the support {0: +, 2: -} and the zero at 1, with signs, the sign
// vector of the row.
PathPoint leverage_point(std::vector<double> x, std::vector<double> signs)
{
    PathPoint point;
    point.t = 1.0;
    point.solution.status = po::qp::Status::Optimal;
    point.solution.x = std::move(x);
    point.solution.gross_exposure_signs = std::move(signs);
    point.solution.gross_exposure_multipliers = {1.0};
    return point;
}

TEST(ParametricPath, SameSegmentIgnoresTheSignVectorRepresentation)
{
    PathPoint earlier = leverage_point({0.6, 0.0, -0.4}, {1.0, 1.0, -1.0});
    PathPoint later = leverage_point({0.6, 0.0, -0.4}, {1.0, -1.0, -1.0});
    EXPECT_TRUE(po::problems::same_segment(earlier, later));
    later.solution.x[2] = 0.4;
    EXPECT_FALSE(po::problems::same_segment(earlier, later));
}

// Three assets, long-only with the budget, a well-conditioned Sigma. From t = 0 to about 1.33 no bound is active
// (vol 0.1642 to 0.2745); there w_0 reaches 0 and the next segment begins.
Data three_assets()
{
    std::vector<double> mu = {0.05, 0.08, 0.12};
    std::vector<double> covariance = {0.04, 0.006, 0.004, 0.006, 0.09, 0.012, 0.004, 0.012, 0.16};
    return Data::create(3, mu, covariance).value();
}

Spec long_only_max_return(double cap)
{
    Spec spec;
    spec.objective = po::Objective::MaxReturn;
    spec.lower = 0.0;
    spec.upper = 1.0;
    spec.max_volatility = cap;
    return spec;
}

// The cap 0.22 is reached inside the first segment; the weights are those po::solve gave before the tangent
// crossing (commit 7287d84), which interpolated two solutions of that segment.
TEST(ParametricPath, TangentCrossingMatchesTheSolveWithoutLeverage)
{
    Data data = three_assets();
    Path path = po::problems::make_path(long_only_max_return(0.22), data);
    PathPoint point = po::problems::solve_at(path, 0.3);
    ASSERT_EQ(point.solution.status, po::qp::Status::Optimal);
    ASSERT_TRUE(point.solution.active_inequalities.empty());
    std::optional<Crossing> crossing = po::problems::tangent_crossing(path, point);
    ASSERT_TRUE(crossing.has_value());
    EXPECT_GT(crossing->t, point.t);
    std::vector<double> before = {0x1.ae46ac2aff0ap-3, 0x1.6fd7173f2bac8p-2, 0x1.b90592ab54ce8p-2};
    for (std::size_t i = 0; i < before.size(); ++i)
    {
        SCOPED_TRACE(i);
        EXPECT_NEAR(crossing->solution.x[i], before[i], 1e-12);
    }
}

// The cap 0.29 is beyond the end of the first segment (vol 0.2745 where w_0 reaches 0): the tangent of the point
// at t = 0.3 meets it with w_0 < 0, off the face, and the certificate refuses it.
TEST(ParametricPath, TangentCrossingRejectsAPointOffItsFace)
{
    Data data = three_assets();
    Path path = po::problems::make_path(long_only_max_return(0.29), data);
    PathPoint point = po::problems::solve_at(path, 0.3);
    ASSERT_EQ(point.solution.status, po::qp::Status::Optimal);
    EXPECT_FALSE(po::problems::tangent_crossing(path, point).has_value());
}

} // namespace
