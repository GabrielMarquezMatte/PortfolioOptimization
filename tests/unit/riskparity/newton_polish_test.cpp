#include "riskparity/newton_polish.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "core/data.hpp"
#include "core/objective.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "fuzz/frontier_reference.hpp"
#include "fuzz/risk_parity/risk_parity_reference.hpp"
#include "random_covariance.hpp"
#include "riskparity/kkt.hpp"
#include "riskparity/polyhedral_descent.hpp"

namespace
{

using po::Data;
using po::Result;
using po::Spec;
using po::riskparity::DescentOptions;
using po::riskparity::DescentProblem;
using po::riskparity::DescentResult;
using po::riskparity::DescentStop;
using po::riskparity::KktPoint;
using po::riskparity::KktProblem;
using po::riskparity::PolishResult;
using po::testing::random_covariance;
using po::testing::random_vector;

constexpr double kBand = 1e-10;
constexpr double kStartTolerance = 1e-6;

Data sample_data(std::size_t n, std::uint64_t seed)
{
    return Data::create(n, random_vector(n, seed + 500, 0.03, 0.15), random_covariance(n, seed, 1e2)).value();
}

// Box [0, 1.2 / n] and an uneven budget, which holds weights at the cap.
Spec capped_spec(std::size_t n)
{
    Spec spec;
    spec.objective = po::Objective::RiskParity;
    spec.lower = 0.0;
    spec.upper = 1.2 / static_cast<double>(n);
    spec.risk_budget = {0.5, 0.2, 0.15, 0.1, 0.05};
    return spec;
}

// A descent stopped at a relative KKT residual of 1e-6: a point near a KKT point, with the active set of its
// last QP in the bound multipliers.
KktPoint near_kkt_point(const Spec& spec, const Data& data)
{
    std::vector<double> start = po::testing::frontier_point(spec, data, spec.target_return).weights;
    DescentProblem problem{spec, data, spec.risk_budget, std::nullopt};
    DescentResult descent =
        po::riskparity::descend(problem, start, DescentOptions{0.05, kStartTolerance, 1000});
    EXPECT_EQ(descent.stop, DescentStop::Converged);
    EXPECT_GT(descent.stationarity, 100.0 * kBand);
    return KktPoint{descent.weights, descent.equality_multipliers, descent.bound_multipliers};
}

Result as_result(const Spec& spec, const KktPoint& point)
{
    Result result;
    result.status = po::SolveStatus::Optimal;
    result.weights = point.weights;
    result.bound_multipliers = point.bound_multipliers;
    result.budget_multiplier = point.equality_multipliers.at(0);
    if (spec.target_return.has_value())
    {
        result.target_return_multiplier = point.equality_multipliers.at(1);
    }
    return result;
}

double merit_of(const KktProblem& problem, const KktPoint& point)
{
    return largest(residual_of(problem, point, evaluate_kkt(problem, point)));
}

// From a KKT point at 1e-6, one or two Newton steps reach the band of cond(Sigma) and the rounding floor of
// the residual; there the next step does not improve it and is refused.
TEST(NewtonPolish, ReachesRoundingFloor)
{
    std::size_t n = 5;
    Data data = sample_data(n, 7);
    for (std::optional<double> target : {std::optional<double>(), std::optional<double>(0.09)})
    {
        SCOPED_TRACE(target.has_value());
        Spec spec = capped_spec(n);
        spec.target_return = target;
        KktProblem problem{spec, data, spec.risk_budget};
        KktPoint start = near_kkt_point(spec, data);
        double start_merit = merit_of(problem, start);
        PolishResult polished = po::riskparity::newton_polish(problem, start);
        ASSERT_TRUE(polished.improved);
        EXPECT_GE(polished.steps, 1);
        EXPECT_LE(polished.steps, 2);
        EXPECT_EQ(polished.residual, merit_of(problem, polished.point));
        EXPECT_LT(polished.residual, 1e-6 * start_merit);
        EXPECT_LE(polished.residual, 1.0);
        po::testing::expect_risk_parity_kkt(spec, data, as_result(spec, polished.point), kBand);
        bool capped = false;
        for (std::size_t i = 0; i < n; ++i)
        {
            capped = capped || polished.point.bound_multipliers[i] < 0.0;
        }
        EXPECT_TRUE(capped);
        PolishResult again = po::riskparity::newton_polish(problem, polished.point);
        EXPECT_FALSE(again.improved);
        EXPECT_EQ(again.steps, 0);
        EXPECT_EQ(again.point.weights, polished.point.weights);
        EXPECT_EQ(again.point.equality_multipliers, polished.point.equality_multipliers);
        EXPECT_EQ(again.point.bound_multipliers, polished.point.bound_multipliers);
        EXPECT_EQ(again.residual, polished.residual);
    }
}

// A start that is not near any KKT point of its active set: the weights of the start with multipliers of the
// wrong size. Whatever the steps do, the result is never worse than the start.
TEST(NewtonPolish, NeverReturnsAWorsePoint)
{
    std::size_t n = 5;
    Data data = sample_data(n, 7);
    Spec spec = capped_spec(n);
    KktProblem problem{spec, data, spec.risk_budget};
    for (double factor : {-50.0, 0.0, 3.0, 1e6})
    {
        SCOPED_TRACE(factor);
        KktPoint start = near_kkt_point(spec, data);
        for (double& multiplier : start.equality_multipliers)
        {
            multiplier *= factor;
        }
        double start_merit = merit_of(problem, start);
        PolishResult polished = po::riskparity::newton_polish(problem, start);
        EXPECT_LE(polished.residual, start_merit);
        EXPECT_EQ(polished.improved, polished.residual < start_merit);
        if (!polished.improved)
        {
            EXPECT_EQ(polished.point.weights, start.weights);
            EXPECT_EQ(polished.steps, 0);
        }
    }
}

// The KKT residual in units of its noise separates a KKT point from a point an error of 1e-9 of a weight
// away from it.
TEST(Kkt, ResidualSeesAnErrorOfTheWeights)
{
    std::size_t n = 5;
    Data data = sample_data(n, 7);
    Spec spec = capped_spec(n);
    spec.target_return = 0.09;
    KktProblem problem{spec, data, spec.risk_budget};
    PolishResult polished = po::riskparity::newton_polish(problem, near_kkt_point(spec, data));
    ASSERT_TRUE(polished.improved);
    KktPoint moved = polished.point;
    std::size_t first_free = n;
    std::size_t second_free = n;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (moved.bound_multipliers[i] != 0.0)
        {
            continue;
        }
        second_free = first_free;
        first_free = i;
    }
    ASSERT_LT(second_free, n);
    moved.weights[first_free] += 1e-9;
    moved.weights[second_free] -= 1e-9;
    po::riskparity::KktResidual residual = residual_of(problem, moved, evaluate_kkt(problem, moved));
    EXPECT_GT(residual.stationarity, 1e3);
    EXPECT_GT(residual.equalities, 1e3);
    EXPECT_EQ(residual.bounds, 0.0);
    KktPoint outside = polished.point;
    outside.weights[first_free] = spec.upper + 1e-9;
    EXPECT_GT(residual_of(problem, outside, evaluate_kkt(problem, outside)).bounds, 1e3);
    KktPoint wrong_sign = polished.point;
    for (double& beta : wrong_sign.bound_multipliers)
    {
        beta = -beta;
    }
    EXPECT_GT(residual_of(problem, wrong_sign, evaluate_kkt(problem, wrong_sign)).complementarity, 1e3);
    KktPoint not_a_number = polished.point;
    not_a_number.weights[0] = std::nan("");
    EXPECT_TRUE(std::isnan(largest(residual_of(problem, not_a_number, evaluate_kkt(problem, not_a_number)))));
}

// Two weights pinned by their bounds at (1 + t, -t), with beta = grad f, so that the stationarity and the
// bounds hold exactly. At t = 1 the point is a KKT point. At t = 1e16 the first weight rounds to t and the
// weights sum to 0, not to 1: a gap of 1, which is within the rounding of weights of that size, so the
// residual alone cannot reject the point. The budget row that is no longer resolved does.
TEST(Kkt, WeightsThatDoNotResolveTheBudgetAreNotAKktPoint)
{
    Data data = Data::create(2, {}, std::vector<double>{0.04, 0.006, 0.006, 0.09}).value();
    for (double size : {1.0, 1e16})
    {
        SCOPED_TRACE(size);
        Spec spec;
        spec.objective = po::Objective::RiskParity;
        spec.risk_budget = {0.5, 0.5};
        spec.asset_lower = {1.0 + size, -size};
        spec.asset_upper = spec.asset_lower;
        KktProblem problem{spec, data, spec.risk_budget};
        KktPoint point{spec.asset_lower, {0.0}, {0.0, 0.0}};
        point.bound_multipliers = evaluate_kkt(problem, point).linearization.gradient;
        double residual = largest(residual_of(problem, point, evaluate_kkt(problem, point)));
        EXPECT_LE(residual, 1.0);
        EXPECT_EQ(is_kkt_point(problem, point), size == 1.0);
    }
    EXPECT_EQ((1.0 + 1e16) + (-1e16), 0.0);
}

} // namespace
