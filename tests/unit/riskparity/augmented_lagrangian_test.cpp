#include "riskparity/augmented_lagrangian.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "core/data.hpp"
#include "core/objective.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "fuzz/frontier_reference.hpp"
#include "fuzz/risk_parity/risk_parity_reference.hpp"
#include "random_covariance.hpp"
#include "riskparity/kkt.hpp"

namespace
{

using po::Data;
using po::Result;
using po::Spec;
using po::riskparity::AugmentedOptions;
using po::riskparity::AugmentedResult;
using po::riskparity::KktProblem;
using po::testing::random_covariance;
using po::testing::random_vector;

constexpr double kBand = 1e-10;
constexpr AugmentedOptions kOptions{1e-11, 60, 2000};

Data sample_data(std::size_t n, std::uint64_t seed)
{
    return Data::create(n, random_vector(n, seed + 500, 0.03, 0.15), random_covariance(n, seed, 1e2)).value();
}

Spec long_only_spec()
{
    Spec spec;
    spec.objective = po::Objective::RiskParity;
    spec.lower = 0.0;
    spec.upper = 1.0;
    spec.risk_budget = {0.3, 0.25, 0.2, 0.15, 0.1};
    return spec;
}

// The minimum-variance point of the polyhedron of spec and its volatility.
struct Minimum
{
    std::vector<double> weights;
    double volatility;
};

Minimum minimum_of(const Spec& spec, const Data& data)
{
    Result minimum = po::testing::frontier_point(spec, data, spec.target_return);
    EXPECT_EQ(minimum.status, po::SolveStatus::Optimal);
    return Minimum{minimum.weights, po::testing::portfolio_volatility(data, minimum.weights)};
}

// The equality multipliers of the result in the order budget, target_return, volatility.
Result as_result(const Spec& spec, const AugmentedResult& augmented)
{
    Result result;
    result.status = po::SolveStatus::Optimal;
    result.weights = augmented.point.weights;
    result.bound_multipliers = augmented.point.bound_multipliers;
    result.budget_multiplier = augmented.point.equality_multipliers.at(0);
    if (spec.target_return.has_value())
    {
        result.target_return_multiplier = augmented.point.equality_multipliers.at(1);
    }
    result.volatility_multiplier = augmented.point.equality_multipliers.back();
    return result;
}

TEST(AugmentedLagrangian, ReachesTheTargetVolatilityAtAKktPoint)
{
    Data data = sample_data(5, 23);
    for (bool with_target_return : {false, true})
    {
        SCOPED_TRACE(with_target_return);
        Spec spec = long_only_spec();
        if (with_target_return)
        {
            spec.target_return = 0.09;
        }
        Minimum minimum = minimum_of(spec, data);
        spec.target_volatility = 1.2 * minimum.volatility;
        KktProblem problem{spec, data, spec.risk_budget};
        AugmentedResult result = po::riskparity::augmented_lagrangian(problem, minimum.weights, kOptions);
        ASSERT_TRUE(result.converged) << "outer " << result.outer_iterations << " penalty " << result.penalty;
        EXPECT_LT(result.outer_iterations, kOptions.max_outer_iterations);
        EXPECT_GT(result.inner_iterations, 0);
        ASSERT_EQ(result.point.equality_multipliers.size(), with_target_return ? 3u : 2u);
        po::testing::expect_risk_parity_kkt(spec, data, as_result(spec, result), kBand);
    }
}

// Close to the minimum volatility the multiplier grows like 1 / sqrt(s / sigma_min - 1) and the multiplier
// method needs a large penalty. Both stay finite and inside their bounds whether or not the iteration ends
// converged, and the iterate is at the target to the size the penalty enforces.
TEST(AugmentedLagrangian, SafeguardedMultiplierStaysBounded)
{
    Data data = sample_data(5, 23);
    Spec spec = long_only_spec();
    Minimum minimum = minimum_of(spec, data);
    double target = minimum.volatility * (1.0 + 1e-6);
    spec.target_volatility = target;
    KktProblem problem{spec, data, spec.risk_budget};
    AugmentedResult result = po::riskparity::augmented_lagrangian(problem, minimum.weights, kOptions);
    ASSERT_EQ(result.point.equality_multipliers.size(), 2u);
    double multiplier = result.point.equality_multipliers[1];
    ASSERT_TRUE(std::isfinite(multiplier));
    EXPECT_LT(multiplier * target, -1.0);
    double gap = po::testing::portfolio_volatility(data, result.point.weights) / target - 1.0;
    EXPECT_LE(std::fabs(multiplier * target),
              po::riskparity::kAugmentedMultiplierBound + result.penalty * std::fabs(gap));
    EXPECT_TRUE(std::isfinite(result.penalty));
    EXPECT_GT(result.penalty, 0.0);
    EXPECT_LE(result.penalty, po::riskparity::kAugmentedPenaltyCeiling);
    EXPECT_LE(result.outer_iterations, kOptions.max_outer_iterations);
    EXPECT_NEAR(po::testing::portfolio_volatility(data, result.point.weights), target, 1e-6 * target);
    for (double weight : result.point.weights)
    {
        EXPECT_TRUE(std::isfinite(weight));
    }
}

TEST(AugmentedLagrangian, IterationCapsAreReported)
{
    Data data = sample_data(5, 23);
    Spec spec = long_only_spec();
    Minimum minimum = minimum_of(spec, data);
    spec.target_volatility = 1.2 * minimum.volatility;
    KktProblem problem{spec, data, spec.risk_budget};
    AugmentedResult outer_capped =
        po::riskparity::augmented_lagrangian(problem, minimum.weights, AugmentedOptions{1e-11, 2, 2000});
    EXPECT_FALSE(outer_capped.converged);
    EXPECT_EQ(outer_capped.outer_iterations, 2);
    AugmentedResult inner_capped =
        po::riskparity::augmented_lagrangian(problem, minimum.weights, AugmentedOptions{1e-11, 3, 1});
    EXPECT_FALSE(inner_capped.converged);
    EXPECT_EQ(inner_capped.outer_iterations, 3);
    EXPECT_EQ(inner_capped.inner_iterations, 3);
}

// Sigma = diag(0.04, 0.09), long only: the largest volatility of the simplex is 0.3, so s = 0.5 is out of
// reach. The iteration parks at a point with vol < s and reports that it did not converge.
TEST(AugmentedLagrangian, UnreachableTargetDoesNotConverge)
{
    Data data = Data::create(2, {}, std::vector<double>{0.04, 0.0, 0.0, 0.09}).value();
    Spec spec;
    spec.objective = po::Objective::RiskParity;
    spec.lower = 0.0;
    spec.risk_budget = {0.5, 0.5};
    spec.target_volatility = 0.5;
    KktProblem problem{spec, data, spec.risk_budget};
    std::vector<double> start = {0.6, 0.4};
    AugmentedResult result = po::riskparity::augmented_lagrangian(problem, start, kOptions);
    EXPECT_FALSE(result.converged);
    EXPECT_GT(result.outer_iterations, 1);
    EXPECT_LE(result.outer_iterations, kOptions.max_outer_iterations);
    EXPECT_GT(result.penalty, 0.0);
    EXPECT_LE(result.penalty, po::riskparity::kAugmentedPenaltyCeiling);
    EXPECT_NEAR(po::testing::portfolio_volatility(data, result.point.weights), 0.3, 1e-9);
    EXPECT_NEAR(result.point.weights[1], 1.0, 1e-9);
}

} // namespace
