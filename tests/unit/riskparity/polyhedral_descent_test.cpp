#include "riskparity/polyhedral_descent.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <vector>

#include "core/data.hpp"
#include "core/objective.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "fuzz/frontier_reference.hpp"
#include "fuzz/risk_parity/risk_parity_reference.hpp"
#include "random_covariance.hpp"
#include "riskparity/spinu.hpp"

namespace
{

using po::Data;
using po::Result;
using po::Spec;
using po::riskparity::DescentOptions;
using po::riskparity::DescentProblem;
using po::riskparity::DescentResult;
using po::riskparity::DescentStop;
using po::riskparity::VolatilityTerm;
using po::testing::random_covariance;
using po::testing::random_vector;

constexpr double kProximal = 0.05;
constexpr double kTolerance = 1e-12;
constexpr double kBand = 1e-10;

Data sample_data(std::size_t n, std::uint64_t seed)
{
    return Data::create(n, random_vector(n, seed + 500, 0.03, 0.15), random_covariance(n, seed, 1e2)).value();
}

Spec risk_parity_spec(double lower, double upper)
{
    Spec spec;
    spec.objective = po::Objective::RiskParity;
    spec.lower = lower;
    spec.upper = upper;
    return spec;
}

DescentOptions options_with(int max_iterations)
{
    return DescentOptions{kProximal, kTolerance, max_iterations};
}

// The descent read as a Result: equality multipliers in the order budget, target_return.
Result as_result(const Spec& spec, const DescentResult& descent)
{
    Result result;
    result.status = po::SolveStatus::MaxIterations;
    if (descent.stop == DescentStop::Converged)
    {
        result.status = po::SolveStatus::Optimal;
    }
    result.weights = descent.weights;
    result.bound_multipliers = descent.bound_multipliers;
    result.budget_multiplier = descent.equality_multipliers.at(0);
    if (spec.target_return.has_value())
    {
        result.target_return_multiplier = descent.equality_multipliers.at(1);
    }
    result.iterations = descent.iterations;
    return result;
}

double volatility_term(const Data& data, const VolatilityTerm& term, const std::vector<double>& weights)
{
    double gap = po::testing::portfolio_volatility(data, weights) - term.target;
    return -term.multiplier * gap + 0.5 * term.penalty * gap * gap;
}

double merit(const DescentProblem& problem, const std::vector<double>& weights)
{
    double value = po::testing::risk_parity_value(problem.data, problem.budget, weights);
    if (problem.volatility.has_value())
    {
        value += volatility_term(problem.data, *problem.volatility, weights);
    }
    return value;
}

// eps sum_k |u_k| |w|_1: how far L moves when the rounding of the weights, of order eps |w_j| each, crosses
// the rows of the polyhedron, whose multipliers make grad L large across them. Every accepted step decreases
// L in exact arithmetic; a computed L may still grow by this much once the decrease is below it.
double row_rounding(const DescentResult& descent)
{
    double multipliers = 0.0;
    for (double multiplier : descent.equality_multipliers)
    {
        multipliers += std::fabs(multiplier);
    }
    for (double multiplier : descent.bound_multipliers)
    {
        multipliers += std::fabs(multiplier);
    }
    double gross = 0.0;
    for (double weight : descent.weights)
    {
        gross += std::fabs(weight);
    }
    return std::numeric_limits<double>::epsilon() * multipliers * gross;
}

TEST(PolyhedralDescent, AtSpinuPointStaysPut)
{
    std::size_t n = 6;
    Data data = sample_data(n, 3);
    std::vector<double> budget = {0.3, 0.2, 0.1, 0.1, 0.2, 0.1};
    po::riskparity::SpinuResult spinu = po::riskparity::spinu(data, budget);
    ASSERT_TRUE(spinu.converged);
    Spec spec = risk_parity_spec(-1.0, 2.0);
    DescentProblem problem{spec, data, budget, std::nullopt};
    DescentResult result = po::riskparity::descend(problem, spinu.weights, options_with(100));
    EXPECT_EQ(result.stop, DescentStop::Converged);
    EXPECT_EQ(result.iterations, 0);
    EXPECT_EQ(result.weights, spinu.weights);
    EXPECT_EQ(result.bound_multipliers, std::vector<double>(n, 0.0));
    ASSERT_EQ(result.equality_multipliers.size(), 1u);
    EXPECT_LE(std::fabs(result.equality_multipliers[0]),
              kTolerance * po::testing::unit_gap_gradient_scale(data, result.weights));
    EXPECT_LE(result.stationarity, kTolerance);
}

TEST(PolyhedralDescent, TargetReturnLineAgainstScan)
{
    using Shape = po::testing::RiskParityShape;
    for (Shape::Box box : {Shape::Box::LongOnly, Shape::Box::Capped})
    {
        for (Shape::Budget budget : {Shape::Budget::Equal, Shape::Budget::Uneven})
        {
            for (std::uint64_t seed = 1; seed <= 10; ++seed)
            {
                SCOPED_TRACE(testing::Message() << "box " << static_cast<int>(box) << " budget "
                                                << static_cast<int>(budget) << " seed " << seed);
                po::testing::GeneratedRiskParity generated =
                    po::testing::generate(Shape{3, 1e4, 1.0, budget, box, true, false, seed});
                const Spec& spec = generated.spec;
                const Data& data = generated.data;
                po::testing::ScanCurve curve = po::testing::return_line_curve(spec, data);
                std::vector<double> start = curve.point(0.5 * (curve.first + curve.last));
                DescentResult result = po::riskparity::descend(
                    DescentProblem{spec, data, spec.risk_budget, std::nullopt}, start, options_with(1000));
                ASSERT_EQ(result.stop, DescentStop::Converged) << "stationarity " << result.stationarity;
                po::testing::expect_risk_parity_kkt(spec, data, as_result(spec, result), kBand);
                po::testing::ScanResult scan = po::testing::scan_reference(spec, data, spec.risk_budget);
                double value = po::testing::risk_parity_value(data, spec.risk_budget, result.weights);
                if (value > scan.value * (1.0 + 1e-9))
                {
                    double distance = 0.0;
                    for (std::size_t i = 0; i < 3; ++i)
                    {
                        distance = std::fmax(distance, std::fabs(result.weights[i] - scan.weights[i]));
                    }
                    EXPECT_GT(distance, 1e-6) << "value " << value << " scan " << scan.value;
                }
            }
        }
    }
}

TEST(PolyhedralDescent, ActiveBoundsKkt)
{
    std::size_t n = 5;
    Data data = sample_data(n, 7);
    Spec spec = risk_parity_spec(0.0, 1.2 / static_cast<double>(n));
    spec.risk_budget = {0.5, 0.2, 0.15, 0.1, 0.05};
    std::vector<double> start(n, 1.0 / static_cast<double>(n));
    DescentResult result = po::riskparity::descend(DescentProblem{spec, data, spec.risk_budget, std::nullopt},
                                                   start, options_with(1000));
    ASSERT_EQ(result.stop, DescentStop::Converged) << "stationarity " << result.stationarity;
    bool capped = false;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (result.bound_multipliers[i] < 0.0)
        {
            capped = true;
            EXPECT_NEAR(result.weights[i], spec.upper, kBand) << "weight " << i;
        }
        if (result.bound_multipliers[i] > 0.0)
        {
            EXPECT_NEAR(result.weights[i], spec.lower, kBand) << "weight " << i;
        }
    }
    EXPECT_TRUE(capped);
    po::testing::expect_risk_parity_kkt(spec, data, as_result(spec, result), kBand);
}

TEST(PolyhedralDescent, MonotoneDescent)
{
    std::size_t n = 8;
    Data data = sample_data(n, 19);
    Spec spec = risk_parity_spec(0.0, 1.2 / static_cast<double>(n));
    spec.risk_budget = {0.3, 0.2, 0.15, 0.1, 0.1, 0.07, 0.05, 0.03};
    std::vector<double> start(n, 1.0 / static_cast<double>(n));
    double target = 0.9 * po::testing::portfolio_volatility(data, start);
    DescentProblem problem{spec, data, spec.risk_budget, VolatilityTerm{0.05, 10.0, target}};
    DescentResult full = po::riskparity::descend(problem, start, options_with(1000));
    ASSERT_EQ(full.stop, DescentStop::Converged) << "stationarity " << full.stationarity;
    ASSERT_GT(full.iterations, 1);
    double initial = merit(problem, start);
    double previous = initial;
    for (int cap = 1; cap <= full.iterations; ++cap)
    {
        SCOPED_TRACE(testing::Message() << "iteration " << cap);
        DescentResult partial = po::riskparity::descend(problem, start, options_with(cap));
        double value = merit(problem, partial.weights);
        EXPECT_LE(value, previous + row_rounding(partial));
        previous = value;
    }
    EXPECT_EQ(merit(problem, full.weights), previous);
    EXPECT_LT(previous, initial);
}

TEST(PolyhedralDescent, IterationCapIsReported)
{
    std::size_t n = 5;
    Data data = sample_data(n, 7);
    Spec spec = risk_parity_spec(0.0, 1.2 / static_cast<double>(n));
    spec.risk_budget = {0.5, 0.2, 0.15, 0.1, 0.05};
    std::vector<double> start(n, 1.0 / static_cast<double>(n));
    DescentResult result = po::riskparity::descend(DescentProblem{spec, data, spec.risk_budget, std::nullopt},
                                                   start, options_with(1));
    EXPECT_EQ(result.stop, DescentStop::IterationCap);
    EXPECT_EQ(result.iterations, 1);
    EXPECT_GT(result.stationarity, kTolerance);
}

// grad f + grad L3 = N'u at the fixed point, with L3 = -multiplier h3 + penalty / 2 h3^2: the stationarity of
// the convention with volatility multiplier multiplier - penalty h3.
TEST(PolyhedralDescent, VolatilityTermFixedPoint)
{
    std::size_t n = 5;
    Data data = sample_data(n, 23);
    Spec spec = risk_parity_spec(0.0, 1.0);
    spec.risk_budget = {0.3, 0.25, 0.2, 0.15, 0.1};
    std::vector<double> start(n, 1.0 / static_cast<double>(n));
    VolatilityTerm term{0.02, 50.0, 0.8 * po::testing::portfolio_volatility(data, start)};
    DescentProblem problem{spec, data, spec.risk_budget, term};
    DescentResult descent = po::riskparity::descend(problem, start, options_with(1000));
    ASSERT_EQ(descent.stop, DescentStop::Converged) << "stationarity " << descent.stationarity;
    EXPECT_LT(merit(problem, descent.weights), merit(problem, start));
    Result result = as_result(spec, descent);
    double gap = po::testing::portfolio_volatility(data, descent.weights) - term.target;
    result.volatility_multiplier = term.multiplier - term.penalty * gap;
    po::testing::ObjectiveGradient gradient =
        po::testing::sum_of(po::testing::risk_parity_gradient_terms(data, spec.risk_budget, descent.weights));
    po::testing::expect_user_bounds(spec, result, kBand);
    po::testing::expect_user_stationarity(data, result, gradient, kBand);
}

Data scaled_data(const Data& data, double scale)
{
    std::vector<double> covariance(data.covariance().begin(), data.covariance().end());
    for (double& entry : covariance)
    {
        entry *= scale;
    }
    return Data::create(data.asset_count(), data.mu(), covariance).value();
}

// c, J and f are homogeneous of degree 0 in Sigma, and so is every decision of the descent: daily and annual
// units give the same iterates, to the rounding of the scaled covariance.
TEST(PolyhedralDescent, ScalingSigmaLeavesIteratesUnchanged)
{
    std::size_t n = 6;
    Data annual = sample_data(n, 29);
    Data daily = scaled_data(annual, 1e-4);
    Spec spec = risk_parity_spec(0.0, 1.2 / static_cast<double>(n));
    spec.risk_budget = {0.3, 0.25, 0.15, 0.15, 0.1, 0.05};
    std::vector<double> start(n, 1.0 / static_cast<double>(n));
    for (int cap : {1, 2, 3, 5, 1000})
    {
        SCOPED_TRACE(testing::Message() << "cap " << cap);
        DescentResult first = po::riskparity::descend(
            DescentProblem{spec, annual, spec.risk_budget, std::nullopt}, start, options_with(cap));
        DescentResult second = po::riskparity::descend(
            DescentProblem{spec, daily, spec.risk_budget, std::nullopt}, start, options_with(cap));
        EXPECT_EQ(first.stop, second.stop);
        EXPECT_EQ(first.iterations, second.iterations);
        for (std::size_t i = 0; i < n; ++i)
        {
            EXPECT_NEAR(first.weights[i], second.weights[i], kBand) << "weight " << i;
        }
    }
}

// Daily units at cond(Sigma) = 1e8, where a proximal weight in the units of Sigma fell below the rounding of
// 2 J'J and the factorization of Q failed.
TEST(PolyhedralDescent, DailyScaleAtCondition1e8)
{
    using Shape = po::testing::RiskParityShape;
    Shape shape{3, 1e8, 1e-4, Shape::Budget::Equal, Shape::Box::None, false, false, 5 * 7919 + 3};
    po::testing::GeneratedRiskParity generated = po::testing::generate(shape);
    const Spec& spec = generated.spec;
    const Data& data = generated.data;
    Result minimum = po::testing::frontier_point(spec, data, std::nullopt);
    ASSERT_EQ(minimum.status, po::SolveStatus::Optimal);
    DescentResult result = po::riskparity::descend(DescentProblem{spec, data, spec.risk_budget, std::nullopt},
                                                   minimum.weights, DescentOptions{kProximal, 1e-7, 1000});
    ASSERT_EQ(result.stop, DescentStop::Converged) << "stationarity " << result.stationarity;
    po::testing::expect_risk_parity_kkt(spec, data, as_result(spec, result), 1e-7);
}

} // namespace
