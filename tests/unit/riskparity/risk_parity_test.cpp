#include "riskparity/risk_parity.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>
#include <vector>

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "fuzz/frontier_reference.hpp"
#include "fuzz/risk_parity/risk_parity_reference.hpp"
#include "random_covariance.hpp"
#include "solve/solve.hpp"

namespace
{

using po::Data;
using po::Failure;
using po::Objective;
using po::Result;
using po::SolveStatus;
using po::Spec;
using po::testing::random_covariance;
using po::testing::random_vector;

Data sample_data(std::size_t n, std::uint64_t seed)
{
    return Data::create(n, random_vector(n, seed + 500, 0.03, 0.15), random_covariance(n, seed, 1e4)).value();
}

Spec risk_parity_spec()
{
    Spec spec;
    spec.objective = Objective::RiskParity;
    return spec;
}

Result solved(const Spec& spec, const Data& data)
{
    std::expected<Result, Failure> result = po::solve(spec, data);
    EXPECT_TRUE(result.has_value());
    return result.value_or(Result{});
}

TEST(RiskParity, BudgetOnlyIsSpinu)
{
    Data data = sample_data(6, 3);
    Spec spec = risk_parity_spec();
    spec.risk_budget = {0.3, 0.2, 0.1, 0.1, 0.2, 0.1};
    Result result = solved(spec, data);
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    EXPECT_EQ(result.budget_multiplier, 0.0);
    EXPECT_EQ(result.bound_multipliers, std::vector<double>(6, 0.0));
    EXPECT_TRUE(std::isnan(result.target_return_multiplier));
    EXPECT_TRUE(std::isnan(result.volatility_multiplier));
    EXPECT_TRUE(std::isnan(result.leverage_multiplier));
    EXPECT_LE(result.residual, 1e-18);
    po::testing::expect_risk_parity_kkt(spec, data, result, 1e-10);
}

TEST(RiskParity, DefaultBudgetIsOneOverN)
{
    Data data = sample_data(4, 5);
    Spec spec = risk_parity_spec();
    EXPECT_EQ(po::riskparity::risk_budget_of(spec, data), std::vector<double>(4, 0.25));
    spec.risk_budget = {0.4, 0.3, 0.2, 0.1};
    EXPECT_EQ(po::riskparity::risk_budget_of(spec, data), spec.risk_budget);
}

TEST(RiskParity, BudgetOnlyInsideTheBoxIsSpinu)
{
    Data data = sample_data(5, 7);
    Spec spec = risk_parity_spec();
    spec.lower = 0.0;
    spec.upper = 1.0;
    Result result = solved(spec, data);
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    EXPECT_EQ(result.weights, solved(risk_parity_spec(), data).weights);
    po::testing::expect_risk_parity_kkt(spec, data, result, 1e-10);
}

TEST(RiskParity, NegativeLowerBoundKeepsPositiveOrthant)
{
    Data data = sample_data(5, 9);
    Spec bounded = risk_parity_spec();
    bounded.lower = -0.5;
    bounded.upper = 1.0;
    Result open = solved(risk_parity_spec(), data);
    Result result = solved(bounded, data);
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    EXPECT_EQ(result.weights, open.weights);
    for (double weight : result.weights)
    {
        EXPECT_GT(weight, 0.0);
    }
}

TEST(RiskParity, SingleAssetIsTheWholeBudget)
{
    Data data = Data::create(1, {}, std::vector<double>{0.04}).value();
    Result result = solved(risk_parity_spec(), data);
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    ASSERT_EQ(result.weights.size(), 1u);
    EXPECT_NEAR(result.weights[0], 1.0, 1e-15);
}

constexpr double kBand = 1e-10;
constexpr double kInfinity = std::numeric_limits<double>::infinity();

// A status without a point: NaN weights, multipliers and residual.
void expect_no_point(const Result& result, SolveStatus status, std::size_t n)
{
    EXPECT_EQ(result.status, status);
    ASSERT_EQ(result.weights.size(), n);
    ASSERT_EQ(result.bound_multipliers.size(), n);
    for (std::size_t i = 0; i < n; ++i)
    {
        EXPECT_TRUE(std::isnan(result.weights[i])) << "weight " << i;
        EXPECT_TRUE(std::isnan(result.bound_multipliers[i])) << "bound multiplier " << i;
    }
    EXPECT_TRUE(std::isnan(result.budget_multiplier));
    EXPECT_TRUE(std::isnan(result.target_return_multiplier));
    EXPECT_TRUE(std::isnan(result.volatility_multiplier));
    EXPECT_TRUE(std::isnan(result.leverage_multiplier));
    EXPECT_TRUE(std::isnan(result.residual));
}

// OPTIMAL with finite multipliers, a KKT point at the band, and the residual equal to f at the weights: the
// residual is the metric, which sums Sigma w plainly and differs from the long double reference by the
// cancellation of that product, within the band of cond(Sigma).
void expect_stationary(const Spec& spec, const Data& data, const Result& result)
{
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    po::testing::expect_risk_parity_kkt(spec, data, result, kBand);
    EXPECT_TRUE(std::isfinite(result.volatility_multiplier) || !spec.target_volatility.has_value());
    std::vector<double> budget = po::riskparity::risk_budget_of(spec, data);
    double value = po::testing::risk_parity_value(data, budget, result.weights);
    EXPECT_NEAR(result.residual, value, kBand * (1.0 + value));
    EXPECT_GT(result.iterations, 0);
}

Spec long_only_spec()
{
    Spec spec = risk_parity_spec();
    spec.lower = 0.0;
    spec.upper = 1.0;
    return spec;
}

double minimum_volatility(const Spec& spec, const Data& data)
{
    Result minimum = po::testing::frontier_point(spec, data, spec.target_return);
    EXPECT_EQ(minimum.status, SolveStatus::Optimal);
    return po::testing::portfolio_volatility(data, minimum.weights);
}

// Two assets with different means: the budget and the target return fix w, their gradients 1 and mu are
// independent, and the multipliers are the solution of [1 mu] lambda = grad f.
TEST(RiskParity, TwoAssetsWithTargetReturnIsThePoint)
{
    std::vector<double> covariance = {0.04, 0.006, 0.006, 0.09};
    Data data = Data::create(2, std::vector<double>{0.05, 0.10}, covariance).value();
    Spec spec = risk_parity_spec();
    spec.risk_budget = {0.3, 0.7};
    spec.target_return = 0.07;
    Result result = solved(spec, data);
    expect_stationary(spec, data, result);
    EXPECT_NEAR(result.weights[0], 0.6, 1e-13);
    EXPECT_NEAR(result.weights[1], 0.4, 1e-13);
    std::vector<double> point = {0.6, 0.4};
    std::vector<double> gradient = po::testing::risk_parity_gradient(data, spec.risk_budget, point);
    double target_return = (gradient[1] - gradient[0]) / 0.05;
    double budget = gradient[0] - 0.05 * target_return;
    double scale = std::fmax(std::fabs(budget), std::fabs(target_return));
    EXPECT_NEAR(result.target_return_multiplier, target_return, 1e-9 * scale);
    EXPECT_NEAR(result.budget_multiplier, budget, 1e-9 * scale);
    EXPECT_EQ(result.bound_multipliers, std::vector<double>(2, 0.0));
    EXPECT_GT(result.residual, 1e-3);
}

// One asset: the only feasible point is w = 1, where c = 1 = b and grad f = 0. With a target return equal to
// its mean the two equalities are the same row; every multiplier 0 is the valid set returned. The descent
// cannot run there (J = 0, so its curvature is 0), so this point is answered directly.
TEST(RiskParity, SingleAssetWithTargetReturnIsThePoint)
{
    Data data = Data::create(1, std::vector<double>{0.05}, std::vector<double>{0.04}).value();
    Spec spec = risk_parity_spec();
    spec.target_return = 0.05;
    for (double upper : {std::numeric_limits<double>::infinity(), 1.0})
    {
        SCOPED_TRACE(upper);
        spec.lower = -std::numeric_limits<double>::infinity();
        spec.upper = upper;
        Result result = solved(spec, data);
        ASSERT_EQ(result.status, SolveStatus::Optimal);
        EXPECT_EQ(result.weights, std::vector<double>{1.0});
        EXPECT_EQ(result.budget_multiplier, 0.0);
        EXPECT_EQ(result.target_return_multiplier, 0.0);
        EXPECT_TRUE(std::isnan(result.volatility_multiplier));
        EXPECT_EQ(result.bound_multipliers, std::vector<double>{0.0});
        EXPECT_EQ(result.residual, 0.0);
        po::testing::expect_risk_parity_kkt(spec, data, result, 1e-10);
    }
}

// The solution is a local minimum from the Spinu start: either the global minimum of the scan, or a point
// away from it with a larger value.
void expect_agrees_with_scan(const Spec& spec, const Data& data, const Result& result)
{
    po::testing::ScanResult scan = po::testing::scan_reference(spec, data, spec.risk_budget);
    ASSERT_FALSE(std::isnan(scan.value));
    double value = po::testing::risk_parity_value(data, spec.risk_budget, result.weights);
    EXPECT_GE(value, scan.value * (1.0 - 1e-9) - 1e-18);
    if (value <= scan.value * (1.0 + 1e-9))
    {
        return;
    }
    double distance = 0.0;
    for (std::size_t i = 0; i < 3; ++i)
    {
        distance = std::fmax(distance, std::fabs(result.weights[i] - scan.weights[i]));
    }
    EXPECT_GT(distance, 1e-6) << "value " << value << " scan " << scan.value;
}

TEST(RiskParity, TargetReturnLineAgainstScan)
{
    using Shape = po::testing::RiskParityShape;
    for (Shape::Box box : {Shape::Box::LongOnly, Shape::Box::Capped})
    {
        for (Shape::Budget budget : {Shape::Budget::Equal, Shape::Budget::Uneven})
        {
            for (std::uint64_t seed = 1; seed <= 10; ++seed)
            {
                SCOPED_TRACE(po::testing::describe(Shape{3, 1e4, 1.0, budget, box, true, false, seed}));
                po::testing::GeneratedRiskParity generated =
                    po::testing::generate(Shape{3, 1e4, 1.0, budget, box, true, false, seed});
                Result result = solved(generated.spec, generated.data);
                expect_stationary(generated.spec, generated.data, result);
                expect_agrees_with_scan(generated.spec, generated.data, result);
            }
        }
    }
}

TEST(RiskParity, TargetVolatilityEllipseAgainstScan)
{
    using Shape = po::testing::RiskParityShape;
    for (Shape::Box box : {Shape::Box::None, Shape::Box::LongOnly})
    {
        for (Shape::Budget budget : {Shape::Budget::Equal, Shape::Budget::Uneven})
        {
            for (std::uint64_t seed = 1; seed <= 10; ++seed)
            {
                SCOPED_TRACE(po::testing::describe(Shape{3, 1e4, 1.0, budget, box, false, true, seed}));
                po::testing::GeneratedRiskParity generated =
                    po::testing::generate(Shape{3, 1e4, 1.0, budget, box, false, true, seed});
                Result result = solved(generated.spec, generated.data);
                expect_stationary(generated.spec, generated.data, result);
                expect_agrees_with_scan(generated.spec, generated.data, result);
            }
        }
    }
}

// s = sigma_min leaves one feasible point, the minimum-variance portfolio, where the gradient of the
// volatility depends on those of the other active rows: no multipliers, and the optimum falls with infinite
// slope in s. A hair above it the problem is regular again, with a large negative volatility multiplier.
TEST(RiskParity, VolatilityEqualMinimumIsMinimumVariance)
{
    Data data = sample_data(5, 13);
    for (Spec spec : {risk_parity_spec(), long_only_spec()})
    {
        SCOPED_TRACE(spec.lower);
        Result minimum = po::testing::frontier_point(spec, data, std::nullopt);
        ASSERT_EQ(minimum.status, SolveStatus::Optimal);
        double sigma_min = po::testing::portfolio_volatility(data, minimum.weights);
        spec.target_volatility = sigma_min;
        Result result = solved(spec, data);
        ASSERT_EQ(result.status, SolveStatus::Optimal);
        EXPECT_EQ(result.weights, minimum.weights);
        EXPECT_EQ(result.volatility_multiplier, -kInfinity);
        EXPECT_TRUE(std::isnan(result.budget_multiplier));
        EXPECT_TRUE(std::isnan(result.target_return_multiplier));
        for (double beta : result.bound_multipliers)
        {
            EXPECT_TRUE(std::isnan(beta));
        }
        std::vector<double> budget = po::riskparity::risk_budget_of(spec, data);
        double value = po::testing::risk_parity_value(data, budget, minimum.weights);
        EXPECT_NEAR(result.residual, value, 1e-14 * (1.0 + value));
        po::testing::expect_risk_parity_kkt(spec, data, result, kBand);

        spec.target_volatility = sigma_min * (1.0 + 1e-9);
        Result above = solved(spec, data);
        expect_stationary(spec, data, above);
        EXPECT_LT(above.volatility_multiplier * sigma_min, -100.0);
        EXPECT_LT(above.residual, result.residual);
        for (std::size_t i = 0; i < above.weights.size(); ++i)
        {
            EXPECT_NEAR(above.weights[i], minimum.weights[i], 1e-3) << "weight " << i;
        }
    }
}

// The same at cond(Sigma) = 1e8, where Sigma w cancels by eight digits: s is the minimum volatility computed
// without that loss, as the solver computes it, and the check must not lose it either.
TEST(RiskParity, VolatilityEqualMinimumAtCondition1e8)
{
    std::size_t n = 10;
    Data data = Data::create(n, random_vector(n, 513, 0.03, 0.15), random_covariance(n, 13, 1e8)).value();
    for (Spec spec : {risk_parity_spec(), long_only_spec()})
    {
        SCOPED_TRACE(spec.lower);
        Result minimum = po::testing::frontier_point(spec, data, std::nullopt);
        ASSERT_EQ(minimum.status, SolveStatus::Optimal);
        spec.target_volatility = po::testing::wide_volatility(data, minimum.weights);
        Result result = solved(spec, data);
        ASSERT_EQ(result.status, SolveStatus::Optimal);
        EXPECT_EQ(result.weights, minimum.weights);
        EXPECT_EQ(result.volatility_multiplier, -kInfinity);
        po::testing::expect_risk_parity_kkt(spec, data, result, 1e-7);
    }
}

TEST(RiskParity, VolatilityBelowMinimumInfeasible)
{
    Data data = sample_data(5, 13);
    for (double factor : {0.5, 1.0 - 1e-9})
    {
        SCOPED_TRACE(factor);
        Spec spec = long_only_spec();
        spec.target_volatility = factor * minimum_volatility(spec, data);
        expect_no_point(solved(spec, data), SolveStatus::Infeasible, 5);
    }
}

TEST(RiskParity, TargetReturnOutOfRangeInfeasible)
{
    Data data = sample_data(5, 13);
    Spec spec = long_only_spec();
    spec.target_return = 0.16;
    expect_no_point(solved(spec, data), SolveStatus::Infeasible, 5);
    spec.target_volatility = 0.2;
    expect_no_point(solved(spec, data), SolveStatus::Infeasible, 5);
}

TEST(RiskParity, BoundsInfeasible)
{
    Data data = sample_data(4, 11);
    Spec spec = risk_parity_spec();
    spec.lower = 0.3;
    spec.upper = 1.0;
    expect_no_point(solved(spec, data), SolveStatus::Infeasible, 4);
}

// Sigma = diag(0.04, 0.09), long only: no portfolio has a volatility above 0.3, which no cheap certificate
// shows. The iteration cannot reach s = 0.5 and the answer is MAX_ITER, without a point.
TEST(RiskParity, UnreachableVolatilityIsMaxIter)
{
    Data data = Data::create(2, {}, std::vector<double>{0.04, 0.0, 0.0, 0.09}).value();
    Spec spec = risk_parity_spec();
    spec.lower = 0.0;
    spec.target_volatility = 0.5;
    Result result = solved(spec, data);
    expect_no_point(result, SolveStatus::MaxIterations, 2);
    EXPECT_GT(result.iterations, 0);
}

// Without bounds the budget and the target return leave a line, and f, homogeneous of degree 0, tends along
// it to its value at the direction of the line. On these two shapes the descent from the Spinu point follows
// f down the line without end (on the second the scan finds its minimum at the end of the scanned segment):
// the weights grow until the budget row is lost in their rounding. That is not a portfolio, and the answer is
// MAX_ITER, without a point. On the last two shapes every KKT residual of those weights is within its own
// rounding, which grows with the weights: only the budget row, no longer resolved, says they are no answer.
TEST(RiskParity, TargetReturnWithoutMinimizerIsMaxIter)
{
    using Shape = po::testing::RiskParityShape;
    constexpr Shape::Box kNone = Shape::Box::None;
    for (const Shape& shape : {Shape{3, 1e4, 1.0, Shape::Budget::Equal, kNone, true, false, 1},
                               Shape{3, 1e4, 1.0, Shape::Budget::Equal, kNone, true, false, 5},
                               Shape{5, 1e6, 1.0, Shape::Budget::Uneven, kNone, true, false, 23762},
                               Shape{10, 1e4, 1.0, Shape::Budget::Uneven, kNone, true, false, 15848}})
    {
        SCOPED_TRACE(po::testing::describe(shape));
        po::testing::GeneratedRiskParity generated = po::testing::generate(shape);
        Result result = solved(generated.spec, generated.data);
        expect_no_point(result, SolveStatus::MaxIterations, shape.n);
        EXPECT_GT(result.iterations, 0);
    }
    po::testing::GeneratedRiskParity generated =
        po::testing::generate(Shape{3, 1e4, 1.0, Shape::Budget::Equal, Shape::Box::None, true, false, 5});
    po::testing::ScanCurve curve = po::testing::return_line_curve(generated.spec, generated.data);
    po::testing::ScanResult scan =
        po::testing::scan_reference(generated.spec, generated.data, generated.spec.risk_budget);
    std::vector<double> first = curve.point(curve.first);
    std::vector<double> last = curve.point(curve.last);
    double to_first = 0.0;
    double to_last = 0.0;
    double length = 0.0;
    for (std::size_t i = 0; i < 3; ++i)
    {
        to_first = std::fmax(to_first, std::fabs(scan.weights[i] - first[i]));
        to_last = std::fmax(to_last, std::fabs(scan.weights[i] - last[i]));
        length = std::fmax(length, std::fabs(last[i] - first[i]));
    }
    // Within two of the 20000 steps of the scan from an end of the segment.
    EXPECT_LE(std::fmin(to_first, to_last), 1e-4 * length);
}

// Equal means make the target-return row a multiple of the budget row. With the common mean as the target the
// row is redundant and the solution is that of Spinu, with f = 0; with any other target no portfolio exists.
TEST(RiskParity, EqualMeansRedundantTarget)
{
    std::size_t n = 5;
    Data data = Data::create(n, std::vector<double>(n, 0.07), random_covariance(n, 9, 1e4)).value();
    Spec spec = risk_parity_spec();
    spec.risk_budget = {0.3, 0.25, 0.2, 0.15, 0.1};
    Result open = solved(spec, data);
    ASSERT_EQ(open.status, SolveStatus::Optimal);
    spec.target_return = 0.07;
    Result result = solved(spec, data);
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    for (std::size_t i = 0; i < n; ++i)
    {
        EXPECT_NEAR(result.weights[i], open.weights[i], 1e-12 * open.weights[i]) << "weight " << i;
    }
    EXPECT_LE(result.residual, 1e-18);
    po::testing::expect_risk_parity_kkt(spec, data, result, kBand);
}

TEST(RiskParity, EqualMeansInconsistentTarget)
{
    std::size_t n = 5;
    Data data = Data::create(n, std::vector<double>(n, 0.07), random_covariance(n, 9, 1e4)).value();
    Spec spec = risk_parity_spec();
    spec.target_return = 0.08;
    expect_no_point(solved(spec, data), SolveStatus::Infeasible, n);
}

// A cap of 1.2 / n and an uneven budget: the Spinu point is outside the box and the solution holds weights at
// the cap, alone and together with the two targets.
TEST(RiskParity, ActiveBoundsKkt)
{
    std::size_t n = 5;
    Data data = sample_data(n, 7);
    Spec spec = risk_parity_spec();
    spec.lower = 0.0;
    spec.upper = 1.2 / static_cast<double>(n);
    spec.risk_budget = {0.5, 0.2, 0.15, 0.1, 0.05};
    Result capped = solved(spec, data);
    expect_stationary(spec, data, capped);
    EXPECT_TRUE(std::ranges::any_of(capped.bound_multipliers, [](double beta) { return beta < 0.0; }));
    EXPECT_GT(capped.residual, 1e-4);
    spec.target_return = po::testing::portfolio_return(data, capped.weights) + 0.002;
    Result with_return = solved(spec, data);
    expect_stationary(spec, data, with_return);
    EXPECT_TRUE(std::ranges::any_of(with_return.bound_multipliers, [](double beta) { return beta != 0.0; }));
    spec.target_volatility = 1.05 * minimum_volatility(spec, data);
    Result with_both = solved(spec, data);
    expect_stationary(spec, data, with_both);
    EXPECT_TRUE(std::ranges::any_of(with_both.bound_multipliers, [](double beta) { return beta != 0.0; }));
}

// The contributions are homogeneous of degree 0, so w'grad f = 0, which at a KKT point reads
// budget = -target_return r - volatility s - sum_i beta_i w_i.
TEST(RiskParity, HomogeneityIdentity)
{
    std::size_t n = 6;
    Data data = sample_data(n, 17);
    Spec spec = risk_parity_spec();
    spec.lower = 0.0;
    spec.upper = 0.3;
    spec.risk_budget = {0.4, 0.2, 0.15, 0.1, 0.1, 0.05};
    Result capped = solved(spec, data);
    ASSERT_EQ(capped.status, SolveStatus::Optimal);
    spec.target_return = po::testing::portfolio_return(data, capped.weights) + 0.003;
    spec.target_volatility = 1.1 * minimum_volatility(spec, data);
    Result result = solved(spec, data);
    expect_stationary(spec, data, result);
    double return_term = result.target_return_multiplier * *spec.target_return;
    double volatility_term = result.volatility_multiplier * *spec.target_volatility;
    double bound_term = po::testing::dot(result.bound_multipliers, result.weights);
    double scale = std::fmax(std::fabs(result.budget_multiplier), std::fabs(return_term));
    scale = std::fmax(scale, std::fmax(std::fabs(volatility_term), std::fabs(bound_term)));
    ASSERT_GT(scale, 0.0);
    EXPECT_NE(bound_term, 0.0);
    EXPECT_NEAR(result.budget_multiplier, -return_term - volatility_term - bound_term, 1e-9 * scale);
}

// Sigma -> t Sigma with s -> sqrt(t) s is a change of units: the same weights, and a volatility multiplier
// divided by sqrt(t).
TEST(RiskParity, ScaleInvariance)
{
    std::size_t n = 6;
    constexpr double kScale = 1e-4;
    Data annual = sample_data(n, 17);
    std::vector<double> covariance(annual.covariance().begin(), annual.covariance().end());
    for (double& entry : covariance)
    {
        entry *= kScale;
    }
    Data daily = Data::create(n, annual.mu(), covariance).value();
    Spec spec = long_only_spec();
    spec.upper = 0.3;
    spec.risk_budget = {0.4, 0.2, 0.15, 0.1, 0.1, 0.05};
    spec.target_volatility = 1.1 * minimum_volatility(spec, annual);
    Result first = solved(spec, annual);
    expect_stationary(spec, annual, first);
    Spec scaled = spec;
    scaled.target_volatility = std::sqrt(kScale) * *spec.target_volatility;
    Result second = solved(scaled, daily);
    expect_stationary(scaled, daily, second);
    for (std::size_t i = 0; i < n; ++i)
    {
        EXPECT_NEAR(first.weights[i], second.weights[i], kBand) << "weight " << i;
    }
    EXPECT_NEAR(second.volatility_multiplier * std::sqrt(kScale), first.volatility_multiplier,
                1e-8 * std::fabs(first.volatility_multiplier));
    EXPECT_NEAR(second.residual, first.residual, 1e-9 * first.residual);
}

double optimal_value(Spec spec, const Data& data, double target_return, double target_volatility)
{
    spec.target_return = target_return;
    spec.target_volatility = target_volatility;
    Result result = solved(spec, data);
    EXPECT_EQ(result.status, SolveStatus::Optimal);
    return po::testing::risk_parity_value(data, po::riskparity::risk_budget_of(spec, data), result.weights);
}

// Each multiplier is the derivative of the optimal f in the right-hand side of its row: central differences
// with a relative step of 1e-6 in r and in s.
TEST(RiskParity, MultipliersAreDerivatives)
{
    constexpr double kStep = 1e-6;
    constexpr double kTolerance = 1e-5;
    std::size_t n = 6;
    for (std::uint64_t seed : {17u, 31u, 43u})
    {
        SCOPED_TRACE(seed);
        Data data = sample_data(n, seed);
        Spec spec = long_only_spec();
        spec.upper = 0.3;
        spec.risk_budget = {0.4, 0.2, 0.15, 0.1, 0.1, 0.05};
        Result capped = solved(spec, data);
        ASSERT_EQ(capped.status, SolveStatus::Optimal);
        double r = po::testing::portfolio_return(data, capped.weights) + 0.003;
        spec.target_return = r;
        double s = 1.1 * minimum_volatility(spec, data);
        spec.target_volatility = s;
        Result result = solved(spec, data);
        expect_stationary(spec, data, result);
        double return_slope = (optimal_value(spec, data, r * (1.0 + kStep), s) -
                               optimal_value(spec, data, r * (1.0 - kStep), s)) /
                              (2.0 * kStep * r);
        double volatility_slope = (optimal_value(spec, data, r, s * (1.0 + kStep)) -
                                   optimal_value(spec, data, r, s * (1.0 - kStep))) /
                                  (2.0 * kStep * s);
        EXPECT_NEAR(result.target_return_multiplier, return_slope,
                    kTolerance * std::fabs(result.target_return_multiplier));
        EXPECT_NEAR(result.volatility_multiplier, volatility_slope,
                    kTolerance * std::fabs(result.volatility_multiplier));
        EXPECT_NE(result.target_return_multiplier, 0.0);
        EXPECT_NE(result.volatility_multiplier, 0.0);
    }
}

} // namespace
