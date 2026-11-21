#include "riskparity/risk_parity.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>
#include <string>
#include <vector>

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "fuzz/frontier_reference.hpp"
#include "fuzz/risk_parity/risk_parity_reference.hpp"
#include "fuzz/risk_parity/risk_parity_refinement.hpp"
#include "metrics/portfolio.hpp"
#include "random_covariance.hpp"
#include "riskparity/kkt.hpp"
#include "riskparity/newton_polish.hpp"
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

// The solution is a local minimum from the Spinu start: either the global minimum of the scan (its value moved to
// the feasible set, scan_value_at_target), or a point away from it with a larger value.
void expect_agrees_with_scan(const Spec& spec, const Data& data, const Result& result)
{
    po::testing::ScanResult scan = po::testing::scan_reference(spec, data, spec.risk_budget);
    ASSERT_FALSE(std::isnan(scan.value));
    double value = po::testing::risk_parity_value(data, spec.risk_budget, result.weights);
    EXPECT_GE(value, scan.value * (1.0 - 1e-9) - 1e-18);
    if (value <= po::testing::scan_value_at_target(spec, data, result, scan) * (1.0 + 1e-9))
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
        // The residual is f at those weights as po::metrics computes it, in double. Against f in long double it
        // is off by the rounding of Sigma w, whose first-order bound is 3e-12 of f at the minimum-variance point
        // without a box (3.3e-14 measured).
        std::vector<double> budget = po::riskparity::risk_budget_of(spec, data);
        EXPECT_EQ(result.residual, po::metrics::risk_budget_residual(data, minimum.weights, budget));
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

po::riskparity::KktPoint kkt_point_of(const Spec& spec, const Result& result)
{
    po::riskparity::KktPoint point{result.weights, {result.budget_multiplier}, result.bound_multipliers};
    if (spec.target_return.has_value())
    {
        point.equality_multipliers.push_back(result.target_return_multiplier);
    }
    if (spec.target_volatility.has_value())
    {
        point.equality_multipliers.push_back(result.volatility_multiplier);
    }
    return point;
}

// tolerance relative to the largest expected multiplier.
void expect_near_multipliers(const std::vector<double>& actual, const std::vector<double>& expected,
                             double tolerance)
{
    double scale = po::testing::largest_magnitude(expected);
    for (std::size_t k = 0; k < expected.size(); ++k)
    {
        EXPECT_NEAR(actual[k], expected[k], tolerance * scale) << "multiplier " << k;
    }
}

// The rounding of the least-squares fit of a polish away from sigma_min, relative to the largest multiplier.
constexpr double kRefitTolerance = 1e-12;

// The least-squares point of a polish of a result and the correction of its multipliers there.
struct Polished
{
    po::riskparity::KktPoint least_squares;
    po::riskparity::MultiplierCorrection correction;
};

// The polish of a result settles the same weights again (no step): its least-squares multipliers are those
// the KKT gate checked, and its correction the multipliers po::solve returned, both to the rounding of the
// least-squares fit, which depends on the multipliers it starts from: tolerance, relative to the largest.
Polished polish_of(const Spec& spec, const Data& data, const Result& result, double tolerance)
{
    std::vector<double> budget = po::riskparity::risk_budget_of(spec, data);
    po::riskparity::KktProblem problem{spec, data, budget};
    po::riskparity::PolishResult polished =
        po::riskparity::newton_polish(problem, kkt_point_of(spec, result));
    EXPECT_EQ(polished.steps, 0);
    Polished settled{polished.least_squares_point,
                     po::riskparity::corrected_multipliers(problem, polished.least_squares_point)};
    expect_near_multipliers(settled.correction.reported_point.equality_multipliers,
                            kkt_point_of(spec, result).equality_multipliers, tolerance);
    expect_near_multipliers(settled.correction.reported_point.bound_multipliers, result.bound_multipliers,
                            tolerance);
    return settled;
}

// The rounding of the least-squares fit of the multipliers at s = sigma_min (1 + delta), relative to the largest.
// There w = w_min + u with u'Sigma u = s^2 - sigma_min^2, about 2 delta sigma_min^2, and Sigma w_min is
// sigma_min^2 1: the gradient of the volatility, Sigma w / s, leaves the budget row 1 only by Sigma u / s, of
// relative size sqrt(2 delta). The two columns of the fit are that close to parallel, so its condition number
// grows like 1 + 1 / sqrt(2 delta), and the fit rounds to n eps times that. Measured: at most 0.47 eps (1 +
// 1 / sqrt(2 delta)) (delta = 1e-2), 0.26 at 1e-6 and 0.08 at 1e-10, all ten times below the band. The band is
// 7.9e-13 at delta = 1e-6, tighter than the 1e-12 it replaces, and 7.9e-11 at 1e-10, where the fit itself
// rounds to 1.2e-12.
double refit_tolerance(std::size_t n, double delta)
{
    return static_cast<double>(n) * std::numeric_limits<double>::epsilon() * (1.0 + 1.0 / std::sqrt(2.0 * delta));
}

// s = sigma_min (1 + delta), where the volatility multiplier grows like 1 / sqrt(delta) and the gradients of
// the budget and of the volatility turn parallel: the Newton correction of the multipliers is returned only
// where the safeguard accepts it. Every result is OPTIMAL and a KKT point at the band, whose returned
// multipliers expect_risk_parity_kkt finds no farther from the refined KKT point than the least-squares ones.
// The safeguard keeps the least-squares multipliers at delta from 1e-7 to 1e-10 without a box and at 1e-6 and
// 1e-7 long only: 6 of the 18, measured. Which ones it keeps turns on single roundings: its contraction test
// (second correction at most half the first) is near its threshold here, at 0.43 long only at 1e-10.
TEST(RiskParity, NearMinimumVolatilityKeepsLeastSquaresMultipliers)
{
    Data data = sample_data(5, 13);
    int kept = 0;
    for (Spec spec : {risk_parity_spec(), long_only_spec()})
    {
        double sigma_min = minimum_volatility(spec, data);
        for (double delta : {1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10})
        {
            SCOPED_TRACE(::testing::Message() << "lower " << spec.lower << " delta " << delta);
            spec.target_volatility = sigma_min * (1.0 + delta);
            Result result = solved(spec, data);
            expect_stationary(spec, data, result);
            Polished polished = polish_of(spec, data, result, refit_tolerance(data.asset_count(), delta));
            kept += polished.correction.accepted ? 0 : 1;
            // Where the correction, though it contracts, would leave a stationarity far above its noise
            // (docs/guide/validation.md): the safeguard keeps the least-squares multipliers.
            if (delta == 1e-10 && std::isinf(spec.lower))
            {
                EXPECT_FALSE(polished.correction.accepted);
            }
        }
    }
    RecordProperty("kept_least_squares", kept);
    EXPECT_EQ(kept, 6);
}

// The path with a target return and no target volatility, at cond(Sigma) = 1e8: the correction of the
// multipliers never leaves them farther from the refined KKT point than the least-squares fit it starts from,
// beyond 1e-14. Measured: the fit is off by up to 1.3e-11 (long-only, uneven budget, seed 2), the corrected
// multipliers by up to 3.1e-14 there; elsewhere the fit is within 5.0e-14 and the correction within 1.2e-14.
TEST(RiskParity, TargetReturnOnlyCorrectionNeverWorsensTheFit)
{
    if (!po::testing::kLongDoubleIsExtended)
    {
        GTEST_SKIP() << "long double is double here: the refinement has no extra precision";
    }
    using Shape = po::testing::RiskParityShape;
    int same_bits = 0;
    double worst = 0.0;
    double worst_least_squares = 0.0;
    for (Shape::Box box : {Shape::Box::LongOnly, Shape::Box::Capped})
    {
        for (Shape::Budget budget : {Shape::Budget::Equal, Shape::Budget::Uneven})
        {
            for (std::uint64_t seed = 1; seed <= 4; ++seed)
            {
                Shape shape{8, 1e8, 1.0, budget, box, true, false, seed};
                SCOPED_TRACE(po::testing::describe(shape));
                po::testing::GeneratedRiskParity generated = po::testing::generate(shape);
                Result result = solved(generated.spec, generated.data);
                ASSERT_EQ(result.status, SolveStatus::Optimal);
                Polished polished = polish_of(generated.spec, generated.data, result, kRefitTolerance);
                Result least_squares = result;
                const po::riskparity::KktPoint& fitted = polished.least_squares;
                least_squares.budget_multiplier = fitted.equality_multipliers.front();
                least_squares.target_return_multiplier = fitted.equality_multipliers.back();
                least_squares.bound_multipliers = fitted.bound_multipliers;
                bool same = least_squares.budget_multiplier == result.budget_multiplier &&
                            least_squares.target_return_multiplier == result.target_return_multiplier &&
                            least_squares.bound_multipliers == result.bound_multipliers;
                double error = po::testing::forward_error(generated.spec, generated.data, result).multipliers;
                double fit_error =
                    po::testing::forward_error(generated.spec, generated.data, least_squares).multipliers;
                EXPECT_LE(error, std::fmax(1e-14, fit_error)) << "fit error " << fit_error;
                same_bits += same ? 1 : 0;
                worst = std::fmax(worst, error);
                worst_least_squares = std::fmax(worst_least_squares, fit_error);
            }
        }
    }
    RecordProperty("same_bits", same_bits);
    RecordProperty("worst", ::testing::PrintToString(worst));
    RecordProperty("worst_least_squares", ::testing::PrintToString(worst_least_squares));
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
// shows. The iteration cannot reach s = 0.5 and the answer is NOT_CONVERGED, without a point.
TEST(RiskParity, UnreachableVolatilityIsNotConverged)
{
    Data data = Data::create(2, {}, std::vector<double>{0.04, 0.0, 0.0, 0.09}).value();
    Spec spec = risk_parity_spec();
    spec.lower = 0.0;
    spec.target_volatility = 0.5;
    Result result = solved(spec, data);
    expect_no_point(result, SolveStatus::NotConverged, 2);
    EXPECT_GT(result.iterations, 0);
}

// The scan's minimum within two of the 20000 steps of the scan from an end of the segment.
void expect_scan_minimum_at_an_end(const po::testing::GeneratedRiskParity& generated)
{
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
    EXPECT_LE(std::fmin(to_first, to_last), 1e-4 * length);
}

// Without bounds the budget and the target return leave a line, and f, homogeneous of degree 0, tends along
// it to its value at the direction of the line. On these two shapes the descent from the Spinu point follows
// f down the line without end (on both the scan finds its minimum at an end of the scanned segment): the
// weights grow until the budget row is lost in their rounding. That is not a portfolio, and the answer is
// NOT_CONVERGED, without a point. On the last two shapes every KKT residual of those weights is within its
// own rounding, which grows with the weights: only the budget row, no longer resolved, says they are no
// answer.
// The seeds are the first of each kind with those properties, found by search.
TEST(RiskParity, TargetReturnWithoutMinimizerIsNotConverged)
{
    using Shape = po::testing::RiskParityShape;
    constexpr Shape::Box kNone = Shape::Box::None;
    for (const Shape& shape : {Shape{3, 1e4, 1.0, Shape::Budget::Equal, kNone, true, false, 11},
                               Shape{3, 1e4, 1.0, Shape::Budget::Equal, kNone, true, false, 129},
                               Shape{5, 1e6, 1.0, Shape::Budget::Uneven, kNone, true, false, 32},
                               Shape{10, 1e4, 1.0, Shape::Budget::Uneven, kNone, true, false, 33}})
    {
        SCOPED_TRACE(po::testing::describe(shape));
        po::testing::GeneratedRiskParity generated = po::testing::generate(shape);
        Result result = solved(generated.spec, generated.data);
        expect_no_point(result, SolveStatus::NotConverged, shape.n);
        EXPECT_GT(result.iterations, 0);
        if (shape.n == 3)
        {
            expect_scan_minimum_at_an_end(generated);
        }
    }
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
// budget = -target_return r - volatility s - sum_i beta_i w_i. At cond(Sigma) 1e4 and 1e6, with the
// multipliers the Newton correction of the polish returns, the bound multipliers of the held weights (box
// [0, 0.3]) included; seed 1 is the first where the box holds a weight at both. The tolerance is the largest
// relative gap measured, 4.8e-17 at 1e4 (4.4e-17 at 1e6), times 10: measured, not derived. On the instances
// before the portable generators the least-squares multipliers at the same weights left gaps of the same order
// as the corrected ones, so the identity checks that the corrected beta are consistent to rounding, not that
// they are nearer the KKT point (RiskParityFuzz.TargetVolatilityBoundMultipliersByCondition does).
TEST(RiskParity, HomogeneityIdentity)
{
    std::size_t n = 6;
    for (double condition : {1e4, 1e6})
    {
        SCOPED_TRACE(condition);
        Data data =
            Data::create(n, random_vector(n, 501, 0.03, 0.15), random_covariance(n, 1, condition)).value();
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
        double gap = result.budget_multiplier + return_term + volatility_term + bound_term;
        RecordProperty("relative_gap_" + std::to_string(static_cast<int>(std::log10(condition))),
                       ::testing::PrintToString(std::fabs(gap) / scale));
        EXPECT_NEAR(result.budget_multiplier, -return_term - volatility_term - bound_term, 4.8e-16 * scale);
    }
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

std::vector<double> scaled_covariance(double scale)
{
    std::vector<double> covariance = {0.04, 0.006, 0.002, 0.006, 0.09, 0.01, 0.002, 0.01, 0.16};
    for (double& entry : covariance)
    {
        entry *= scale;
    }
    return covariance;
}

// Scale invariance at the ends of the doubles: w'Sigma w of 1e160 and 1e-200 is past 1e154 and 1e-154, where a
// product of two variances leaves the doubles. A cap of 0.4 holds the first weight; the weights are those at
// Sigma unscaled.
TEST(RiskParity, ExtremeScalesWithBox)
{
    Spec spec = long_only_spec();
    spec.upper = 0.4;
    Result unscaled = solved(spec, Data::create(3, {}, scaled_covariance(1.0)).value());
    ASSERT_EQ(unscaled.status, SolveStatus::Optimal);
    ASSERT_LT(unscaled.bound_multipliers[0], 0.0);
    for (double scale : {1e160, 1e-200})
    {
        SCOPED_TRACE(scale);
        Data data = Data::create(3, {}, scaled_covariance(scale)).value();
        Result result = solved(spec, data);
        expect_stationary(spec, data, result);
        for (std::size_t i = 0; i < 3; ++i)
        {
            EXPECT_NEAR(result.weights[i], unscaled.weights[i], kBand) << "weight " << i;
        }
    }
}

// A budget entry of 1e-200 squares to 0 in the barrier curvature b_i / y_i^2 unless divided twice.
TEST(RiskParity, TinyBudgetEntryIsOptimal)
{
    Data data = Data::create(3, {}, scaled_covariance(1.0)).value();
    Spec spec = risk_parity_spec();
    spec.risk_budget = {1e-200, 0.5, 0.5};
    Result result = solved(spec, data);
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    EXPECT_TRUE(std::ranges::all_of(result.weights, [](double weight) { return weight > 0.0; }));
    EXPECT_LE(result.residual, 1e-18);
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
