#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "fuzz/risk_parity/risk_parity_reference.hpp"
#include "fuzz/risk_parity/risk_parity_refinement.hpp"
#include "solve/solve.hpp"

namespace
{

using po::Result;
using po::SolveStatus;
using Shape = po::testing::RiskParityShape;

double band_of(double condition)
{
    if (condition <= 1e4)
    {
        return 1e-10;
    }
    if (condition <= 1e6)
    {
        return 1e-9;
    }
    return 1e-7;
}

bool is_all_nan(const Result& result)
{
    bool all = std::isnan(result.budget_multiplier) && std::isnan(result.target_return_multiplier) &&
               std::isnan(result.volatility_multiplier) && std::isnan(result.residual);
    auto is_nan = [](double value) { return std::isnan(value); };
    return all && std::ranges::all_of(result.weights, is_nan) &&
           std::ranges::all_of(result.bound_multipliers, is_nan);
}

// Without a point: NOT_CONVERGED or INFEASIBLE, every number NaN. With one: KKT at the band of cond(Sigma).
void expect_never_silently_wrong(const po::Spec& spec, const po::Data& data, const Result& result,
                                 double condition)
{
    if (result.status != SolveStatus::Optimal)
    {
        EXPECT_TRUE(result.status == SolveStatus::NotConverged || result.status == SolveStatus::Infeasible);
        EXPECT_TRUE(is_all_nan(result));
        return;
    }
    po::testing::expect_risk_parity_kkt(spec, data, result, band_of(condition));
}

double largest_distance(const std::vector<double>& a, const std::vector<double>& b)
{
    double distance = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        distance = std::fmax(distance, std::fabs(a[i] - b[i]));
    }
    return distance;
}

// The OPTIMAL results of AgainstScanReference that are another local minimum of the curve, above the global
// one of the scan, as measured (docs/guide/validation.md): more is a regression.
constexpr int kOtherLocalMinima = 18;

// No feasible scan point of the curve within 20 steps of the result (of 20001 points; 40 of the 40001 on the
// ellipse, the same share of its parameter) has a smaller f: the result is a local minimum along the curve, not a
// saddle or a maximum of it.
bool is_local_minimum_on_curve(const po::Spec& spec, const po::Data& data, const Result& result, double value)
{
    po::testing::ScanCurve curve = po::testing::scan_curve_of(spec, data).value();
    int points = static_cast<int>(curve.points);
    int window = 20 * (points - 1) / static_cast<int>(po::testing::kScanPoints - 1);
    double step = (curve.last - curve.first) / (points - 1);
    int nearest = 0;
    double nearest_distance = std::numeric_limits<double>::infinity();
    for (int k = 0; k < points; ++k)
    {
        double distance = largest_distance(curve.point(curve.first + step * k), result.weights);
        nearest = distance < nearest_distance ? k : nearest;
        nearest_distance = std::fmin(distance, nearest_distance);
    }
    for (int k = std::max(0, nearest - window); k <= std::min(points - 1, nearest + window); ++k)
    {
        po::testing::ScanResult sample =
            po::testing::scan_point(spec, data, spec.risk_budget, curve, curve.first + step * k);
        if (sample.value < value * (1.0 - 1e-9))
        {
            return false;
        }
    }
    return true;
}

struct ScanTally
{
    int compared = 0;
    int other_minima = 0;
};

void compare_with_scan(const Shape& shape, ScanTally& tally)
{
    SCOPED_TRACE(po::testing::describe(shape));
    po::testing::GeneratedRiskParity generated = po::testing::generate(shape);
    const po::Spec& spec = generated.spec;
    Result result = po::solve(spec, generated.data).value();
    expect_never_silently_wrong(spec, generated.data, result, shape.condition);
    po::testing::ScanResult scan = po::testing::scan_reference(spec, generated.data, spec.risk_budget);
    if (result.status != SolveStatus::Optimal || std::isnan(scan.value))
    {
        return;
    }
    ++tally.compared;
    double value = po::testing::risk_parity_value(generated.data, spec.risk_budget, result.weights);
    if (value <= po::testing::scan_value_at_target(spec, generated.data, result, scan) * (1.0 + 1e-9))
    {
        return;
    }
    ++tally.other_minima;
    EXPECT_GT(largest_distance(result.weights, scan.weights), 1e-6)
        << "value " << value << " scan " << scan.value;
    EXPECT_TRUE(is_local_minimum_on_curve(spec, generated.data, result, value));
}

// cond(Sigma) 1e2 to 1e8, the three boxes, a target return or a target volatility, equal and uneven budgets.
std::vector<Shape> scan_shapes()
{
    std::vector<Shape> shapes;
    for (double condition : {1e2, 1e4, 1e6, 1e8})
    {
        for (Shape::Box box : {Shape::Box::None, Shape::Box::LongOnly, Shape::Box::Capped})
        {
            for (std::uint64_t seed = 1; seed <= 8; ++seed)
            {
                Shape::Budget budget = seed % 2 == 0 ? Shape::Budget::Uneven : Shape::Budget::Equal;
                bool target_return = seed <= 4;
                shapes.push_back(
                    Shape{3, condition, 1.0, budget, box, target_return, !target_return, 101 + seed});
            }
        }
    }
    return shapes;
}

// n = 3 with a target return or a target volatility, where the feasible set is a curve that scan_reference
// scans for its global minimum. f of the result is at most the scan's (1 + 1e-9), unless the result is
// another local minimum: a KKT point at the band (expect_risk_parity_kkt), away from the scan's point, with
// no smaller f on the curve near it. Those are counted.
TEST(RiskParityFuzz, AgainstScanReference)
{
    ScanTally tally;
    for (const Shape& shape : scan_shapes())
    {
        compare_with_scan(shape, tally);
    }
    RecordProperty("compared", tally.compared);
    RecordProperty("other_minima", tally.other_minima);
    EXPECT_GE(tally.compared, 60);
    EXPECT_LE(tally.other_minima, kOtherLocalMinima);
}

// b_i log-spaced from 1e-6 to 1, normalized.
std::vector<double> extreme_budget(std::size_t n)
{
    std::vector<double> budget(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        budget[i] = po::testing::exp10(-6.0 + 6.0 * static_cast<double>(i) / static_cast<double>(n - 1));
    }
    double total = std::accumulate(budget.begin(), budget.end(), 0.0);
    for (double& entry : budget)
    {
        entry /= total;
    }
    return budget;
}

// Budgets from 1e-6 to 1 with cond(Sigma) up to 1e8, at two scales, the budget row alone and with the Capped
// box. Every OPTIMAL is KKT at the band; any other status has no point.
TEST(RiskParityFuzz, UnevenBudgetsAndConditioning)
{
    constexpr std::size_t kSizes[] = {5, 20};
    constexpr double kConditions[] = {1e2, 1e4, 1e6, 1e8};
    constexpr double kScales[] = {1.0, 1e-4};
    constexpr Shape::Box kBoxes[] = {Shape::Box::None, Shape::Box::Capped};
    int optimal = 0;
    // Every combination of the four lists above and seeds 1 to 3, one index each.
    for (std::size_t k = 0; k < 2 * 4 * 2 * 2 * 3; ++k)
    {
        std::size_t n = kSizes[k / 48];
        double condition = kConditions[k / 12 % 4];
        Shape::Box box = kBoxes[k / 3 % 2];
        Shape shape{n, condition, kScales[k / 6 % 2], Shape::Budget::Uneven, box, false, false, 212 + k % 3};
        SCOPED_TRACE(po::testing::describe(shape));
        po::testing::GeneratedRiskParity generated = po::testing::generate(shape);
        generated.spec.risk_budget = extreme_budget(n);
        Result result = po::solve(generated.spec, generated.data).value();
        expect_never_silently_wrong(generated.spec, generated.data, result, condition);
        optimal += result.status == SolveStatus::Optimal;
    }
    RecordProperty("optimal", optimal);
    EXPECT_EQ(optimal, 96);
}

// The descents of the polyhedral-descent sweep that the polish did not turn into a KKT point on the instances before
// plan 8, through po::solve: OPTIMAL and KKT at the band, or no point (docs/guide/validation.md). The shape of seed
// 31696 was NOT_CONVERGED at both scales then; on the present instances both are OPTIMAL.
TEST(RiskParityFuzz, UnrepairedDescentShapesThroughSolve)
{
    for (double scale : {1.0, 1e-4})
    {
        Shape shape{20, 1e8, scale, Shape::Budget::Equal, Shape::Box::Capped, false, false, 31696};
        SCOPED_TRACE(po::testing::describe(shape));
        po::testing::GeneratedRiskParity generated = po::testing::generate(shape);
        Result result = po::solve(generated.spec, generated.data).value();
        RecordProperty("status_" + std::to_string(static_cast<int>(scale == 1.0)),
                       static_cast<int>(result.status));
        expect_never_silently_wrong(generated.spec, generated.data, result, shape.condition);
    }
}

// The guard of the weights of this class: they are not held to the band of cond(Sigma), which does not apply
// at the floor, but their error must not grow past this measured ceiling (the largest of 229 results is
// 2.8e-5; docs/guide/validation.md). Measured, not derived.
constexpr double kWeightErrorCeiling = 1e-4;
// A sanity bound only: weight_error_bound holds but is at least 5e2 times above the error (the largest error over
// the bound, in 229 results, is 1.9e-3), so it guards nothing the ceiling does not; it would catch a refinement or a
// bound that came apart.
constexpr double kSanityBoundFactor = 4.0;

// What the forward error of one OPTIMAL result is, against the refined KKT point of its active set.
struct WeightCheck
{
    po::testing::ForwardError forward;
    po::testing::WeightErrorBound bound;
    double value;
    bool floor_limited;
};

WeightCheck weight_check(const po::testing::GeneratedRiskParity& generated, const Result& result, double band)
{
    const po::Spec& spec = generated.spec;
    WeightCheck check{po::testing::forward_error(spec, generated.data, result),
                      po::testing::weight_error_bound(spec, generated.data, result),
                      po::testing::risk_parity_value(generated.data, spec.risk_budget, result.weights),
                      false};
    check.floor_limited =
        po::testing::risk_parity_stationarity(spec, generated.data, result, band).relative > band;
    return check;
}

Shape thinnest_class(std::size_t n, std::uint64_t seed)
{
    return Shape{n, 1e8, 1.0, Shape::Budget::Uneven, Shape::Box::None, false, true, 7919 * seed + n};
}

// The thinnest forward-error class of the sweep: no box, a target volatility, cond(Sigma) = 1e8, n = 20.
// Seeds 1 to 4 have weight errors from 9.8e-14 to 6.2e-12 (seed 4, 31696, the largest); over 199 seeds of n = 20
// the worst is 1.3e-5, seed 135, at a stationarity far below its rounding floor. The value is right at the band and
// the weights below their measured ceiling, kWeightErrorCeiling; the band of the weights does not apply to a result
// at its floor.
TEST(RiskParityFuzz, WeightsAtCondition1e8BelowTheirMeasuredCeiling)
{
    if (!po::testing::kLongDoubleIsExtended)
    {
        GTEST_SKIP() << "long double is double here: the refinement has no extra precision";
    }
    for (std::uint64_t seed = 1; seed <= 4; ++seed)
    {
        Shape shape = thinnest_class(20, seed);
        SCOPED_TRACE(po::testing::describe(shape));
        po::testing::GeneratedRiskParity generated = po::testing::generate(shape);
        Result result = po::solve(generated.spec, generated.data).value();
        ASSERT_EQ(result.status, SolveStatus::Optimal);
        po::testing::expect_risk_parity_kkt(generated.spec, generated.data, result, 1e-7);
        WeightCheck check = weight_check(generated, result, 1e-7);
        RecordProperty("seed_" + std::to_string(shape.seed),
                       ::testing::PrintToString(std::vector<double>{check.forward.weights, check.bound.floor,
                                                                    check.bound.residual, check.forward.value,
                                                                    check.value}));
        EXPECT_LE(check.forward.reference_over_floor, 1e-2);
        EXPECT_LE(check.forward.value, 1e-7 * (1.0 + check.value));
        EXPECT_LE(check.forward.weights, kWeightErrorCeiling);
        EXPECT_LE(check.forward.weights, kSanityBoundFactor * check.bound.residual)
            << "sanity bound " << check.bound.residual;
    }
}

// The ceilings of the multiplier errors of the target-volatility path at one cond(Sigma): the largest
// measured on the cases of the test that uses them times 10, measured, not derived
// (docs/guide/validation.md).
struct MultiplierCeiling
{
    double condition;
    // forward_error(...).multipliers: max_k |lambda_k - lambda_ref_k| / max_k |lambda_ref_k|.
    double equalities;
    // bound_multiplier_error; 0 where the test has no bound multiplier.
    double bounds;
};

// The multipliers of the target-volatility path against those of the KKT point refined in long double from
// the result, by cond(Sigma): no box, equal budget, n = 5, 8, 12 and 20. The polish keeps the Newton
// correction of the multipliers at the weights it returns. Measured largest: 8.1e-13 at 1e4 (n = 8, seed 405),
// 2.4e-12 at 1e6 and 9.3e-11 at 1e8 (both n = 8, seed 404). On the instances before the portable generators they
// were 2.0e-14, 5.2e-11 and 5.3e-10, all at n = 5 (4.8e-10, 3.3e-7 and 3.7e-2 with the least-squares multipliers
// the polish returned before); there, at 1e8, the scaled KKT matrix of n = 12 (seed 404) and n = 20 (seeds 403
// to 405) failed the rank test of the QR, and with that test as a veto they kept least-squares multipliers wrong
// by up to 4.1e-2, corrected within 2.9e-11.
constexpr MultiplierCeiling kMultiplierCeilings[] = {
    {1e4, 8.2e-12, 0.0}, {1e6, 2.4e-11, 0.0}, {1e8, 9.4e-10, 0.0}};

TEST(RiskParityFuzz, TargetVolatilityMultipliersByCondition)
{
    if (!po::testing::kLongDoubleIsExtended)
    {
        GTEST_SKIP() << "long double is double here: the refinement has no extra precision";
    }
    for (MultiplierCeiling ceiling : kMultiplierCeilings)
    {
        double worst = 0.0;
        for (std::size_t n : {5u, 8u, 12u, 20u})
        {
            for (std::uint64_t seed = 1; seed <= 4; ++seed)
            {
                Shape shape{n, ceiling.condition, 1.0, Shape::Budget::Equal, Shape::Box::None,
                            false, true, 401 + seed};
                SCOPED_TRACE(po::testing::describe(shape));
                po::testing::GeneratedRiskParity generated = po::testing::generate(shape);
                Result result = po::solve(generated.spec, generated.data).value();
                ASSERT_EQ(result.status, SolveStatus::Optimal);
                po::testing::ForwardError forward =
                    po::testing::forward_error(generated.spec, generated.data, result);
                EXPECT_LE(forward.reference_over_floor, 1e-2);
                EXPECT_LE(forward.multipliers, ceiling.equalities);
                worst = std::fmax(worst, forward.multipliers);
            }
        }
        RecordProperty("worst_" + std::to_string(static_cast<int>(std::log10(ceiling.condition))),
                       ::testing::PrintToString(worst));
    }
}

// The bound multipliers of the held weights of a result (beta != 0) against those of the KKT point refined in
// long double from it, where beta_ref_j = (grad f - sum_k lambda_k grad h_k)_j is what row j leaves:
// max_j abs(beta_j - beta_ref_j) over max(abs(lambda_budget*), abs(lambda_vol*) s, max_j abs(beta_ref_j)),
// the size of a term of the stationarity. NaN where a Newton system of the refinement is singular, 0 with no
// weight held.
double bound_multiplier_error(const po::Spec& spec, const po::Data& data, const Result& result)
{
    constexpr int kSteps = 6;
    std::vector<double> budget = po::testing::budget_of(spec, data.asset_count());
    po::testing::WideKkt kkt = po::testing::wide_kkt_of(spec, result);
    for (int step = 0; step < kSteps; ++step)
    {
        if (!po::testing::wide_newton_step(spec, data, budget, kkt))
        {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    po::testing::WideEqualities equalities = po::testing::wide_equalities(spec, data, kkt.weights);
    std::vector<long double> reference = po::testing::wide_stationarity(data, budget, kkt, equalities);
    long double scale = std::fmax(std::fabs(kkt.multipliers.front()),
                                  std::fabs(kkt.multipliers.back()) * *spec.target_volatility);
    long double distance = 0.0L;
    for (std::size_t j = 0; j < result.weights.size(); ++j)
    {
        if (result.bound_multipliers[j] != 0.0)
        {
            scale = std::fmax(scale, std::fabs(reference[j]));
            distance = std::fmax(distance, std::fabs(result.bound_multipliers[j] - reference[j]));
        }
    }
    return static_cast<double>(distance / scale);
}

// The largest errors over the cases of one cond(Sigma), and those of the least-squares fit at the same
// weights (least_squares_fit, which leaves out the H dw term of beta), where the free weights determine it.
struct BoundTally
{
    int held = 0;
    double equalities = 0.0;
    double bounds = 0.0;
    double fitted_bounds = 0.0;
};

void tally_bound_case(const Shape& shape, const MultiplierCeiling& ceiling, BoundTally& tally)
{
    SCOPED_TRACE(po::testing::describe(shape));
    po::testing::GeneratedRiskParity generated = po::testing::generate(shape);
    const po::Spec& spec = generated.spec;
    Result result = po::solve(spec, generated.data).value();
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    po::testing::expect_risk_parity_kkt(spec, generated.data, result, band_of(shape.condition));
    po::testing::ForwardError forward = po::testing::forward_error(spec, generated.data, result);
    EXPECT_LE(forward.reference_over_floor, 1e-2);
    double bounds = bound_multiplier_error(spec, generated.data, result);
    EXPECT_LE(forward.multipliers, ceiling.equalities);
    EXPECT_LE(bounds, ceiling.bounds);
    tally.held += static_cast<int>(
        std::ranges::count_if(result.bound_multipliers, [](double beta) { return beta != 0.0; }));
    tally.equalities = std::fmax(tally.equalities, forward.multipliers);
    tally.bounds = std::fmax(tally.bounds, bounds);
    std::optional<Result> fit = po::testing::least_squares_fit(spec, generated.data, result);
    if (fit.has_value())
    {
        double fitted = bound_multiplier_error(spec, generated.data, *fit);
        tally.fitted_bounds = std::fmax(tally.fitted_bounds, fitted);
    }
}

// The same with a box, LongOnly or Capped (lower 0, upper 1.2 / n), both budgets, n = 5 and 8, 4 seeds: the
// held weights and their beta, sign and complementarity included (expect_risk_parity_kkt), and beta by value
// against the refined KKT point. Measured largest (equalities / beta) on the instances of plan 8: 1.4e-14 / 8.6e-15
// at 1e4, 2.1e-13 / 2.6e-12 at 1e6, 2.0e-12 / 1.6e-10 at 1e8; the beta of the least-squares fit at the same weights,
// without the H dw term, reach 5.7e-12, 1.1e-8 and 3.1e-7. The ceilings below are those of the earlier instances
// (the largest measured there, 1.9e-14 / 1.8e-14, 2.6e-12 / 3.5e-13 and 1.4e-10 / 1.5e-10, times 10): they still
// hold, with 1.4 times to spare for beta at 1e6.
constexpr MultiplierCeiling kBoxMultiplierCeilings[] = {
    {1e4, 1.9e-13, 1.8e-13}, {1e6, 2.6e-11, 3.6e-12}, {1e8, 1.4e-9, 1.5e-9}};

TEST(RiskParityFuzz, TargetVolatilityBoundMultipliersByCondition)
{
    if (!po::testing::kLongDoubleIsExtended)
    {
        GTEST_SKIP() << "long double is double here: the refinement has no extra precision";
    }
    for (MultiplierCeiling ceiling : kBoxMultiplierCeilings)
    {
        BoundTally tally;
        for (std::size_t k = 0; k < 2 * 2 * 2 * 4; ++k)
        {
            Shape::Box box = k / 16 == 0 ? Shape::Box::LongOnly : Shape::Box::Capped;
            Shape::Budget budget = k / 8 % 2 == 0 ? Shape::Budget::Equal : Shape::Budget::Uneven;
            std::size_t n = k / 4 % 2 == 0 ? 5 : 8;
            tally_bound_case(Shape{n, ceiling.condition, 1.0, budget, box, false, true, 402 + k % 4}, ceiling,
                             tally);
        }
        std::string exponent = std::to_string(static_cast<int>(std::log10(ceiling.condition)));
        RecordProperty("held_" + exponent, tally.held);
        RecordProperty("equalities_" + exponent, ::testing::PrintToString(tally.equalities));
        RecordProperty("bounds_" + exponent, ::testing::PrintToString(tally.bounds));
        RecordProperty("fitted_bounds_" + exponent, ::testing::PrintToString(tally.fitted_bounds));
        EXPECT_GE(tally.held, 20);
        EXPECT_GT(tally.fitted_bounds, 10.0 * ceiling.bounds);
    }
}

void print_quantiles(const char* name, std::vector<double> values)
{
    std::ranges::sort(values);
    if (values.empty())
    {
        return;
    }
    auto at = [&](double q)
    { return values[static_cast<std::size_t>(q * static_cast<double>(values.size() - 1))]; };
    std::printf("%s: count %zu, median %.2g, 90%% %.2g, 99%% %.2g, max %.2g\n", name, values.size(), at(0.5),
                at(0.9), at(0.99), values.back());
}

// The distribution behind docs/guide/validation.md: the weight error of the thinnest class over many seeds,
// n = 20 (200 seeds) and n = 50 (30 seeds, cut at 240 s of solves).
TEST(RiskParityFuzz, DISABLED_MeasureWeightErrorDistribution)
{
    for (auto [n, seeds] : {std::pair<std::size_t, std::uint64_t>{20, 200}, {50, 30}})
    {
        std::vector<double> errors;
        std::vector<double> floor_limited;
        std::vector<double> over_bound;
        std::vector<double> floor_bound;
        std::vector<double> values;
        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
        for (std::uint64_t seed = 1; seed <= seeds; ++seed)
        {
            if (std::chrono::steady_clock::now() - start > std::chrono::seconds(240))
            {
                std::printf("n %zu cut after %llu seeds\n", n, static_cast<unsigned long long>(seed - 1));
                break;
            }
            po::testing::GeneratedRiskParity generated = po::testing::generate(thinnest_class(n, seed));
            Result result = po::solve(generated.spec, generated.data).value();
            if (result.status != SolveStatus::Optimal)
            {
                std::printf("seed %llu status %d\n", static_cast<unsigned long long>(seed),
                            static_cast<int>(result.status));
                continue;
            }
            WeightCheck check = weight_check(generated, result, 1e-7);
            errors.push_back(check.forward.weights);
            if (check.floor_limited)
            {
                floor_limited.push_back(check.forward.weights);
            }
            over_bound.push_back(check.forward.weights / check.bound.residual);
            floor_bound.push_back(check.bound.floor);
            values.push_back(check.forward.value / (1e-7 * (1.0 + check.value)));
        }
        std::printf("n %zu\n", n);
        print_quantiles("weight error, all OPTIMAL", errors);
        std::printf("above 1e-7: %td\n",
                    std::ranges::count_if(errors, [](double error) { return error > 1e-7; }));
        print_quantiles("weight error, above the band", floor_limited);
        print_quantiles("weight error over its residual bound", over_bound);
        print_quantiles("bound from the floor", floor_bound);
        print_quantiles("value error over band (1 + f)", values);
    }
}

} // namespace
