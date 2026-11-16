#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <chrono>
#include <cstdio>
#include <limits>
#include <numeric>
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

// Without a point: MAX_ITER or INFEASIBLE, every number NaN. With one: KKT at the band of cond(Sigma).
void expect_never_silently_wrong(const po::Spec& spec, const po::Data& data, const Result& result,
                                 double condition)
{
    if (result.status != SolveStatus::Optimal)
    {
        EXPECT_TRUE(result.status == SolveStatus::MaxIterations || result.status == SolveStatus::Infeasible);
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
constexpr int kOtherLocalMinima = 6;

// The scan's value moved to the target to first order: its point sits on the curve only to the rounding of
// the curve (a gap of up to 2e-10 of s measured), and f changes by lambda_k h_k off the curve, with the
// multipliers of the result. Without this a scan point just inside the ellipse beats an exact solution.
double scan_value_at_target(const po::Spec& spec, const po::Data& data, const Result& result,
                            const po::testing::ScanResult& scan)
{
    double total = std::accumulate(scan.weights.begin(), scan.weights.end(), 0.0);
    double value = scan.value - result.budget_multiplier * (total - 1.0);
    if (spec.target_return.has_value())
    {
        double mean = po::testing::compensated_dot_of<double>(3, data.mu(), scan.weights);
        value -= result.target_return_multiplier * (mean - *spec.target_return);
    }
    if (spec.target_volatility.has_value())
    {
        double volatility = po::testing::wide_volatility(data, scan.weights);
        value -= result.volatility_multiplier * (volatility - *spec.target_volatility);
    }
    return value;
}

// No feasible point of the curve within 20 scan steps of the result has a smaller f: the result is a local
// minimum along the curve, not a saddle or a maximum of it.
bool is_local_minimum_on_curve(const po::Spec& spec, const po::Data& data, const Result& result, double value)
{
    constexpr int kPoints = 20001;
    constexpr int kWindow = 20;
    po::testing::ScanCurve curve = po::testing::scan_curve_of(spec, data).value();
    double step = (curve.last - curve.first) / (kPoints - 1);
    int nearest = 0;
    double nearest_distance = std::numeric_limits<double>::infinity();
    for (int k = 0; k < kPoints; ++k)
    {
        double distance = largest_distance(curve.point(curve.first + step * k), result.weights);
        nearest = distance < nearest_distance ? k : nearest;
        nearest_distance = std::fmin(distance, nearest_distance);
    }
    for (int k = std::max(0, nearest - kWindow); k <= std::min(kPoints - 1, nearest + kWindow); ++k)
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
    if (value <= scan_value_at_target(spec, generated.data, result, scan) * (1.0 + 1e-9))
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
        budget[i] = std::pow(10.0, -6.0 + 6.0 * static_cast<double>(i) / static_cast<double>(n - 1));
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

// The descents of the polyhedral-descent sweep that the polish does not turn into a KKT point, through
// po::solve: OPTIMAL and KKT at the band, or no point (docs/guide/validation.md).
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

// The weights of a result at its rounding floor are determined only to how much the KKT matrix of its active
// set amplifies that floor: the bound of weight_error_bound, which this factor covers with room for the
// second-order terms it leaves out (the largest error over the bound measured: docs/guide/validation.md).
constexpr double kBoundSafety = 4.0;
// The weight error of this class is not held to the band of cond(Sigma), which does not apply at the floor,
// but it must not grow: the largest of 229 measured results is 3.7e-5 (docs/guide/validation.md).
constexpr double kWeightErrorCeiling = 1e-4;

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
// Seed 4 (31696) is a result whose weights are 2.3e-5 from the refined KKT point while its stationarity is
// at 0.22 of its rounding floor. The value is right at the band and the weights within their a posteriori
// bound; the band of the weights does not apply to a result at its floor.
TEST(RiskParityFuzz, WeightsAtCondition1e8WithinTheirAPosterioriBound)
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
        EXPECT_LE(check.forward.weights, kBoundSafety * check.bound.residual)
            << "bound " << check.bound.residual;
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
