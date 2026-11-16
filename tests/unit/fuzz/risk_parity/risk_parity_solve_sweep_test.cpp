#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "fuzz/risk_parity/risk_parity_reference.hpp"
#include "fuzz/risk_parity/risk_parity_refinement.hpp"
#include "problems/constraints.hpp"
#include "qp/goldfarb_idnani.hpp"
#include "qp/problem.hpp"
#include "riskparity/augmented_lagrangian.hpp"
#include "riskparity/kkt.hpp"
#include "riskparity/newton_polish.hpp"
#include "riskparity/spinu.hpp"
#include "solve/solve.hpp"

namespace
{

using po::Result;
using po::SolveStatus;
using Shape = po::testing::RiskParityShape;

// OPTIMAL results of OptimalRunsAreKkt out of its 576 shapes, as measured (docs/guide/validation.md). The
// solver is deterministic: fewer is a regression.
constexpr int kOptimalFloor = 564;
// The reference of a floor-limited result is at most this fraction of the rounding floor of the result from
// stationarity: the point the forward error is measured against is that much closer to a KKT point.
constexpr double kReferenceOverFloor = 1e-2;

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

struct SolveRun
{
    Shape shape;
    // False when the linear rows of the shape have no feasible point and the shape was not generated.
    bool generated = false;
    Result result;
    // The failures of expect_risk_parity_kkt at the band of cond(Sigma); only checked for an OPTIMAL result.
    std::vector<std::string> kkt_failures;
    // Of an OPTIMAL result with finite multipliers: its stationarity by the two measures of the check, and
    // whether it is above the band, so that the check passes it by the rounding floor alone.
    po::testing::RiskParityStationarity stationarity{0.0, 0.0};
    bool floor_limited = false;
    // Of a floor-limited result, against the KKT point of its active set refined in long double.
    po::testing::ForwardError forward{0.0, 0.0, 0.0, 0.0};
    double value = 0.0;
    bool no_point_is_nan = true;
    double seconds = 0.0;
};

bool is_all_nan(const Result& result)
{
    bool all = std::isnan(result.budget_multiplier) && std::isnan(result.target_return_multiplier) &&
               std::isnan(result.volatility_multiplier) && std::isnan(result.residual);
    auto is_nan = [](double value) { return std::isnan(value); };
    return all && std::ranges::all_of(result.weights, is_nan) &&
           std::ranges::all_of(result.bound_multipliers, is_nan);
}

void check_optimal(const po::testing::GeneratedRiskParity& generated, SolveRun& run)
{
    const po::Spec& spec = generated.spec;
    double band = band_of(run.shape.condition);
    run.kkt_failures = po::testing::risk_parity_kkt_failures(spec, generated.data, run.result, band);
    if (!std::isfinite(run.result.budget_multiplier))
    {
        return;
    }
    run.stationarity = po::testing::risk_parity_stationarity(spec, generated.data, run.result, band);
    run.floor_limited = !(run.stationarity.relative <= band);
    run.value = po::testing::risk_parity_value(generated.data, spec.risk_budget, run.result.weights);
    if (run.floor_limited && po::testing::kLongDoubleIsExtended)
    {
        run.forward = po::testing::forward_error(spec, generated.data, run.result);
    }
}

SolveRun run_shape(const Shape& shape)
{
    SolveRun run;
    run.shape = shape;
    try
    {
        po::testing::GeneratedRiskParity generated = po::testing::generate(shape);
        run.generated = true;
        std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
        run.result = po::solve(generated.spec, generated.data).value();
        run.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        if (run.result.status != SolveStatus::Optimal)
        {
            run.no_point_is_nan = is_all_nan(run.result);
            return run;
        }
        check_optimal(generated, run);
    }
    catch (const std::runtime_error&)
    {
        run.generated = false;
    }
    return run;
}

// Every combination of cond(Sigma) in {1e2, 1e4, 1e6, 1e8}, scale {1, 1e-4}, equal and uneven budgets, the
// three boxes, target_return and target_volatility (192 shapes), for each n and seed.
std::vector<SolveRun> sweep(const std::vector<std::size_t>& sizes, std::uint64_t seeds)
{
    std::vector<SolveRun> runs;
    for (std::size_t n : sizes)
    {
        for (std::uint64_t seed = 1; seed <= seeds; ++seed)
        {
            for (const Shape& base : po::testing::risk_parity_shapes(192, 1))
            {
                Shape shape = base;
                shape.n = n;
                shape.seed = seed * 7919 + n;
                runs.push_back(run_shape(shape));
            }
        }
    }
    return runs;
}

// |f(w) - f(w_ref)| over band (1 + f): at most 1 when the value of the result is right at the band.
double value_error_over_band(const SolveRun& run)
{
    return run.forward.value / (band_of(run.shape.condition) * (1.0 + run.value));
}

struct Tally
{
    int runs = 0;
    int optimal = 0;
    int infeasible = 0;
    int max_iterations = 0;
    int kkt_failures = 0;
    int floor_limited = 0;
    long iterations = 0;
    int most_iterations = 0;
    double largest_over_floor = 0.0;
    double largest_weight_error = 0.0;
    double largest_value_error = 0.0;
    double largest_reference = 0.0;
    double seconds = 0.0;
    double longest = 0.0;
};

void count(const SolveRun& run, Tally& tally)
{
    if (!run.generated)
    {
        return;
    }
    ++tally.runs;
    tally.seconds += run.seconds;
    tally.longest = std::fmax(tally.longest, run.seconds);
    tally.infeasible += run.result.status == SolveStatus::Infeasible;
    tally.max_iterations += run.result.status == SolveStatus::MaxIterations;
    if (run.result.status != SolveStatus::Optimal)
    {
        return;
    }
    ++tally.optimal;
    tally.kkt_failures += !run.kkt_failures.empty();
    tally.iterations += run.result.iterations;
    tally.most_iterations = std::max(tally.most_iterations, run.result.iterations);
    if (!run.floor_limited)
    {
        return;
    }
    ++tally.floor_limited;
    tally.largest_over_floor = po::testing::worst_of(tally.largest_over_floor, run.stationarity.over_floor);
    tally.largest_weight_error = po::testing::worst_of(tally.largest_weight_error, run.forward.weights);
    tally.largest_value_error = po::testing::worst_of(tally.largest_value_error, value_error_over_band(run));
    tally.largest_reference =
        po::testing::worst_of(tally.largest_reference, run.forward.reference_over_floor);
}

void print_table(const char* first_column, const std::map<std::string, Tally>& table)
{
    std::printf("| %s | optimal | max iter | infeasible | KKT failed | floor-limited | stationarity / floor "
                "| weight error | value error / band | reference / floor | iterations | seconds (longest) "
                "|\n",
                first_column);
    for (const auto& [key, tally] : table)
    {
        double mean = 0.0;
        if (tally.optimal > 0)
        {
            mean = static_cast<double>(tally.iterations) / tally.optimal;
        }
        std::printf("| %s | %d / %d | %d | %d | %d | %d | %.2f | %.2g | %.2g | %.2g | %.1f (%d) | "
                    "%.2f (%.3f) |\n",
                    key.c_str(), tally.optimal, tally.runs, tally.max_iterations, tally.infeasible,
                    tally.kkt_failures, tally.floor_limited, tally.largest_over_floor,
                    tally.largest_weight_error, tally.largest_value_error, tally.largest_reference, mean,
                    tally.most_iterations, tally.seconds, tally.longest);
    }
}

const char* box_name(const Shape& shape)
{
    constexpr const char* kBoxes[] = {"None", "LongOnly", "Capped"};
    return kBoxes[static_cast<int>(shape.box)];
}

std::string constraints_name(const Shape& shape)
{
    std::string name = box_name(shape);
    if (shape.target_return)
    {
        name += " + return";
    }
    if (shape.target_volatility)
    {
        name += " + volatility";
    }
    return name;
}

void print_failures(const SolveRun& run)
{
    if (run.kkt_failures.empty())
    {
        return;
    }
    std::printf("KKT failure: %s\n", po::testing::describe(run.shape).GetString().c_str());
    for (const std::string& failure : run.kkt_failures)
    {
        std::printf("%s\n", failure.c_str());
    }
}

void print_sweep(const std::vector<SolveRun>& runs)
{
    std::map<std::string, Tally> by_constraints;
    std::map<std::string, Tally> by_condition;
    std::map<std::string, Tally> by_size;
    std::map<std::string, Tally> floor_limited;
    for (const SolveRun& run : runs)
    {
        char key[96];
        count(run, by_constraints[constraints_name(run.shape)]);
        std::snprintf(key, sizeof(key), "%g", run.shape.condition);
        count(run, by_condition[key]);
        std::snprintf(key, sizeof(key), "%3zu", run.shape.n);
        count(run, by_size[key]);
        print_failures(run);
        if (run.floor_limited)
        {
            const char* box = box_name(run.shape);
            std::snprintf(key, sizeof(key), "%g | %3zu | %s", run.shape.condition, run.shape.n, box);
            count(run, floor_limited[key]);
        }
    }
    print_table("constraints", by_constraints);
    print_table("cond", by_condition);
    print_table("n", by_size);
    print_table("cond | n | box, the floor-limited results", floor_limited);
}

// A result that the check passes by the rounding floor of its stationarity alone must still be the right
// point: against the KKT point of its active set refined in long double, the weights and f agree at the band.
void expect_forward_error_at_the_band(const SolveRun& run)
{
    if (!po::testing::kLongDoubleIsExtended)
    {
        return;
    }
    EXPECT_LE(run.forward.reference_over_floor, kReferenceOverFloor);
    EXPECT_LE(run.forward.reference_gap, band_of(run.shape.condition));
    EXPECT_LE(run.forward.weights, band_of(run.shape.condition));
    EXPECT_LE(value_error_over_band(run), 1.0);
}

// No result is silently wrong. Every status other than OPTIMAL comes without a point. Every OPTIMAL passes
// expect_risk_parity_kkt at the band of its cond(Sigma), and those it passes by the rounding floor of the
// stationarity are also checked by their forward error.
TEST(RiskParitySolveSweep, OptimalRunsAreKkt)
{
    std::vector<SolveRun> runs = sweep({3, 5, 10}, 1);
    int optimal = 0;
    int floor_limited = 0;
    for (const SolveRun& run : runs)
    {
        SCOPED_TRACE(po::testing::describe(run.shape));
        if (!run.generated)
        {
            continue;
        }
        EXPECT_TRUE(run.no_point_is_nan);
        EXPECT_NE(run.result.status, SolveStatus::Unbounded);
        optimal += run.result.status == SolveStatus::Optimal;
        floor_limited += run.floor_limited;
        for (const std::string& failure : run.kkt_failures)
        {
            ADD_FAILURE() << failure;
        }
        if (run.floor_limited)
        {
            expect_forward_error_at_the_band(run);
        }
    }
    EXPECT_EQ(runs.size(), 576u);
    EXPECT_GE(optimal, kOptimalFloor);
    // The grid has results above the band, or the forward-error check above checked nothing.
    EXPECT_GT(floor_limited, 0);
}

// The measurement of docs/guide/validation.md: run with --gtest_also_run_disabled_tests.
TEST(RiskParitySolveSweep, DISABLED_MeasureSweep)
{
    print_sweep(sweep({3, 5, 10, 20, 50}, 3));
}

// The start of the solver: the point of the polyhedron nearest to the Spinu solution.
std::vector<double> solver_start(const po::testing::GeneratedRiskParity& generated)
{
    std::size_t n = generated.data.asset_count();
    std::vector<double> spinu = po::riskparity::spinu(generated.data, generated.spec.risk_budget).weights;
    po::qp::Problem problem;
    problem.n = n;
    problem.hessian.assign(n * n, 0.0);
    problem.linear.resize(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        problem.hessian[i * n + i] = 1.0;
        problem.linear[i] = -spinu[i];
    }
    po::problems::add_constraint_rows(generated.spec, generated.data, problem);
    return po::qp::solve(problem).x;
}

struct AugmentedTally
{
    int runs = 0;
    int converged = 0;
    int polished = 0;
    long outer = 0;
    int most_outer = 0;
    long inner = 0;
    int most_inner = 0;
    int most_steps = 0;
    double largest_penalty = 0.0;
    double largest_multiplier = 0.0;
};

void count_augmented(const Shape& shape, const po::riskparity::AugmentedOptions& options,
                     AugmentedTally& tally)
{
    po::testing::GeneratedRiskParity generated = po::testing::generate(shape);
    po::riskparity::KktProblem problem{generated.spec, generated.data, generated.spec.risk_budget};
    po::riskparity::AugmentedResult result =
        po::riskparity::augmented_lagrangian(problem, solver_start(generated), options);
    po::riskparity::PolishResult polish = po::riskparity::newton_polish(problem, result.point);
    ++tally.runs;
    tally.polished += po::riskparity::is_kkt_point(problem, polish.point);
    tally.most_steps = std::max(tally.most_steps, polish.steps);
    if (!result.converged)
    {
        return;
    }
    ++tally.converged;
    tally.outer += result.outer_iterations;
    tally.most_outer = std::max(tally.most_outer, result.outer_iterations);
    tally.inner += result.inner_iterations;
    tally.most_inner = std::max(tally.most_inner, result.inner_iterations);
    tally.largest_penalty = std::fmax(tally.largest_penalty, result.penalty);
    double multiplier = result.point.equality_multipliers.back() * *generated.spec.target_volatility;
    tally.largest_multiplier = std::fmax(tally.largest_multiplier, std::fabs(multiplier));
}

// The augmented Lagrangian alone on the shapes with a target volatility, at the two efforts of the solver,
// each followed by the polish: the numbers behind its constants in docs/guide/validation.md.
std::map<std::string, AugmentedTally> augmented_sweep(const po::riskparity::AugmentedOptions& options)
{
    std::map<std::string, AugmentedTally> table;
    for (std::size_t n : {3, 5, 10, 20})
    {
        for (const Shape& base : po::testing::risk_parity_shapes(192, 1))
        {
            Shape shape = base;
            shape.n = n;
            shape.seed = 7919 + n;
            if (shape.target_volatility)
            {
                count_augmented(shape, options, table[constraints_name(shape)]);
            }
        }
    }
    return table;
}

TEST(RiskParitySolveSweep, DISABLED_MeasureAugmentedLagrangian)
{
    for (po::riskparity::AugmentedOptions options :
         {po::riskparity::AugmentedOptions{1e-6, 50, 200}, po::riskparity::AugmentedOptions{1e-11, 50, 2000}})
    {
        std::printf("tolerance %g, at most %d iterations per descent\n", options.tolerance,
                    options.max_inner_iterations);
        std::printf("| constraints | converged | KKT after the polish | outer | descent iterations | "
                    "polish steps, most | largest rho s^2 | largest abs(lambda s) |\n");
        for (const auto& [key, tally] : augmented_sweep(options))
        {
            double runs = std::max(tally.converged, 1);
            std::printf("| %s | %d / %d | %d | %.1f (%d) | %.0f (%d) | %d | %.0e | %.3g |\n", key.c_str(),
                        tally.converged, tally.runs, tally.polished, tally.outer / runs, tally.most_outer,
                        tally.inner / runs, tally.most_inner, tally.most_steps, tally.largest_penalty,
                        tally.largest_multiplier);
        }
    }
}

} // namespace
