#include <gtest/gtest-spi.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "fuzz/frontier_reference.hpp"
#include "fuzz/risk_parity/risk_parity_reference.hpp"
#include "riskparity/kkt.hpp"
#include "riskparity/newton_polish.hpp"
#include "riskparity/polyhedral_descent.hpp"

namespace
{

using po::riskparity::DescentOptions;
using po::riskparity::DescentProblem;
using po::riskparity::DescentResult;
using po::riskparity::DescentStop;
using Shape = po::testing::RiskParityShape;

constexpr int kMaxIterations = 2000;
// Converged runs of ConvergedRunsAreKkt out of 640, as measured; the descent is deterministic, so fewer is a
// regression.
constexpr int kConvergedFloor = 625;
constexpr double kToleranceFraction = 0.1;

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

struct DescentRun
{
    Shape shape;
    DescentResult result;
    // The first failure of expect_risk_parity_kkt at the band of cond(Sigma), empty when it passed; only
    // checked for a converged run.
    std::string kkt_failure;
};

po::Result as_result(const po::Spec& spec, const DescentResult& descent)
{
    po::Result result;
    result.status = po::SolveStatus::Optimal;
    result.weights = descent.weights;
    result.bound_multipliers = descent.bound_multipliers;
    result.budget_multiplier = descent.equality_multipliers.at(0);
    if (spec.target_return.has_value())
    {
        result.target_return_multiplier = descent.equality_multipliers.at(1);
    }
    return result;
}

std::string kkt_failure_of(const po::testing::GeneratedRiskParity& generated, const DescentResult& result,
                           double band)
{
    ::testing::TestPartResultArray failures;
    {
        ::testing::ScopedFakeTestPartResultReporter reporter(
            ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD, &failures);
        po::testing::expect_risk_parity_kkt(generated.spec, generated.data, as_result(generated.spec, result),
                                            band);
    }
    std::string failure;
    if (failures.size() > 0)
    {
        failure = failures.GetTestPartResult(0).message();
    }
    return failure;
}

// One descent from the minimum-variance point of the polyhedron to a tenth of the band of cond(Sigma): the
// KKT residual of the descent (double, multipliers of the QP) and that of expect_risk_parity_kkt (long
// double) agree to a few percent, not to the last digit, so a stop at the band itself can land just above it.
DescentRun run_shape(const Shape& shape, double proximal)
{
    po::testing::GeneratedRiskParity generated = po::testing::generate(shape);
    const po::Spec& spec = generated.spec;
    po::Result minimum = po::testing::frontier_point(spec, generated.data, spec.target_return);
    EXPECT_EQ(minimum.status, po::SolveStatus::Optimal) << po::testing::describe(shape);
    double band = band_of(shape.condition);
    DescentProblem problem{spec, generated.data, spec.risk_budget, std::nullopt};
    DescentOptions options{proximal, kToleranceFraction * band, kMaxIterations};
    DescentRun run{shape, po::riskparity::descend(problem, minimum.weights, options), {}};
    if (run.result.stop == DescentStop::Converged)
    {
        run.kkt_failure = kkt_failure_of(generated, run.result, band);
    }
    return run;
}

// Every combination of cond(Sigma) in {1e2, 1e4, 1e6, 1e8}, scale {1, 1e-4}, equal and uneven budgets, the
// three boxes and target_return, except Box None with target_return (an unbounded line, where f can decrease
// towards its value at infinity), for each n and seed.
std::vector<DescentRun> sweep(const std::vector<std::size_t>& sizes, std::uint64_t seeds, double proximal)
{
    std::vector<DescentRun> runs;
    for (std::size_t n : sizes)
    {
        for (std::uint64_t seed = 1; seed <= seeds; ++seed)
        {
            for (const Shape& base : po::testing::risk_parity_shapes(96, 1))
            {
                Shape shape = base;
                if (shape.box == Shape::Box::None && shape.target_return)
                {
                    continue;
                }
                shape.n = n;
                shape.seed = seed * 7919 + n;
                runs.push_back(run_shape(shape, proximal));
            }
        }
    }
    return runs;
}

struct Tally
{
    int runs = 0;
    int converged = 0;
    int capped = 0;
    int no_decrease = 0;
    int failed = 0;
    int kkt_failures = 0;
    long iterations = 0;
    int most_iterations = 0;
};

void count(const DescentRun& run, Tally& tally)
{
    ++tally.runs;
    switch (run.result.stop)
    {
    case DescentStop::Converged:
        ++tally.converged;
        tally.iterations += run.result.iterations;
        tally.most_iterations = std::max(tally.most_iterations, run.result.iterations);
        break;
    case DescentStop::IterationCap:
        ++tally.capped;
        break;
    case DescentStop::NoDecrease:
        ++tally.no_decrease;
        break;
    case DescentStop::NotPositiveDefinite:
    case DescentStop::QpFailed:
        ++tally.failed;
        break;
    }
    if (!run.kkt_failure.empty())
    {
        ++tally.kkt_failures;
    }
}

void print_table(const char* title, const std::map<std::string, Tally>& table)
{
    std::printf("%s\n", title);
    for (const auto& [key, tally] : table)
    {
        double mean = 0.0;
        if (tally.converged > 0)
        {
            mean = static_cast<double>(tally.iterations) / tally.converged;
        }
        std::printf("| %s | %d / %d | %d | %d | %d | %d | %.1f (%d) |\n", key.c_str(), tally.converged,
                    tally.runs, tally.capped, tally.no_decrease, tally.failed, tally.kkt_failures, mean,
                    tally.most_iterations);
    }
}

std::string box_name(Shape::Box box)
{
    switch (box)
    {
    case Shape::Box::None:
        return "None";
    case Shape::Box::LongOnly:
        return "LongOnly";
    case Shape::Box::Capped:
        return "Capped";
    }
    return "";
}

void print_sweep(const std::vector<DescentRun>& runs)
{
    std::map<std::string, Tally> by_condition;
    std::map<std::string, Tally> by_size;
    for (const DescentRun& run : runs)
    {
        char key[96];
        std::snprintf(key, sizeof(key), "%g | %g | %s", run.shape.condition, run.shape.scale,
                      box_name(run.shape.box).c_str());
        count(run, by_condition[key]);
        std::snprintf(key, sizeof(key), "%3zu", run.shape.n);
        count(run, by_size[key]);
    }
    print_table("| cond | scale | box | converged | cap | no decrease | failed | KKT failed | iterations |",
                by_condition);
    print_table("| n | converged | cap | no decrease | failed | KKT failed | iterations |", by_size);
}

// Every converged run is a KKT point at the band of its cond(Sigma), no run fails to factor Q or to solve the
// QP, and the share of converged runs stays at the measured level (docs/guide/validation.md).
TEST(RiskParityDescentFuzz, ConvergedRunsAreKkt)
{
    std::vector<DescentRun> runs = sweep({3, 5, 10, 20}, 2, 0.05);
    int converged = 0;
    for (const DescentRun& run : runs)
    {
        SCOPED_TRACE(po::testing::describe(run.shape));
        EXPECT_TRUE(run.kkt_failure.empty()) << run.kkt_failure;
        EXPECT_NE(run.result.stop, DescentStop::NotPositiveDefinite);
        EXPECT_NE(run.result.stop, DescentStop::QpFailed);
        if (run.result.stop == DescentStop::Converged)
        {
            ++converged;
        }
    }
    EXPECT_EQ(runs.size(), 640u);
    EXPECT_GE(converged, kConvergedFloor);
}

// At a converged LongOnly point with weights at 0 the identity w'grad f = 0 holds at the band with the
// reference gradient, and fails when one of its three terms is off by 0.1%: its scale is not so loose that it
// accepts a wrong gradient.
TEST(RiskParityDescentFuzz, HomogeneityIdentityCatchesAWrongGradient)
{
    Shape shape{10, 1e6, 1.0, Shape::Budget::Uneven, Shape::Box::LongOnly, false, false, 3 * 7919 + 10};
    DescentRun run = run_shape(shape, 0.05);
    ASSERT_EQ(run.result.stop, DescentStop::Converged);
    EXPECT_TRUE(run.kkt_failure.empty()) << run.kkt_failure;
    const std::vector<double>& weights = run.result.weights;
    ASSERT_TRUE(std::ranges::any_of(weights, [](double weight) { return std::fabs(weight) <= 1e-12; }));
    po::testing::GeneratedRiskParity generated = po::testing::generate(shape);
    std::vector<std::vector<double>> terms =
        po::testing::risk_parity_gradient_terms(generated.data, generated.spec.risk_budget, weights);
    double band = band_of(shape.condition);
    po::testing::expect_homogeneity_identity(weights, po::testing::sum_of(terms), band);
    for (double& entry : terms[1])
    {
        entry *= 1.001;
    }
    po::testing::ObjectiveGradient wrong = po::testing::sum_of(terms);
    EXPECT_NONFATAL_FAILURE(po::testing::expect_homogeneity_identity(weights, wrong, band), "tolerance");
}

// The measurement of docs/guide/validation.md: run with --gtest_also_run_disabled_tests.
TEST(RiskParityDescentFuzz, DISABLED_MeasureSweep)
{
    for (double proximal : {0.05, 0.005})
    {
        std::printf("proximal %g\n", proximal);
        print_sweep(sweep({3, 5, 10, 20, 50}, 5, proximal));
    }
}

// Whether the Newton polish of the solver turns the descents that did not converge into KKT points, as
// is_kkt_point decides: the count of docs/guide/validation.md.
TEST(RiskParityDescentFuzz, DISABLED_MeasurePolishOfUnconverged)
{
    int unconverged = 0;
    int repaired = 0;
    for (const DescentRun& run : sweep({3, 5, 10, 20, 50}, 5, po::riskparity::kDescentProximal))
    {
        if (run.result.stop == DescentStop::Converged)
        {
            continue;
        }
        ++unconverged;
        po::testing::GeneratedRiskParity generated = po::testing::generate(run.shape);
        const po::Spec& spec = generated.spec;
        po::riskparity::KktProblem problem{spec, generated.data, spec.risk_budget};
        po::riskparity::KktPoint point{run.result.weights, run.result.equality_multipliers,
                                       run.result.bound_multipliers};
        po::riskparity::PolishResult polished = po::riskparity::newton_polish(problem, point);
        bool kkt = po::riskparity::is_kkt_point(problem, polished.point);
        repaired += kkt;
        std::printf("%s | stop %d | polish steps %d | KKT %d\n",
                    po::testing::describe(run.shape).GetString().c_str(),
                    static_cast<int>(run.result.stop), polished.steps, static_cast<int>(kkt));
    }
    std::printf("unconverged %d, KKT after the polish %d\n", unconverged, repaired);
}

} // namespace
