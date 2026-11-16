#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "core/result.hpp"
#include "fuzz/random_problems.hpp"
#include "fuzz/ratio_reference.hpp"
#include "solve/solve.hpp"
#include "user_kkt.hpp"

namespace
{

using po::Data;
using po::Objective;
using po::Result;
using po::SolveStatus;
using po::Spec;
using po::testing::ProblemShape;
using po::testing::RatioReference;
using po::testing::describe;

constexpr std::array kObjectives = {Objective::MaxSharpe, Objective::MaxDiversification};

struct Counts
{
    std::size_t optimal = 0;
    std::size_t infeasible = 0;
    std::size_t unbounded = 0;
    std::size_t binding = 0;
};

// No point: every weight, bound multiplier, the budget and leverage multipliers and the residual are NaN.
void expect_no_point(const Result& result)
{
    for (double weight : result.weights)
    {
        EXPECT_TRUE(std::isnan(weight));
    }
    for (double multiplier : result.bound_multipliers)
    {
        EXPECT_TRUE(std::isnan(multiplier));
    }
    EXPECT_TRUE(std::isnan(result.budget_multiplier));
    EXPECT_TRUE(std::isnan(result.leverage_multiplier));
    EXPECT_TRUE(std::isnan(result.residual));
}

// The status of the reference; when OPTIMAL, the ratio of the reference within the band, the cash weight
// exactly 0 and KKT in the user's variables.
void check_case(const Spec& spec, const Data& data, double tolerance, Counts& counts)
{
    std::vector<double> v = po::testing::ratio_vector(data, spec.objective);
    RatioReference reference = po::testing::max_ratio_reference(spec, data, v);
    Result result = po::solve(spec, data).value();
    EXPECT_EQ(result.status, reference.status);
    if (result.status != SolveStatus::Optimal || reference.status != SolveStatus::Optimal)
    {
        counts.infeasible += result.status == SolveStatus::Infeasible ? 1 : 0;
        counts.unbounded += result.status == SolveStatus::Unbounded ? 1 : 0;
        if (result.status != SolveStatus::Optimal)
        {
            expect_no_point(result);
        }
        return;
    }
    ++counts.optimal;
    counts.binding += result.leverage_multiplier < 0.0 ? 1 : 0;
    if (data.risk_free().has_value())
    {
        EXPECT_EQ(result.weights.back(), 0.0);
    }
    EXPECT_NEAR(po::testing::ratio_of(data, v, result.weights), reference.ratio, tolerance * reference.ratio);
    po::testing::expect_user_kkt(po::testing::with_cash_fixed(spec, data), data, result,
                                 po::testing::ratio_gradient(data, v, result), tolerance);
}

void check_shapes(const std::vector<ProblemShape>& shapes, Counts& counts)
{
    for (const ProblemShape& shape : shapes)
    {
        SCOPED_TRACE(describe(shape));
        po::testing::GeneratedProblem problem = po::testing::generate(shape);
        double tolerance = po::testing::condition_tolerance(shape.condition);
        for (Objective objective : kObjectives)
        {
            SCOPED_TRACE(static_cast<int>(objective));
            Spec spec = problem.spec;
            spec.objective = objective;
            check_case(spec, problem.data, tolerance, counts);
        }
    }
}

// MAX_SHARPE and MAX_DIVERSIFICATION over fuzz_shapes(300, 3), without and with the leverage of each shape,
// against the frontier reference: the same status, and when OPTIMAL the same ratio and KKT, all within the
// band. Of the 1200 cases, 1039 are OPTIMAL (91 binding the leverage), 133 INFEASIBLE and 28 UNBOUNDED
// (docs/guide/validation.md).
TEST(RatioFuzz, AgainstFrontierReference)
{
    std::vector<ProblemShape> shapes = po::testing::fuzz_shapes(300, 3);
    Counts counts;
    check_shapes(shapes, counts);
    check_shapes(po::testing::with_leverage(shapes), counts);
    EXPECT_GE(counts.optimal, 1000u);
    EXPECT_GE(counts.infeasible, 100u);
    EXPECT_GE(counts.unbounded, 20u);
    EXPECT_GE(counts.binding, 80u);
}

} // namespace
