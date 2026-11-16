#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "core/result.hpp"
#include "fuzz/frontier_reference.hpp"
#include "fuzz/random_problems.hpp"
#include "problems/solve.hpp"
#include "random_covariance.hpp"
#include "user_kkt.hpp"

namespace
{

using po::Data;
using po::Result;
using po::SolveStatus;
using po::Spec;
using po::testing::GeneratedProblem;
using po::testing::ProblemShape;

// KKT residuals and agreement with the reference, relative, by cond(Sigma) of the covariance solved; above 1e8,
// the tolerance of 1e8.
double condition_tolerance(double condition)
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

testing::Message describe(const ProblemShape& shape)
{
    return testing::Message() << "n " << shape.n << " cond " << shape.condition << " scale " << shape.scale
                              << " means " << static_cast<int>(shape.means) << " box "
                              << static_cast<int>(shape.box) << " cash " << shape.cash << " seed " << shape.seed;
}

Spec max_return(Spec bounds, double cap)
{
    bounds.objective = po::Objective::MaxReturn;
    bounds.max_volatility = cap;
    return bounds;
}

double largest_asset_volatility(const Data& data)
{
    std::size_t n = data.asset_count();
    double largest = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        largest = std::max(largest, std::sqrt(data.covariance()[i * n + i]));
    }
    return largest;
}

double minimum_volatility(const Spec& bounds, const Data& data)
{
    return po::testing::portfolio_volatility(data, po::testing::frontier_point(bounds, data, std::nullopt).weights);
}

// f = -(mu - rf 1)'w, written here independently of the builder.
po::testing::ObjectiveGradient max_return_gradient(const Data& data)
{
    std::vector<double> gradient = po::testing::excess_mu(data);
    for (double& entry : gradient)
    {
        entry = -entry;
    }
    return po::testing::sum_of({gradient});
}

// Returns are compared relative to the larger of the reference and the largest mean, so that a reference
// near 0 is not compared with its own rounding.
double return_scale(const Data& data, double reference)
{
    return std::max(std::fabs(reference), po::testing::largest_magnitude(data.mu()));
}

void expect_against_reference(const GeneratedProblem& problem, double cap, double tolerance)
{
    double reference = po::testing::max_return_reference(problem.spec, problem.data, cap);
    Result result = po::solve(max_return(problem.spec, cap), problem.data).value();
    if (std::isnan(reference))
    {
        EXPECT_EQ(result.status, SolveStatus::Infeasible);
        return;
    }
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    double achieved = po::testing::portfolio_return(problem.data, result.weights);
    EXPECT_NEAR(achieved, reference, tolerance * return_scale(problem.data, reference));
    po::testing::expect_user_kkt(max_return(problem.spec, cap), problem.data, result,
                                 max_return_gradient(problem.data), tolerance);
    EXPECT_LE(po::testing::portfolio_volatility(problem.data, result.weights), cap * (1.0 + tolerance));
}

// The one known residual of the battery (docs/guide/validation.md): cash, nearly equal means, cond(Sigma) = 1e8
// at the daily scale, where the cap is reached near t = 1e12 and Q(t) is lost in its own rounding. It stops at
// MAX_ITER; no other case may.
struct KnownResidual
{
    std::uint64_t seed;
    double fraction;
};
constexpr KnownResidual kKnownResidual{11394049047870518593ULL, 0.9};

TEST(MaxReturnFuzz, AgainstFrontierReference)
{
    for (const ProblemShape& shape : po::testing::fuzz_shapes(400, 1))
    {
        if (shape.box == ProblemShape::Box::Fixed)
        {
            continue;
        }
        SCOPED_TRACE(describe(shape));
        GeneratedProblem problem = po::testing::generate(shape);
        double minimum = minimum_volatility(problem.spec, problem.data);
        double largest = largest_asset_volatility(problem.data);
        for (double fraction : {0.1, 0.3, 0.5, 0.7, 0.9})
        {
            SCOPED_TRACE(fraction);
            double cap = minimum + fraction * (largest - minimum);
            if (shape.seed == kKnownResidual.seed && fraction == kKnownResidual.fraction)
            {
                EXPECT_EQ(po::solve(max_return(problem.spec, cap), problem.data).value().status,
                          SolveStatus::MaxIterations);
                continue;
            }
            expect_against_reference(problem, cap, condition_tolerance(shape.condition));
        }
    }
}

// Log-uniform in [1e-10, 1e-6].
double riskless_variance(std::uint64_t seed)
{
    return std::pow(10.0, po::testing::random_vector(1, seed, -10.0, -6.0)[0]);
}

// cond(Sigma) of the risky block, whose eigenvalues random_covariance spreads from 0.04 to 0.04 cond (times the
// scale), with the riskless variance appended.
double augmented_condition(const ProblemShape& shape, double variance)
{
    double smallest = 0.04 * shape.scale;
    double largest = smallest * shape.condition;
    return std::max(largest, variance) / std::min(smallest, variance);
}

// A nearly riskless asset of the given variance and a mean of 0.4 times the smallest one, appended after the
// risky assets of problem (before cash), long-only within any box.
GeneratedProblem with_riskless_asset(const GeneratedProblem& problem, double variance)
{
    std::size_t n = problem.data.asset_count();
    std::size_t m = n + 1;
    std::vector<double> covariance(m * m, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        std::copy_n(problem.data.covariance().begin() + static_cast<std::ptrdiff_t>(i * n), n,
                    covariance.begin() + static_cast<std::ptrdiff_t>(i * m));
    }
    covariance[m * m - 1] = variance;
    std::vector<double> mu(problem.data.mu().begin(), problem.data.mu().end());
    mu.push_back(0.4 * *std::ranges::min_element(mu));
    Data data = Data::create(m, mu, covariance).value();
    if (problem.data.risk_free().has_value())
    {
        data.set_risk_free(*problem.data.risk_free());
    }
    Spec spec = problem.spec;
    if (!spec.asset_lower.empty())
    {
        spec.asset_lower.insert(spec.asset_lower.begin() + static_cast<std::ptrdiff_t>(n), 0.0);
        spec.asset_upper.insert(spec.asset_upper.begin() + static_cast<std::ptrdiff_t>(n), 1.0);
    }
    return GeneratedProblem{spec, data};
}

void expect_within_bounds(const Spec& spec, const Result& result, double tolerance)
{
    for (std::size_t k = 0; k < result.weights.size(); ++k)
    {
        SCOPED_TRACE(k);
        EXPECT_GE(result.weights[k], po::testing::lower_of(spec, k) - tolerance);
        EXPECT_LE(result.weights[k], po::testing::upper_of(spec, k) + tolerance);
    }
}

// OPTIMAL only when the reference agrees, inside the bounds and the cap; MAX_ITER otherwise.
void expect_correct_or_stopped(const GeneratedProblem& problem, double cap, double tolerance)
{
    Result result = po::solve(max_return(problem.spec, cap), problem.data).value();
    if (result.status == SolveStatus::MaxIterations)
    {
        return;
    }
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    double reference = po::testing::max_return_reference(problem.spec, problem.data, cap);
    double achieved = po::testing::portfolio_return(problem.data, result.weights);
    EXPECT_NEAR(achieved, reference, tolerance * return_scale(problem.data, reference));
    expect_within_bounds(problem.spec, result, tolerance);
    EXPECT_LE(po::testing::portfolio_volatility(problem.data, result.weights), cap * (1.0 + tolerance));
}

TEST(MaxReturnFuzz, NeverSilentlyWrong)
{
    for (const ProblemShape& shape : po::testing::fuzz_shapes(400, 2))
    {
        if (shape.means != ProblemShape::Means::NearlyEqual)
        {
            continue;
        }
        SCOPED_TRACE(describe(shape));
        double variance = riskless_variance(shape.seed + 7000);
        GeneratedProblem problem = with_riskless_asset(po::testing::generate(shape), variance);
        double tolerance = condition_tolerance(augmented_condition(shape, variance));
        double minimum = minimum_volatility(problem.spec, problem.data);
        double largest = largest_asset_volatility(problem.data);
        std::vector<double> caps;
        for (double fraction : {0.1, 0.5, 0.9})
        {
            caps.push_back(minimum + fraction * (largest - minimum));
        }
        for (double factor : {1.01, 2.0, 10.0})
        {
            caps.push_back(factor * minimum);
        }
        for (double cap : caps)
        {
            SCOPED_TRACE(cap);
            if (cap > 0.0)
            {
                expect_correct_or_stopped(problem, cap, tolerance);
            }
        }
    }
}

} // namespace
