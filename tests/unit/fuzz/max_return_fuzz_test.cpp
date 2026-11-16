#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
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
using po::testing::condition_tolerance;
using po::testing::describe;

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

// INFEASIBLE exactly when the reference is NaN; otherwise OPTIMAL with the return of the reference, KKT and
// vol <= s. Returns the result when it is OPTIMAL.
std::optional<Result> checked_against_reference(const GeneratedProblem& problem, double cap, double tolerance)
{
    double reference = po::testing::max_return_reference(problem.spec, problem.data, cap);
    Result result = po::solve(max_return(problem.spec, cap), problem.data).value();
    if (std::isnan(reference))
    {
        EXPECT_EQ(result.status, SolveStatus::Infeasible);
        return std::nullopt;
    }
    EXPECT_EQ(result.status, SolveStatus::Optimal);
    if (result.status != SolveStatus::Optimal)
    {
        return std::nullopt;
    }
    double achieved = po::testing::portfolio_return(problem.data, result.weights);
    EXPECT_NEAR(achieved, reference, tolerance * return_scale(problem.data, reference));
    po::testing::expect_user_kkt(max_return(problem.spec, cap), problem.data, result,
                                 max_return_gradient(problem.data), tolerance);
    EXPECT_LE(po::testing::portfolio_volatility(problem.data, result.weights), cap * (1.0 + tolerance));
    return result;
}

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
            checked_against_reference(problem, cap, condition_tolerance(shape.condition));
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

double risky_gross_exposure(const Data& data, const Result& result)
{
    double total = 0.0;
    for (double weight : std::span<const double>(result.weights).first(data.asset_count()))
    {
        total += std::fabs(weight);
    }
    return total;
}

// A leverage below the gross exposure of the MAX_RETURN optimum at cap without leverage: a fraction in [0.2, 0.9],
// drawn from seed, of the way to it from the least the budget allows (1 without cash, 0 with cash); none when that
// optimum is not OPTIMAL or within 1e-6 of the least. It binds, or the cap is below the minimum volatility under it.
std::optional<double> binding_leverage(const GeneratedProblem& problem, double cap, std::uint64_t seed)
{
    Spec spec = problem.spec;
    spec.leverage.reset();
    Result free = po::solve(max_return(spec, cap), problem.data).value();
    if (free.status != SolveStatus::Optimal)
    {
        return std::nullopt;
    }
    double least = 1.0;
    if (problem.data.risk_free().has_value())
    {
        least = 0.0;
    }
    double gross = risky_gross_exposure(problem.data, free);
    if (!(gross - least > 1e-6))
    {
        return std::nullopt;
    }
    return least + po::testing::random_vector(1, seed, 0.2, 0.9)[0] * (gross - least);
}

struct LeverageCount
{
    std::size_t leverage = 0;
    std::size_t both = 0;
};

void count_binding(const std::optional<Result>& result, LeverageCount& count)
{
    if (!result.has_value() || !(result->leverage_multiplier < 0.0))
    {
        return;
    }
    ++count.leverage;
    if (result->volatility_multiplier < 0.0)
    {
        ++count.both;
    }
}

// gcc-release runs LeverageFuzz in full. gcc-debug, where it took 117 s, runs the first 256 shapes (every combination
// once) at three of the five caps, 38 s, where 407 cases bind the leverage, 278 with the cap
// (docs/guide/validation.md).
#ifdef NDEBUG
constexpr std::size_t kLeverageShapes = 400;
constexpr std::array kCapFractions = {0.1, 0.3, 0.5, 0.7, 0.9};
constexpr std::size_t kLeverageBinds = 1000;
constexpr std::size_t kBothBind = 700;
#endif
#ifndef NDEBUG
constexpr std::size_t kLeverageShapes = 256;
constexpr std::array kCapFractions = {0.1, 0.5, 0.9};
constexpr std::size_t kLeverageBinds = 400;
constexpr std::size_t kBothBind = 270;
#endif

// The caps of AgainstFrontierReference (kCapFractions of them), from the minimum volatility under the leverage of
// spec; the largest asset volatility alone when no portfolio is feasible under it.
std::vector<double> caps_of(const GeneratedProblem& problem)
{
    double largest = largest_asset_volatility(problem.data);
    std::vector<double> caps;
    Result minimum = po::testing::frontier_point(problem.spec, problem.data, std::nullopt);
    if (minimum.status != SolveStatus::Optimal)
    {
        caps.push_back(largest);
        return caps;
    }
    double lowest = po::testing::portfolio_volatility(problem.data, minimum.weights);
    for (double fraction : kCapFractions)
    {
        caps.push_back(lowest + fraction * (largest - lowest));
    }
    return caps;
}

// Every shape of fuzz_shapes(kLeverageShapes, 3) without the fixed box, at the leverage of with_leverage (uniform in [0.8, 3]:
// infeasible, inactive or binding) and at a leverage below the gross exposure of the optimum without it, against
// the frontier reference whose MIN_VARIANCE has the same leverage. In gcc-release, of the 2314 cases, 1057 bind the
// leverage, 763 of them with the cap (docs/guide/validation.md).
TEST(MaxReturnFuzz, LeverageFuzz)
{
    LeverageCount count;
    for (const ProblemShape& shape : po::testing::with_leverage(po::testing::fuzz_shapes(kLeverageShapes, 3)))
    {
        if (shape.box == ProblemShape::Box::Fixed)
        {
            continue;
        }
        SCOPED_TRACE(describe(shape));
        GeneratedProblem problem = po::testing::generate(shape);
        double tolerance = condition_tolerance(shape.condition);
        for (double cap : caps_of(problem))
        {
            SCOPED_TRACE(cap);
            count_binding(checked_against_reference(problem, cap, tolerance), count);
            GeneratedProblem binding = problem;
            binding.spec.leverage = binding_leverage(problem, cap, shape.seed + 1000);
            if (binding.spec.leverage.has_value())
            {
                SCOPED_TRACE(testing::Message() << "binding leverage " << *binding.spec.leverage);
                count_binding(checked_against_reference(binding, cap, tolerance), count);
            }
        }
    }
    EXPECT_GE(count.leverage, kLeverageBinds);
    EXPECT_GE(count.both, kBothBind);
}

} // namespace
