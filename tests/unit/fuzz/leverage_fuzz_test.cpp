#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "core/result.hpp"
#include "fuzz/random_problems.hpp"
#include "fuzz/sign_pattern_reference.hpp"
#include "problems/mean_variance.hpp"
#include "solve/solve.hpp"
#include "problems/tracking_error.hpp"
#include "qp/problem.hpp"
#include "random_covariance.hpp"
#include "user_kkt.hpp"

namespace
{

using po::Data;
using po::Objective;
using po::Result;
using po::SolveStatus;
using po::Spec;
using po::testing::GeneratedProblem;
using po::testing::ProblemShape;
using po::testing::condition_tolerance;
using po::testing::describe;

// The QP objectives checked here against the sign patterns; MAX_RETURN is checked by MaxReturnFuzz.LeverageFuzz.
constexpr std::array kObjectives = {Objective::MinVariance, Objective::Utility, Objective::Importance,
                                    Objective::Kelly, Objective::TrackingError};

// Returns R = L' + rf over T = n periods, with L the factor of Sigma, so that R~ = L' and R~'R~ = Sigma: the
// TRACKING_ERROR QP has the condition number of the shape. The benchmark is R~ w_b + rf plus noise, with w_b
// uniform in [-1, 1] shifted to 1'w_b = 1, whose gross exposure is often above the leverage.
Data with_returns(Data data, double scale, std::uint64_t seed)
{
    std::size_t n = data.asset_count();
    double rate = data.risk_free().value_or(0.0);
    std::span<const double> factor = data.covariance_factor();
    std::vector<double> target = po::testing::random_vector(n, seed, -1.0, 1.0);
    double shift = 1.0;
    for (double weight : target)
    {
        shift -= weight;
    }
    for (double& weight : target)
    {
        weight += shift / static_cast<double>(n);
    }
    std::vector<double> noise = po::testing::random_vector(n, seed + 1, -0.01, 0.01);
    std::vector<double> returns(n * n);
    std::vector<double> benchmark(n, rate);
    for (std::size_t t = 0; t < n; ++t)
    {
        benchmark[t] += noise[t] * std::sqrt(scale);
        for (std::size_t i = 0; i < n; ++i)
        {
            returns[t * n + i] = factor[i * n + t] + rate;
            benchmark[t] += factor[i * n + t] * target[i];
        }
    }
    EXPECT_TRUE(data.set_returns(returns, benchmark, n).has_value());
    return data;
}

// grad f = quadratic Sigma w - linear (mu - rf 1) for the four mean-variance objectives at their default
// parameters, written here independently of the builder.
po::testing::ObjectiveGradient mean_variance_gradient(const Data& data, const Result& result,
                                                      Objective objective)
{
    double rate = data.risk_free().value_or(0.0);
    double quadratic = 2.0;
    double linear = 0.0;
    if (objective == Objective::Utility)
    {
        quadratic = 5.0;
        linear = 1.0;
    }
    if (objective == Objective::Importance)
    {
        quadratic = 1.0;
        linear = 0.5;
    }
    if (objective == Objective::Kelly)
    {
        quadratic = 1.0 / (1.0 + rate);
        linear = 1.0;
    }
    std::vector<double> curvature = po::testing::covariance_times(data, result.weights);
    std::vector<double> excess = po::testing::excess_mu(data);
    for (std::size_t i = 0; i < curvature.size(); ++i)
    {
        curvature[i] *= quadratic;
        excess[i] *= -linear;
    }
    return po::testing::sum_of({curvature, excess});
}

// grad f = 2 R~'(R~ w + rf 1 - b) as its two terms 2 R~'R~ w and -2 R~'(b - rf 1).
po::testing::ObjectiveGradient tracking_gradient(const Data& data, const Result& result)
{
    std::size_t n = data.asset_count();
    double rate = data.risk_free().value_or(0.0);
    std::vector<double> curvature(n, 0.0);
    std::vector<double> linear(n, 0.0);
    for (std::size_t t = 0; t < data.observation_count(); ++t)
    {
        std::span<const double> row = data.returns().subspan(t * n, n);
        double fitted = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            fitted += (row[i] - rate) * result.weights[i];
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            curvature[i] += 2.0 * (row[i] - rate) * fitted;
            linear[i] -= 2.0 * (row[i] - rate) * (data.benchmark()[t] - rate);
        }
    }
    return po::testing::sum_of({curvature, linear});
}

po::testing::ObjectiveGradient objective_gradient(const Data& data, const Result& result, Objective objective)
{
    if (objective == Objective::TrackingError)
    {
        return tracking_gradient(data, result);
    }
    return mean_variance_gradient(data, result, objective);
}

// The builder's QP of spec without the leverage, whose objective the sign patterns minimize.
po::qp::Problem problem_without_family(Spec spec, const Data& data)
{
    spec.leverage.reset();
    if (spec.objective == Objective::TrackingError)
    {
        return po::problems::tracking_error_problem(spec, data).value().problem;
    }
    po::problems::MeanVariance scales = po::problems::mean_variance_of(spec, data);
    return po::problems::mean_variance_problem(spec, data, scales).problem;
}

// The value is compared relative to the larger of the magnitude of its terms and |g|max |w|max over every user
// weight, cash included, so that an optimum at x = 0 (all cash) is not compared with its own rounding.
double value_scale(const po::qp::Problem& problem, const Result& result)
{
    std::span<const double> x = std::span<const double>(result.weights).first(problem.n);
    double linear =
        po::testing::largest_magnitude(problem.linear) * po::testing::largest_magnitude(result.weights);
    return std::fmax(po::testing::qp_objective_scale(problem, x), linear);
}

// INFEASIBLE exactly when no sign pattern is Optimal; otherwise OPTIMAL with the value of the reference and KKT
// in the user's variables. Returns whether the leverage binds (lambda_lev < 0).
bool binds_after_checks(const Spec& spec, const Data& data, double tolerance)
{
    po::qp::Problem problem = problem_without_family(spec, data);
    std::optional<double> reference = po::testing::sign_pattern_optimum(problem, *spec.leverage);
    Result result = po::solve(spec, data).value();
    if (!reference.has_value())
    {
        EXPECT_EQ(result.status, SolveStatus::Infeasible);
        return false;
    }
    EXPECT_EQ(result.status, SolveStatus::Optimal);
    if (result.status != SolveStatus::Optimal)
    {
        return false;
    }
    std::span<const double> x = std::span<const double>(result.weights).first(data.asset_count());
    EXPECT_NEAR(po::testing::qp_objective(problem, x), *reference, tolerance * value_scale(problem, result));
    po::testing::ObjectiveGradient gradient = objective_gradient(data, result, spec.objective);
    po::testing::expect_user_kkt(spec, data, result, gradient, tolerance);
    return result.leverage_multiplier < 0.0;
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

// A leverage below the gross exposure of the optimum without leverage: a fraction in [0.2, 0.9], drawn from seed,
// of the way to it from the least the budget allows (1 without cash, 0 with cash); none when that optimum is within
// 1e-6 of the least (a long-only box without cash, all cash). It binds, or it is infeasible where the bounds force
// more (the cash cap of the tight box, the fixed weights).
std::optional<double> binding_leverage(Spec spec, const Data& data, std::uint64_t seed)
{
    spec.leverage.reset();
    Result free = po::solve(spec, data).value();
    double least = 1.0;
    if (data.risk_free().has_value())
    {
        least = 0.0;
    }
    double gross = risky_gross_exposure(data, free);
    if (free.status != SolveStatus::Optimal || !(gross - least > 1e-6))
    {
        return std::nullopt;
    }
    return least + po::testing::random_vector(1, seed, 0.2, 0.9)[0] * (gross - least);
}

// Every case at the leverage of its shape (uniform in [0.8, 3]: infeasible, inactive or binding) and, where one
// exists, at a leverage below the optimum's gross exposure: most shapes have a box under which the leverage of the
// shape never binds. Of the 1979 cases, 699 bind (docs/guide/validation.md).
TEST(LeverageFuzz, AgainstSignPatterns)
{
    std::size_t binding = 0;
    for (const ProblemShape& shape : po::testing::with_leverage(po::testing::fuzz_shapes(256, 4, 8)))
    {
        SCOPED_TRACE(describe(shape));
        GeneratedProblem problem = po::testing::generate(shape);
        Data data = with_returns(problem.data, shape.scale, shape.seed + 900);
        double tolerance = condition_tolerance(shape.condition);
        for (Objective objective : kObjectives)
        {
            SCOPED_TRACE(static_cast<int>(objective));
            Spec spec = problem.spec;
            spec.objective = objective;
            binding += binds_after_checks(spec, data, tolerance) ? 1 : 0;
            spec.leverage = binding_leverage(spec, data, shape.seed + 1000 + static_cast<std::uint64_t>(objective));
            if (spec.leverage.has_value())
            {
                SCOPED_TRACE(testing::Message() << "binding leverage " << *spec.leverage);
                binding += binds_after_checks(spec, data, tolerance) ? 1 : 0;
            }
        }
    }
    EXPECT_GE(binding, 600u);
}

} // namespace
