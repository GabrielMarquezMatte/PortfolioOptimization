#include "riskparity/risk_parity.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>

#include "metrics/portfolio.hpp"
#include "problems/constraints.hpp"
#include "problems/mean_variance.hpp"
#include "qp/goldfarb_idnani.hpp"
#include "qp/problem.hpp"
#include "riskparity/augmented_lagrangian.hpp"
#include "riskparity/kkt.hpp"
#include "riskparity/newton_polish.hpp"
#include "riskparity/objective.hpp"
#include "riskparity/polyhedral_descent.hpp"
#include "riskparity/spinu.hpp"

namespace po::riskparity
{
namespace
{

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
// How far the descent (or the augmented Lagrangian) is taken before each polish: the relative KKT residual
// asked of it and the iterations each descent may use.
struct Effort
{
    double tolerance;
    int max_descent_iterations;
};

// First a short descent to a point near enough for Newton; then, when that polish did not end at a KKT
// point, the tightest and longest the descent is asked for.
constexpr std::array<Effort, 2> kStages = {Effort{1e-6, 200}, Effort{1e-11, 2000}};
constexpr int kMaxOuterIterations = 50;

// The descent of one stage read as a candidate KKT point; not solved when a subproblem could not be.
struct Stage
{
    KktPoint point;
    int iterations;
    bool solved;
};

bool inside_box(const Spec& spec, std::span<const double> weights)
{
    for (std::size_t i = 0; i < weights.size(); ++i)
    {
        bool above_lower = weights[i] >= problems::lower_of(spec, i);
        bool below_upper = weights[i] <= problems::upper_of(spec, i);
        if (!above_lower || !below_upper)
        {
            return false;
        }
    }
    return true;
}

Result no_point(std::size_t n, SolveStatus status, int iterations)
{
    Result result;
    result.status = status;
    result.weights.assign(n, kNaN);
    result.bound_multipliers.assign(n, kNaN);
    result.budget_multiplier = kNaN;
    result.residual = kNaN;
    result.iterations = iterations;
    return result;
}

Result at_weights(const Data& data, std::vector<double> weights, std::span<const double> budget,
                  int iterations)
{
    Result result;
    result.status = SolveStatus::Optimal;
    result.residual = metrics::risk_budget_residual(data, weights, budget);
    result.iterations = iterations;
    result.weights = std::move(weights);
    return result;
}

Result spinu_point(const Data& data, std::vector<double> weights, std::span<const double> budget,
                   int iterations)
{
    Result result = at_weights(data, std::move(weights), budget, iterations);
    result.budget_multiplier = 0.0;
    result.bound_multipliers.assign(result.weights.size(), 0.0);
    return result;
}

// At s = sigma_min the minimum-variance portfolio is the only feasible point and the gradient of the
// volatility there is a combination of those of the other active rows: in general no multipliers satisfy the
// stationarity, and the optimal f falls with infinite slope as s grows. The same convention as MAX_RETURN on
// a vertical frontier: a volatility multiplier of -infinity and no value for the others.
Result minimum_volatility_point(const Data& data, std::vector<double> weights, std::span<const double> budget,
                                int iterations)
{
    Result result = at_weights(data, std::move(weights), budget, iterations);
    result.bound_multipliers.assign(result.weights.size(), kNaN);
    result.volatility_multiplier = -std::numeric_limits<double>::infinity();
    return result;
}

Result kkt_point(const KktProblem& problem, KktPoint point, int iterations)
{
    Result result = at_weights(problem.data, std::move(point.weights), problem.budget, iterations);
    result.budget_multiplier = point.equality_multipliers.front();
    if (problem.spec.target_return.has_value())
    {
        result.target_return_multiplier = point.equality_multipliers[1];
    }
    if (problem.spec.target_volatility.has_value())
    {
        result.volatility_multiplier = point.equality_multipliers.back();
    }
    result.bound_multipliers = std::move(point.bound_multipliers);
    return result;
}

// One asset without a target volatility: the minimum-variance portfolio w = 1 is the only feasible point, where
// c = 1 = b and grad f = 0, so every multiplier 0 is a valid set (the budget and a target return are then the
// same row). The descent cannot run there: J = 0 leaves it no curvature.
Result single_asset(const KktProblem& problem, std::vector<double> weights, int iterations)
{
    Result result = at_weights(problem.data, std::move(weights), problem.budget, iterations);
    result.budget_multiplier = 0.0;
    if (problem.spec.target_return.has_value())
    {
        result.target_return_multiplier = 0.0;
    }
    result.bound_multipliers.assign(result.weights.size(), 0.0);
    return result;
}

// The minimum-variance portfolio over the budget, the target return and the bounds: it decides whether the
// linear rows are feasible and gives the smallest volatility they allow.
Result minimum_variance(const Spec& spec, const Data& data)
{
    Spec linear = spec;
    linear.objective = Objective::MinVariance;
    return problems::solve_mean_variance(linear, data);
}

double volatility_of(const Data& data, std::span<const double> weights)
{
    return std::sqrt(compensated_dot(weights, covariance_times(data, weights)));
}

// The point of the polyhedron nearest to target: min 1/2 |x - target|^2 over the budget, the target return
// and the bounds.
qp::Solution nearest_in_polyhedron(const Spec& spec, const Data& data, std::span<const double> target)
{
    std::size_t n = target.size();
    qp::Problem problem;
    problem.n = n;
    problem.hessian.assign(n * n, 0.0);
    problem.linear.resize(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        problem.hessian[i * n + i] = 1.0;
        problem.linear[i] = -target[i];
    }
    problem.hessian_factor = problem.hessian;
    problems::add_constraint_rows(spec, data, problem);
    return qp::solve(problem);
}

// A descent whose last subproblem could not be solved leaves NaN multipliers.
bool has_multipliers(const KktPoint& point)
{
    return std::ranges::none_of(point.equality_multipliers, [](double value) { return std::isnan(value); });
}

Stage descent_stage(const KktProblem& problem, std::span<const double> start, const Effort& effort)
{
    DescentProblem descent{problem.spec, problem.data, problem.budget, std::nullopt};
    DescentOptions options{kDescentProximal, effort.tolerance, effort.max_descent_iterations};
    DescentResult result = descend(descent, start, options);
    KktPoint point{std::move(result.weights), std::move(result.equality_multipliers),
                   std::move(result.bound_multipliers)};
    bool solved = has_multipliers(point);
    return Stage{std::move(point), result.iterations, solved};
}

Stage augmented_stage(const KktProblem& problem, std::span<const double> start, const Effort& effort)
{
    AugmentedOptions options{effort.tolerance, kMaxOuterIterations, effort.max_descent_iterations};
    AugmentedResult result = augmented_lagrangian(problem, start, options);
    bool solved = has_multipliers(result.point);
    return Stage{std::move(result.point), result.inner_iterations, solved};
}

Stage stage_of(const KktProblem& problem, std::span<const double> start, const Effort& effort)
{
    if (problem.spec.target_volatility.has_value())
    {
        return augmented_stage(problem, start, effort);
    }
    return descent_stage(problem, start, effort);
}

// Each stage takes the descent to its tolerance, from where the last one stopped, and polishes the point by
// Newton on the KKT system. Only a KKT point to the noise of the problem is an answer, judged with the
// least-squares multipliers of the polish; the answer carries its corrected multipliers. When no stage ends
// at one, the answer is MaxIterations, without a point.
Result stationary_point(const KktProblem& problem, std::vector<double> start, int iterations)
{
    for (const Effort& effort : kStages)
    {
        Stage stage = stage_of(problem, start, effort);
        iterations += stage.iterations;
        if (!stage.solved)
        {
            break;
        }
        PolishResult polished = newton_polish(problem, stage.point);
        iterations += polished.steps;
        if (is_kkt_point(problem, polished.least_squares_point))
        {
            MultiplierCorrection correction = corrected_multipliers(problem, polished.least_squares_point);
            return kkt_point(problem, std::move(correction.reported_point), iterations);
        }
        start = std::move(stage.point.weights);
    }
    return no_point(start.size(), SolveStatus::MaxIterations, iterations);
}

// The general path, from a minimum-variance portfolio that exists: the certificates of the target
// volatility, then the descent from the point of the polyhedron nearest to the Spinu solution.
Result constrained(const KktProblem& problem, const Result& minimum, const SpinuResult& spinu_result)
{
    std::size_t n = problem.data.asset_count();
    int iterations = spinu_result.iterations;
    if (problem.spec.target_volatility.has_value())
    {
        double target = *problem.spec.target_volatility;
        double smallest = volatility_of(problem.data, minimum.weights);
        if (target < smallest * (1.0 - problems::kMinimumVolatilityTolerance))
        {
            return no_point(n, SolveStatus::Infeasible, iterations);
        }
        if (std::fabs(target - smallest) <= problems::kMinimumVolatilityTolerance * target)
        {
            return minimum_volatility_point(problem.data, minimum.weights, problem.budget, iterations);
        }
    }
    if (n == 1 && !problem.spec.target_volatility.has_value())
    {
        return single_asset(problem, minimum.weights, iterations);
    }
    qp::Solution start = nearest_in_polyhedron(problem.spec, problem.data, spinu_result.weights);
    if (start.status != qp::Status::Optimal)
    {
        return no_point(n, SolveStatus::MaxIterations, iterations);
    }
    return stationary_point(problem, std::move(start.x), iterations);
}

} // namespace

std::vector<double> risk_budget_of(const Spec& spec, const Data& data)
{
    if (spec.risk_budget.empty())
    {
        return std::vector<double>(data.asset_count(), 1.0 / static_cast<double>(data.asset_count()));
    }
    return spec.risk_budget;
}

std::expected<Result, Failure> solve_risk_parity(const Spec& spec, const Data& data)
{
    std::vector<double> budget = risk_budget_of(spec, data);
    SpinuResult spinu_result = spinu(data, budget);
    if (!spinu_result.converged)
    {
        return no_point(data.asset_count(), SolveStatus::MaxIterations, spinu_result.iterations);
    }
    bool has_target = spec.target_return.has_value() || spec.target_volatility.has_value();
    if (!has_target && inside_box(spec, spinu_result.weights))
    {
        return spinu_point(data, std::move(spinu_result.weights), budget, spinu_result.iterations);
    }
    Result minimum = minimum_variance(spec, data);
    if (minimum.status != SolveStatus::Optimal)
    {
        return no_point(data.asset_count(), minimum.status, spinu_result.iterations);
    }
    return constrained(KktProblem{spec, data, budget}, minimum, spinu_result);
}

} // namespace po::riskparity
