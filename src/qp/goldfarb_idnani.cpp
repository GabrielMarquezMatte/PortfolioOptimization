#include "qp/goldfarb_idnani.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "linalg/cholesky.hpp"
#include "qp/active_set.hpp"

namespace po::qp
{

namespace
{

constexpr double kRelativeTolerance = 1e-12;
constexpr double kInfinity = std::numeric_limits<double>::infinity();
constexpr double kNotANumber = std::numeric_limits<double>::quiet_NaN();

struct Constraints
{
    std::size_t n = 0;
    std::size_t equality_count = 0;
    std::vector<double> rows;
    std::vector<double> rhs;
    std::vector<double> row_norms;
};

// G x + g = N u holds for the active normals N and multipliers u of set.active.
struct State
{
    std::vector<double> x;
    ActiveSet set;
    std::size_t next_equality = 0;
    int iterations = 0;
    int iteration_limit = 0;
};

struct Step
{
    double length = kInfinity;
    bool adds = true;
    std::size_t drop_position = 0;
};

enum class Progress
{
    Added,
    Infeasible,
    MaxIterations
};

double max_abs(std::span<const double> values)
{
    double largest = 0.0;
    for (double value : values)
    {
        largest = std::max(largest, std::fabs(value));
    }
    return largest;
}

std::span<const double> row(const Constraints& constraints, std::size_t index)
{
    return std::span<const double>(constraints.rows).subspan(index * constraints.n, constraints.n);
}

double slack(const Constraints& constraints, const State& state, std::size_t index)
{
    return dot(row(constraints, index), state.x) - constraints.rhs[index];
}

double tolerance(const Constraints& constraints, std::size_t index, double x_scale)
{
    double row_scale = constraints.row_norms[index] * x_scale;
    return kRelativeTolerance * std::max({1.0, std::fabs(constraints.rhs[index]), row_scale});
}

Constraints gather_constraints(const Problem& problem)
{
    Constraints constraints;
    constraints.n = problem.n;
    constraints.equality_count = problem.equality_rhs.size();
    constraints.rows = problem.equality_rows;
    constraints.rows.insert(constraints.rows.end(), problem.inequality_rows.begin(),
                            problem.inequality_rows.end());
    constraints.rhs = problem.equality_rhs;
    constraints.rhs.insert(constraints.rhs.end(), problem.inequality_rhs.begin(),
                           problem.inequality_rhs.end());
    constraints.row_norms.resize(constraints.rhs.size());
    for (std::size_t index = 0; index < constraints.rhs.size(); ++index)
    {
        std::span<const double> normal = row(constraints, index);
        constraints.row_norms[index] = std::sqrt(dot(normal, normal));
    }
    return constraints;
}

State initial_state(const Problem& problem, std::span<const double> factor)
{
    std::size_t n = problem.n;
    State state;
    state.x.resize(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        state.x[i] = -problem.linear[i];
    }
    linalg::solve_lower(factor, n, state.x);
    linalg::solve_lower_transposed(factor, n, state.x);
    state.set = empty_active_set(factor, n);
    std::size_t constraint_count = problem.equality_rhs.size() + problem.inequality_rhs.size();
    state.iteration_limit = static_cast<int>(10 * (n + constraint_count) + 10);
    return state;
}

bool is_active(const State& state, std::size_t index)
{
    return std::ranges::any_of(state.set.active, [index](const ActiveConstraint& constraint)
                               { return constraint.index == index; });
}

std::optional<ActiveConstraint> most_violated_inequality(const Constraints& constraints, const State& state)
{
    std::optional<ActiveConstraint> candidate;
    double x_scale = max_abs(state.x);
    double worst = 0.0;
    for (std::size_t index = constraints.equality_count; index < constraints.rhs.size(); ++index)
    {
        double value = slack(constraints, state, index);
        if (value >= -tolerance(constraints, index, x_scale) || value >= worst || is_active(state, index))
        {
            continue;
        }
        worst = value;
        candidate = ActiveConstraint{index, 1.0, 0.0};
    }
    return candidate;
}

std::optional<ActiveConstraint> next_candidate(const Constraints& constraints, State& state)
{
    if (state.next_equality == constraints.equality_count)
    {
        return most_violated_inequality(constraints, state);
    }
    std::size_t index = state.next_equality;
    ++state.next_equality;
    double sign = 1.0;
    if (slack(constraints, state, index) >= 0.0)
    {
        sign = -1.0;
    }
    return ActiveConstraint{index, sign, 0.0};
}

double primal_step_length(const Constraints& constraints, const State& state,
                          const ActiveConstraint& candidate)
{
    double violation = candidate.sign * slack(constraints, state, candidate.index);
    double curvature = candidate.sign * dot(state.set.primal_direction, row(constraints, candidate.index));
    return -violation / curvature;
}

Step choose_step(const Constraints& constraints, const State& state, double primal_length)
{
    Step step;
    step.length = primal_length;
    for (std::size_t position = 0; position < state.set.active.size(); ++position)
    {
        const ActiveConstraint& constraint = state.set.active[position];
        double rate = state.set.dual_direction[position];
        if (constraint.index < constraints.equality_count || !(rate > 0.0))
        {
            continue;
        }
        double ratio = constraint.multiplier / rate;
        if (ratio < step.length)
        {
            step = Step{ratio, false, position};
        }
    }
    return step;
}

void take_step(State& state, ActiveConstraint& candidate, double length, bool moves_primal)
{
    for (std::size_t position = 0; position < state.set.active.size(); ++position)
    {
        state.set.active[position].multiplier -= length * state.set.dual_direction[position];
    }
    candidate.multiplier += length;
    if (!moves_primal)
    {
        return;
    }
    for (std::size_t i = 0; i < state.x.size(); ++i)
    {
        state.x[i] += length * state.set.primal_direction[i];
    }
}

bool is_independent_of_active(const ActiveSet& set)
{
    std::span<const double> projection = set.projection;
    std::span<const double> outside_active = projection.subspan(set.active.size());
    double outside_norm = std::sqrt(dot(outside_active, outside_active));
    return outside_norm > kRelativeTolerance * std::sqrt(dot(projection, projection));
}

Progress add_constraint(const Constraints& constraints, State& state, ActiveConstraint candidate)
{
    while (true)
    {
        compute_directions(state.set, row(constraints, candidate.index), candidate.sign);
        double violation = slack(constraints, state, candidate.index);
        double limit = tolerance(constraints, candidate.index, max_abs(state.x));
        bool moves_primal = is_independent_of_active(state.set);
        if (!moves_primal && std::fabs(violation) <= limit)
        {
            return Progress::Added;
        }
        double primal_length = kInfinity;
        if (moves_primal)
        {
            primal_length = primal_step_length(constraints, state, candidate);
        }
        Step step = choose_step(constraints, state, primal_length);
        if (step.length == kInfinity)
        {
            return Progress::Infeasible;
        }
        if (state.iterations == state.iteration_limit)
        {
            return Progress::MaxIterations;
        }
        ++state.iterations;
        take_step(state, candidate, step.length, moves_primal);
        if (step.adds)
        {
            append_active(state.set, candidate);
            return Progress::Added;
        }
        drop_active(state.set, step.drop_position);
    }
}

void refine_active(const Constraints& constraints, State& state)
{
    std::vector<double> residual(state.set.active.size());
    for (std::size_t k = 0; k < residual.size(); ++k)
    {
        const ActiveConstraint& constraint = state.set.active[k];
        residual[k] = -constraint.sign * slack(constraints, state, constraint.index);
    }
    refine_toward_active(state.set, residual, state.x);
}

Solution optimal_solution(const Constraints& constraints, State& state)
{
    refine_active(constraints, state);
    Solution solution;
    solution.status = Status::Optimal;
    solution.x = std::move(state.x);
    solution.equality_multipliers.assign(constraints.equality_count, 0.0);
    solution.inequality_multipliers.assign(constraints.rhs.size() - constraints.equality_count, 0.0);
    for (const ActiveConstraint& constraint : state.set.active)
    {
        double multiplier = constraint.sign * constraint.multiplier;
        if (constraint.index < constraints.equality_count)
        {
            solution.equality_multipliers[constraint.index] = multiplier;
            continue;
        }
        solution.inequality_multipliers[constraint.index - constraints.equality_count] = multiplier;
    }
    solution.iterations = state.iterations;
    return solution;
}

Solution failed_solution(const Problem& problem, Progress progress, int iterations)
{
    Solution solution;
    solution.status = Status::Infeasible;
    if (progress == Progress::MaxIterations)
    {
        solution.status = Status::MaxIterations;
    }
    solution.x.assign(problem.n, kNotANumber);
    solution.equality_multipliers.assign(problem.equality_rhs.size(), kNotANumber);
    solution.inequality_multipliers.assign(problem.inequality_rhs.size(), kNotANumber);
    solution.iterations = iterations;
    return solution;
}

Solution solve_factored(const Problem& problem, std::span<const double> factor)
{
    Constraints constraints = gather_constraints(problem);
    State state = initial_state(problem, factor);
    Progress progress = Progress::Added;
    while (progress == Progress::Added)
    {
        std::optional<ActiveConstraint> candidate = next_candidate(constraints, state);
        if (!candidate)
        {
            return optimal_solution(constraints, state);
        }
        progress = add_constraint(constraints, state, *candidate);
    }
    return failed_solution(problem, progress, state.iterations);
}

void validate_sizes(const Problem& problem)
{
    std::size_t n = problem.n;
    if (problem.hessian.size() != n * n || problem.linear.size() != n)
    {
        throw std::invalid_argument("qp: the hessian must be n x n and the linear term must have n entries");
    }
    if (!problem.hessian_factor.empty() && problem.hessian_factor.size() != n * n)
    {
        throw std::invalid_argument("qp: the hessian factor must be empty or n x n");
    }
    if (problem.equality_rows.size() != problem.equality_rhs.size() * n)
    {
        throw std::invalid_argument("qp: the equality rows must have n entries per right-hand side");
    }
    if (problem.inequality_rows.size() != problem.inequality_rhs.size() * n)
    {
        throw std::invalid_argument("qp: the inequality rows must have n entries per right-hand side");
    }
}

} // namespace

Solution solve(const Problem& problem)
{
    validate_sizes(problem);
    if (!problem.hessian_factor.empty())
    {
        return solve_factored(problem, problem.hessian_factor);
    }
    std::optional<std::vector<double>> factor = linalg::cholesky(problem.hessian, problem.n);
    if (!factor)
    {
        throw std::invalid_argument("qp: the hessian is not positive definite");
    }
    return solve_factored(problem, *factor);
}

} // namespace po::qp
