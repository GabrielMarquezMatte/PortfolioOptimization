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

// The generated rows (a - s)'x >= b of the gross exposure family, with their sign vectors s (n entries per row
// each): those of the active set, and while a candidate is added also the candidate and the rows its steps drop.
struct GeneratedRows
{
    std::vector<double> rows;
    std::vector<double> signs;
    std::vector<double> norms;
};

// The explicit rows, equalities first, and after them the generated rows: row index explicit_count + k is the
// generated row k.
struct Constraints
{
    std::size_t n = 0;
    std::size_t equality_count = 0;
    std::size_t explicit_count = 0;
    std::vector<double> rows;
    std::vector<double> rhs;
    std::vector<double> row_norms;
    const GrossExposure* family = nullptr;
    GeneratedRows generated;
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

bool is_generated(const Constraints& constraints, std::size_t index)
{
    return index >= constraints.explicit_count;
}

std::span<const double> row(const Constraints& constraints, std::size_t index)
{
    std::size_t n = constraints.n;
    if (is_generated(constraints, index))
    {
        std::size_t k = index - constraints.explicit_count;
        return std::span<const double>(constraints.generated.rows).subspan(k * n, n);
    }
    return std::span<const double>(constraints.rows).subspan(index * n, n);
}

std::span<const double> signs_of(const Constraints& constraints, std::size_t index)
{
    std::size_t k = index - constraints.explicit_count;
    return std::span<const double>(constraints.generated.signs).subspan(k * constraints.n, constraints.n);
}

double rhs_of(const Constraints& constraints, std::size_t index)
{
    if (is_generated(constraints, index))
    {
        return constraints.family->rhs;
    }
    return constraints.rhs[index];
}

double row_norm_of(const Constraints& constraints, std::size_t index)
{
    if (is_generated(constraints, index))
    {
        return constraints.generated.norms[index - constraints.explicit_count];
    }
    return constraints.row_norms[index];
}

double slack(const Constraints& constraints, const State& state, std::size_t index)
{
    return dot(row(constraints, index), state.x) - rhs_of(constraints, index);
}

double tolerance(const Constraints& constraints, std::size_t index, double x_scale)
{
    return feasibility_tolerance(rhs_of(constraints, index), row_norm_of(constraints, index), x_scale);
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
    constraints.explicit_count = constraints.rhs.size();
    constraints.row_norms.resize(constraints.rhs.size());
    for (std::size_t index = 0; index < constraints.rhs.size(); ++index)
    {
        std::span<const double> normal = row(constraints, index);
        constraints.row_norms[index] = std::sqrt(dot(normal, normal));
    }
    if (problem.gross_exposure.has_value())
    {
        constraints.family = &*problem.gross_exposure;
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
    if (problem.gross_exposure.has_value())
    {
        // At most n generated rows are active at once: the family counts as n rows.
        constraint_count += n;
    }
    state.iteration_limit = static_cast<int>(10 * (n + constraint_count) + 10);
    return state;
}

bool is_active(const State& state, std::size_t index)
{
    return std::ranges::any_of(state.set.active, [index](const ActiveConstraint& constraint)
                               { return constraint.index == index; });
}

struct Violation
{
    std::optional<std::size_t> index;
    double slack = 0.0;
};

Violation most_violated_explicit(const Constraints& constraints, const State& state)
{
    Violation worst;
    double x_scale = max_abs(state.x);
    for (std::size_t index = constraints.equality_count; index < constraints.explicit_count; ++index)
    {
        double value = slack(constraints, state, index);
        bool within_tolerance = value >= -tolerance(constraints, index, x_scale);
        if (within_tolerance || value >= worst.slack || is_active(state, index))
        {
            continue;
        }
        worst = Violation{index, value};
    }
    return worst;
}

// s'x is largest over {-1, 1}^n at s = sign(x), taken as +1 where x_i = 0, so (a - s)'x >= b is the generated row
// most violated at x and its slack a'x - ||x||_1 - b is that of the family.
std::vector<double> sign_vector(std::span<const double> x)
{
    std::vector<double> signs(x.size(), 1.0);
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        if (x[i] < 0.0)
        {
            signs[i] = -1.0;
        }
    }
    return signs;
}

// The generated row (a - s)'x >= b duplicates an active row: a generated row of the same signs, or an explicit row
// with the same normal and right-hand side (an equality that holds a row of the family, as the MAX_RETURN path
// does). Its violation is then the rounding of that row, and the two normals would be dependent.
bool duplicates_active(const Constraints& constraints, const State& state, std::span<const double> signs,
                       std::span<const double> normal)
{
    return std::ranges::any_of(state.set.active,
                               [&](const ActiveConstraint& constraint)
                               {
                                   if (is_generated(constraints, constraint.index))
                                   {
                                       return std::ranges::equal(signs_of(constraints, constraint.index), signs);
                                   }
                                   return rhs_of(constraints, constraint.index) == constraints.family->rhs &&
                                          std::ranges::equal(row(constraints, constraint.index), normal);
                               });
}

std::size_t append_generated(Constraints& constraints, std::span<const double> normal,
                             std::span<const double> signs)
{
    GeneratedRows& generated = constraints.generated;
    generated.rows.insert(generated.rows.end(), normal.begin(), normal.end());
    generated.signs.insert(generated.signs.end(), signs.begin(), signs.end());
    generated.norms.push_back(std::sqrt(dot(normal, normal)));
    return constraints.explicit_count + generated.norms.size() - 1;
}

// The family competes with the explicit rows by the same rule: its slack against the tolerance of the generated
// row. On a tie the explicit row is taken. A generated row that wins is appended after the explicit ones.
std::optional<ActiveConstraint> most_violated_inequality(Constraints& constraints, const State& state)
{
    Violation worst = most_violated_explicit(constraints, state);
    if (constraints.family != nullptr)
    {
        std::vector<double> signs = sign_vector(state.x);
        std::vector<double> normal = constraints.family->offset;
        for (std::size_t i = 0; i < normal.size(); ++i)
        {
            normal[i] -= signs[i];
        }
        double value = dot(normal, state.x) - constraints.family->rhs;
        double norm = std::sqrt(dot(normal, normal));
        double limit = feasibility_tolerance(constraints.family->rhs, norm, max_abs(state.x));
        if (value < -limit && value < worst.slack && !duplicates_active(constraints, state, signs, normal))
        {
            return ActiveConstraint{append_generated(constraints, normal, signs), 1.0, 0.0};
        }
    }
    if (!worst.index.has_value())
    {
        return std::nullopt;
    }
    return ActiveConstraint{*worst.index, 1.0, 0.0};
}

std::optional<ActiveConstraint> next_candidate(Constraints& constraints, State& state)
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

// A generated row that left the active set, or that was redundant and never entered it, is discarded, and the
// family can generate it again. The active generated rows keep their order and are numbered again after the
// explicit rows.
void discard_inactive_generated(Constraints& constraints, State& state)
{
    if (constraints.generated.norms.empty())
    {
        return;
    }
    std::size_t n = constraints.n;
    GeneratedRows previous = std::move(constraints.generated);
    constraints.generated = GeneratedRows{};
    for (ActiveConstraint& constraint : state.set.active)
    {
        if (!is_generated(constraints, constraint.index))
        {
            continue;
        }
        std::size_t k = constraint.index - constraints.explicit_count;
        std::span<const double> normal = std::span<const double>(previous.rows).subspan(k * n, n);
        std::span<const double> signs = std::span<const double>(previous.signs).subspan(k * n, n);
        constraint.index = append_generated(constraints, normal, signs);
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

// The refinement moves x onto the active rows by their residuals, which carry the rounding of every step from
// the unconstrained start: far from it, the move can leave an inactive row violated. Every row, and the family, is
// checked again at the final x, and a violated one is added like any other candidate.
std::optional<ActiveConstraint> violated_after_refinement(Constraints& constraints, State& state)
{
    refine_active(constraints, state);
    return most_violated_inequality(constraints, state);
}

struct GeneratedEntry
{
    std::span<const double> signs;
    double multiplier;
};

// The active generated rows in lexicographic order of their signs, so that equal active sets report equally.
void report_generated(const Constraints& constraints, const State& state, Solution& solution)
{
    std::vector<GeneratedEntry> entries;
    for (const ActiveConstraint& constraint : state.set.active)
    {
        if (is_generated(constraints, constraint.index))
        {
            entries.push_back(GeneratedEntry{signs_of(constraints, constraint.index), constraint.multiplier});
        }
    }
    std::ranges::sort(entries, [](const GeneratedEntry& first, const GeneratedEntry& second)
                      { return std::ranges::lexicographical_compare(first.signs, second.signs); });
    for (const GeneratedEntry& entry : entries)
    {
        solution.gross_exposure_signs.insert(solution.gross_exposure_signs.end(), entry.signs.begin(),
                                             entry.signs.end());
        solution.gross_exposure_multipliers.push_back(entry.multiplier);
    }
}

Solution optimal_solution(const Constraints& constraints, State& state)
{
    Solution solution;
    solution.status = Status::Optimal;
    solution.equality_multipliers.assign(constraints.equality_count, 0.0);
    solution.inequality_multipliers.assign(constraints.explicit_count - constraints.equality_count, 0.0);
    for (const ActiveConstraint& constraint : state.set.active)
    {
        double multiplier = constraint.sign * constraint.multiplier;
        if (constraint.index < constraints.equality_count)
        {
            solution.equality_multipliers[constraint.index] = multiplier;
            continue;
        }
        if (is_generated(constraints, constraint.index))
        {
            continue;
        }
        solution.inequality_multipliers[constraint.index - constraints.equality_count] = multiplier;
        solution.active_inequalities.push_back(constraint.index - constraints.equality_count);
    }
    std::ranges::sort(solution.active_inequalities);
    report_generated(constraints, state, solution);
    solution.x = std::move(state.x);
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
            candidate = violated_after_refinement(constraints, state);
        }
        if (!candidate)
        {
            return optimal_solution(constraints, state);
        }
        progress = add_constraint(constraints, state, *candidate);
        discard_inactive_generated(constraints, state);
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
    if (problem.gross_exposure.has_value() && problem.gross_exposure->offset.size() != n)
    {
        throw std::invalid_argument("qp: the gross exposure offset must have n entries");
    }
}

} // namespace

double feasibility_tolerance(double rhs, double row_norm, double x_scale)
{
    return kRelativeTolerance * std::max({1.0, std::fabs(rhs), row_norm * x_scale});
}

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
