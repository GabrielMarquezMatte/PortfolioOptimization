#include "qp/goldfarb_idnani.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
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

constexpr double kInfinity = std::numeric_limits<double>::infinity();
constexpr double kNotANumber = std::numeric_limits<double>::quiet_NaN();
constexpr double kFaceRoundingFactor = 4.0;

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
    // The column of the only nonzero entry of each explicit row (a bound), or n for a row with more than one.
    std::vector<std::size_t> unit_columns;
    const GrossExposure* family = nullptr;
    GeneratedRows generated;
    // Buffers kept between candidates so that the family does not allocate in the loop: the sign vector and the
    // normal of the most violated row of the family, and the rows discard_inactive_generated renumbers into.
    std::vector<double> family_signs;
    std::vector<double> family_normal;
    GeneratedRows renumbered;
};

// How a violated candidate dependent on the active rows is taken: always by a dual step, or passed over when it
// holds on the face of the active rows (holds_on_active_face).
enum class DependentCandidates
{
    DualStep,
    PassOverHeldOnFace
};

// G x + g = N u holds for the active normals N and multipliers u of set.active. passed_over holds the candidates
// passed over since the last step (pass_over), which the scan does not pick again until a step changes the face.
struct State
{
    std::vector<double> x;
    ActiveSet set;
    std::size_t next_equality = 0;
    int iterations = 0;
    int iteration_limit = 0;
    DependentCandidates dependent = DependentCandidates::DualStep;
    std::vector<std::size_t> passed_over;
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

bool is_unit(const Constraints& constraints, std::size_t index)
{
    return !is_generated(constraints, index) && constraints.unit_columns[index] < constraints.n;
}

double unit_entry(const Constraints& constraints, std::size_t index)
{
    return constraints.rows[index * constraints.n + constraints.unit_columns[index]];
}

// The slack of an explicit row, read at the one nonzero entry of a unit row when x is finite (unit_row_dot).
double scanned_slack(const Constraints& constraints, const State& state, std::size_t index, bool finite_point)
{
    if (!finite_point || !is_unit(constraints, index))
    {
        return slack(constraints, state, index);
    }
    std::size_t column = constraints.unit_columns[index];
    return unit_row_dot(state.x, column, unit_entry(constraints, index)) - constraints.rhs[index];
}

// The column of the only nonzero entry of normal, 0 for a row of zeros, and n for a row with more than one.
std::size_t unit_column_of(std::span<const double> normal)
{
    std::size_t column = 0;
    std::size_t nonzeros = 0;
    for (std::size_t i = 0; i < normal.size(); ++i)
    {
        if (normal[i] != 0.0)
        {
            column = i;
            ++nonzeros;
        }
    }
    if (nonzeros > 1)
    {
        return normal.size();
    }
    return column;
}

Constraints gather_constraints(const Problem& problem)
{
    Constraints constraints;
    constraints.n = problem.n;
    constraints.equality_count = problem.equality_rhs.size();
    constraints.rows.reserve(problem.equality_rows.size() + problem.inequality_rows.size());
    constraints.rows.insert(constraints.rows.end(), problem.equality_rows.begin(), problem.equality_rows.end());
    constraints.rows.insert(constraints.rows.end(), problem.inequality_rows.begin(),
                            problem.inequality_rows.end());
    constraints.rhs.reserve(problem.equality_rhs.size() + problem.inequality_rhs.size());
    constraints.rhs.insert(constraints.rhs.end(), problem.equality_rhs.begin(), problem.equality_rhs.end());
    constraints.rhs.insert(constraints.rhs.end(), problem.inequality_rhs.begin(),
                           problem.inequality_rhs.end());
    constraints.explicit_count = constraints.rhs.size();
    constraints.row_norms.resize(constraints.rhs.size());
    constraints.unit_columns.resize(constraints.rhs.size());
    for (std::size_t index = 0; index < constraints.rhs.size(); ++index)
    {
        std::span<const double> normal = row(constraints, index);
        constraints.row_norms[index] = std::sqrt(dot(normal, normal));
        constraints.unit_columns[index] = unit_column_of(normal);
    }
    if (problem.gross_exposure.has_value())
    {
        constraints.family = &*problem.gross_exposure;
    }
    return constraints;
}

State initial_state(const Problem& problem, std::span<const double> factor, ActiveSet set)
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
    state.set = std::move(set);
    std::size_t constraint_count = problem.equality_rhs.size() + problem.inequality_rhs.size();
    if (problem.gross_exposure.has_value())
    {
        // At most n generated rows are active at once, but a generated row can leave and be generated again: the
        // family counts as 2n rows.
        constraint_count += 2 * n;
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
    bool finite_point = std::ranges::all_of(state.x, [](double value) { return std::isfinite(value); });
    for (std::size_t index = constraints.equality_count; index < constraints.explicit_count; ++index)
    {
        double value = scanned_slack(constraints, state, index, finite_point);
        bool within_tolerance = value >= -tolerance(constraints, index, x_scale);
        if (within_tolerance || value >= worst.slack || is_active(state, index) ||
            std::ranges::find(state.passed_over, index) != state.passed_over.end())
        {
            continue;
        }
        worst = Violation{index, value};
    }
    return worst;
}

// s'x is largest over {-1, 1}^n at s = sign(x), taken as +1 where x_i = 0, so (a - s)'x >= b is the generated row
// most violated at x and its slack a'x - ||x||_1 - b is that of the family.
void fill_sign_vector(std::span<const double> x, std::vector<double>& signs)
{
    signs.assign(x.size(), 1.0);
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        if (x[i] < 0.0)
        {
            signs[i] = -1.0;
        }
    }
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

std::size_t append_generated(const Constraints& constraints, GeneratedRows& generated,
                             std::span<const double> normal, std::span<const double> signs)
{
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
        std::vector<double>& signs = constraints.family_signs;
        std::vector<double>& normal = constraints.family_normal;
        fill_sign_vector(state.x, signs);
        normal.assign(constraints.family->offset.begin(), constraints.family->offset.end());
        for (std::size_t i = 0; i < normal.size(); ++i)
        {
            normal[i] -= signs[i];
        }
        double value = dot(normal, state.x) - constraints.family->rhs;
        double norm = std::sqrt(dot(normal, normal));
        double limit = feasibility_tolerance(constraints.family->rhs, norm, max_abs(state.x));
        if (value < -limit && value < worst.slack && !duplicates_active(constraints, state, signs, normal))
        {
            return ActiveConstraint{append_generated(constraints, constraints.generated, normal, signs), 1.0, 0.0};
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

// The directions of a unit row read one entry of each basis column (compute_unit_directions).
void compute_candidate_directions(const Constraints& constraints, State& state,
                                  const ActiveConstraint& candidate)
{
    if (!is_unit(constraints, candidate.index))
    {
        compute_directions(state.set, row(constraints, candidate.index), candidate.sign);
        return;
    }
    compute_unit_directions(state.set, constraints.unit_columns[candidate.index],
                            unit_entry(constraints, candidate.index), candidate.sign);
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

// For a candidate whose normal is the combination dual_direction of the active normals: its slack less the part
// that the residuals of the active rows explain, which is its slack on their face.
double slack_on_active_face(const Constraints& constraints, const State& state, const ActiveConstraint& candidate)
{
    double value = candidate.sign * slack(constraints, state, candidate.index);
    for (std::size_t k = 0; k < state.set.active.size(); ++k)
    {
        const ActiveConstraint& constraint = state.set.active[k];
        value -= state.set.dual_direction[k] * constraint.sign * slack(constraints, state, constraint.index);
    }
    return value;
}

// The face slack b'r - b_c of a candidate on the face of the active rows is known only to the rounding of those
// rows: row k holds at a vertex to about eps max(|b_k|, ||n_k|| |x|max), the rounding of its data and of its slack,
// and the face slack carries it |r_k| times. Where the face is ill-conditioned (a candidate nearly dependent on
// one active row, |r| of order 1e5 or more) the sum exceeds the tolerance of the candidate, and a candidate that
// holds at the vertex reads as violated. kFaceRoundingFactor stands for the few units of rounding of each slack:
// at feasible degenerate vertices the face slack stayed below 0.6 eps sum, and violated candidates above 8e3.
double face_rounding(const Constraints& constraints, const State& state)
{
    double x_scale = max_abs(state.x);
    double sum = 0.0;
    for (std::size_t k = 0; k < state.set.active.size(); ++k)
    {
        std::size_t index = state.set.active[k].index;
        double scale = std::max(std::fabs(rhs_of(constraints, index)), row_norm_of(constraints, index) * x_scale);
        sum += std::fabs(state.set.dual_direction[k]) * scale;
    }
    return kFaceRoundingFactor * std::numeric_limits<double>::epsilon() * sum;
}

// The active rows hold at x only to the rounding of every step from the unconstrained start, which does not shrink
// with x: from a start far larger than the final x it can exceed the tolerance there. A candidate dependent on them
// (the other bound of a fixed weight, a row of the family that negates the budget) is then violated by that
// rounding alone, and a dual step on it drops a row by a ratio over the rounding of dual_direction, or declares the
// problem infeasible. Before any such step its multiplier is zero, so G x + g = N u holds without it.
bool holds_on_active_face(const Constraints& constraints, const State& state, const ActiveConstraint& candidate,
                          double limit)
{
    if (state.dependent != DependentCandidates::PassOverHeldOnFace || candidate.multiplier != 0.0)
    {
        return false;
    }
    double face_slack = slack_on_active_face(constraints, state, candidate);
    double face_limit = limit + face_rounding(constraints, state);
    if (candidate.index < constraints.equality_count)
    {
        return std::fabs(face_slack) <= face_limit;
    }
    return face_slack >= -face_limit;
}

// The candidate is not added. The refinement moves x onto the active rows, and the candidate is not picked again
// until a step changes them: its slack less the part the active residuals explain does not move with x, so on a
// face where the rounding exceeds its tolerance (face_rounding) no refinement brings it within. Each pass counts as
// an iteration, so that passing over the candidates of one face cannot loop forever.
Progress pass_over(const Constraints& constraints, State& state, const ActiveConstraint& candidate)
{
    if (state.iterations == state.iteration_limit)
    {
        return Progress::MaxIterations;
    }
    ++state.iterations;
    refine_active(constraints, state);
    state.passed_over.push_back(candidate.index);
    return Progress::Added;
}

Progress add_constraint(const Constraints& constraints, State& state, ActiveConstraint candidate)
{
    while (true)
    {
        compute_candidate_directions(constraints, state, candidate);
        double violation = slack(constraints, state, candidate.index);
        double limit = tolerance(constraints, candidate.index, max_abs(state.x));
        bool moves_primal = is_independent_of_active(state.set);
        if (!moves_primal && std::fabs(violation) <= limit)
        {
            return Progress::Added;
        }
        if (!moves_primal && holds_on_active_face(constraints, state, candidate, limit))
        {
            return pass_over(constraints, state, candidate);
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
        state.passed_over.clear();
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
    GeneratedRows& previous = constraints.generated;
    GeneratedRows& kept = constraints.renumbered;
    kept.rows.clear();
    kept.signs.clear();
    kept.norms.clear();
    for (ActiveConstraint& constraint : state.set.active)
    {
        if (!is_generated(constraints, constraint.index))
        {
            continue;
        }
        std::size_t k = constraint.index - constraints.explicit_count;
        std::span<const double> normal = std::span<const double>(previous.rows).subspan(k * n, n);
        std::span<const double> signs = std::span<const double>(previous.signs).subspan(k * n, n);
        constraint.index = append_generated(constraints, kept, normal, signs);
    }
    std::swap(previous, kept);
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
    std::ptrdiff_t active_inequalities = std::ranges::count_if(
        state.set.active, [&](const ActiveConstraint& constraint)
        { return constraint.index >= constraints.equality_count && !is_generated(constraints, constraint.index); });
    solution.active_inequalities.reserve(static_cast<std::size_t>(active_inequalities));
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

Solution run_dual_method(const Problem& problem, std::span<const double> factor, ActiveSet set,
                         DependentCandidates dependent)
{
    Constraints constraints = gather_constraints(problem);
    State state = initial_state(problem, factor, std::move(set));
    state.dependent = dependent;
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

// A failed solve is run again from J0, the same start, passing over the candidates that hold on the face of the
// active rows; the iterations of both runs are reported. If that run fails too, the first failure stands.
void retry_passing_over(const Problem& problem, std::span<const double> factor, Solution& solution)
{
    Solution retried = run_dual_method(problem, factor, empty_active_set(factor, problem.n),
                                       DependentCandidates::PassOverHeldOnFace);
    if (retried.status != Status::Optimal)
    {
        return;
    }
    retried.iterations += solution.iterations;
    solution = std::move(retried);
}

// Passing over changes the path of solves that the dual step alone already finishes, and their solutions are
// committed bit for bit (tests/unit/solve/determinism_test.cpp): only a solve that fails takes it.
Solution solve_factored(const Problem& problem, std::span<const double> factor, ActiveSet set)
{
    Solution solution = run_dual_method(problem, factor, std::move(set), DependentCandidates::DualStep);
    if (solution.status != Status::Optimal)
    {
        retry_passing_over(problem, factor, solution);
    }
    return solution;
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
    std::size_t n = problem.n;
    if (!problem.hessian_factor.empty())
    {
        return solve_factored(problem, problem.hessian_factor, empty_active_set(problem.hessian_factor, n));
    }
    std::optional<std::vector<double>> factor = linalg::cholesky(problem.hessian, n);
    if (!factor)
    {
        throw std::invalid_argument("qp: the hessian is not positive definite");
    }
    return solve_factored(problem, *factor, empty_active_set(*factor, n));
}

Solution solve(const Problem& problem, BasisCache& cache)
{
    if (problem.hessian_factor.empty())
    {
        return solve(problem);
    }
    validate_sizes(problem);
    std::span<const double> factor = problem.hessian_factor;
    bool same_factor = cache.factor.size() == factor.size() &&
                       std::memcmp(cache.factor.data(), factor.data(), factor.size_bytes()) == 0;
    if (!same_factor)
    {
        // Build both in locals first: a throw (bad_alloc) must not leave the new key with the old basis.
        std::vector<double> key(factor.begin(), factor.end());
        std::vector<double> basis = initial_basis(factor, problem.n);
        cache.factor = std::move(key);
        cache.basis = std::move(basis);
    }
    return solve_factored(problem, factor, active_set_from(cache.basis, problem.n));
}

} // namespace po::qp
