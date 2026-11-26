#include "riskparity/newton_polish.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "linalg/qr.hpp"
#include "problems/constraints.hpp"
#include "qp/goldfarb_idnani.hpp"

namespace po::riskparity
{
namespace
{

constexpr int kMaxSteps = 20;
// A residual at or below its noise is a KKT point to the rounding of the weights: no step improves on it.
constexpr double kNoiseFloor = 1.0;

enum class Side
{
    Free,
    Lower,
    Upper
};

// A point and the bounds its weights are held at. The bound multiplier of a held weight has the sign of its
// side; that of a free weight is 0.
struct Position
{
    KktPoint point;
    std::vector<Side> sides;
};

// A position with its KKT system and its residual in units of noise.
struct Candidate
{
    Position position;
    KktEvaluation evaluation;
    double residual;
};

// The free weights and the equalities whose gradients over them are independent.
struct Unknowns
{
    std::vector<std::size_t> free;
    std::vector<std::size_t> equalities;
};

// A square system, size x size (row-major), and its right-hand side. The row of each kept equality is
// scaled by its entry of scales, which scales its unknown by the inverse.
struct LinearSystem
{
    std::size_t size;
    std::vector<double> matrix;
    std::vector<double> rhs;
    std::vector<double> scales;
};

std::vector<Side> sides_named_by(std::span<const double> bound_multipliers)
{
    std::vector<Side> sides(bound_multipliers.size(), Side::Free);
    for (std::size_t j = 0; j < sides.size(); ++j)
    {
        if (bound_multipliers[j] > 0.0)
        {
            sides[j] = Side::Lower;
        }
        if (bound_multipliers[j] < 0.0)
        {
            sides[j] = Side::Upper;
        }
    }
    return sides;
}

void hold_at_bounds(const Spec& spec, std::span<const Side> sides, std::vector<double>& weights)
{
    for (std::size_t j = 0; j < weights.size(); ++j)
    {
        if (sides[j] == Side::Lower)
        {
            weights[j] = problems::lower_of(spec, j);
        }
        if (sides[j] == Side::Upper)
        {
            weights[j] = problems::upper_of(spec, j);
        }
    }
}

bool has_sign_of(double multiplier, Side side)
{
    if (side == Side::Lower)
    {
        return multiplier > 0.0;
    }
    return side == Side::Upper && multiplier < 0.0;
}

double norm_of(std::span<const double> values)
{
    double total = 0.0;
    for (double value : values)
    {
        total += value * value;
    }
    return std::sqrt(total);
}

double largest_magnitude(std::span<const double> values)
{
    double largest = 0.0;
    for (double value : values)
    {
        largest = std::fmax(largest, std::fabs(value));
    }
    return largest;
}

std::vector<double> over_free(std::span<const double> row, std::span<const std::size_t> free)
{
    std::vector<double> entries(free.size());
    for (std::size_t a = 0; a < free.size(); ++a)
    {
        entries[a] = row[free[a]];
    }
    return entries;
}

std::vector<double> equality_row(const KktEvaluation& evaluation, std::size_t k,
                                 std::span<const std::size_t> free)
{
    std::size_t n = evaluation.stationarity.size();
    return over_free(std::span<const double>(evaluation.rows).subspan(k * n, n), free);
}

void remove_component(std::span<const double> direction, std::vector<double>& vector)
{
    double projection = 0.0;
    for (std::size_t a = 0; a < vector.size(); ++a)
    {
        projection += direction[a] * vector[a];
    }
    for (std::size_t a = 0; a < vector.size(); ++a)
    {
        vector[a] -= projection * direction[a];
    }
}

// The equalities kept are those whose gradient over the free weights is independent of the earlier ones, by
// the rule of the QP solver: the part outside their span exceeds its relative tolerance. A dependent row
// (equal means with a target return) holds whenever the kept ones do, and its multiplier is left as it is.
void keep_independent_equalities(const KktEvaluation& evaluation, Unknowns& unknowns)
{
    std::vector<std::vector<double>> basis;
    for (std::size_t k = 0; k < evaluation.gaps.size(); ++k)
    {
        std::vector<double> rest = equality_row(evaluation, k, unknowns.free);
        double norm = norm_of(rest);
        for (const std::vector<double>& direction : basis)
        {
            remove_component(direction, rest);
        }
        double outside = norm_of(rest);
        if (!(outside > qp::kRelativeTolerance * norm))
        {
            continue;
        }
        for (double& entry : rest)
        {
            entry /= outside;
        }
        basis.push_back(std::move(rest));
        unknowns.equalities.push_back(k);
    }
}

Unknowns unknowns_of(std::span<const Side> sides, const KktEvaluation& evaluation)
{
    Unknowns unknowns;
    for (std::size_t j = 0; j < sides.size(); ++j)
    {
        if (sides[j] == Side::Free)
        {
            unknowns.free.push_back(j);
        }
    }
    keep_independent_equalities(evaluation, unknowns);
    return unknowns;
}

// The multipliers of the kept equalities that leave the smallest stationarity residual over the free weights,
//   min over lambda of |(grad f - sum_k lambda_k grad h_k) over the free weights|,
// by a Householder QR of the gradients, each scaled to a largest entry of 1. The stationarity of evaluation
// follows the new multipliers. Nothing changes when the gradients are dependent to working precision.
void estimate_multipliers(const Unknowns& unknowns, KktEvaluation& evaluation, KktPoint& point)
{
    std::size_t rows = unknowns.free.size();
    std::size_t cols = unknowns.equalities.size();
    if (cols == 0 || rows < cols)
    {
        return;
    }
    std::vector<double> gradients(rows * cols);
    std::vector<double> scales(cols);
    for (std::size_t b = 0; b < cols; ++b)
    {
        std::vector<double> row = equality_row(evaluation, unknowns.equalities[b], unknowns.free);
        scales[b] = largest_magnitude(row);
        for (std::size_t a = 0; a < rows; ++a)
        {
            gradients[a * cols + b] = row[a] / scales[b];
        }
    }
    linalg::Qr qr = linalg::householder_qr(gradients, rows, cols);
    if (!linalg::has_full_column_rank(qr))
    {
        return;
    }
    std::vector<double> change = linalg::least_squares(qr, over_free(evaluation.stationarity, unknowns.free));
    std::size_t n = evaluation.stationarity.size();
    for (std::size_t b = 0; b < cols; ++b)
    {
        std::size_t k = unknowns.equalities[b];
        double multiplier_change = change[b] / scales[b];
        point.equality_multipliers[k] += multiplier_change;
        for (std::size_t j = 0; j < n; ++j)
        {
            evaluation.stationarity[j] -= multiplier_change * evaluation.rows[k * n + j];
        }
    }
}

// The multipliers at the weights of the position. Those of the equalities are the least-squares ones over
// the free weights. That of a held weight is what its stationarity row leaves, beta_j = (grad f - sum_k
// lambda_k grad h_k)_j; one with the sign of the other side says that f decreases into the box, so the bound
// is released and the row keeps its residual.
Candidate settled(const KktProblem& problem, Position position)
{
    KktPoint& point = position.point;
    std::ranges::fill(point.bound_multipliers, 0.0);
    KktEvaluation evaluation = evaluate_kkt(problem, point);
    estimate_multipliers(unknowns_of(position.sides, evaluation), evaluation, point);
    for (std::size_t j = 0; j < position.sides.size(); ++j)
    {
        if (!has_sign_of(evaluation.stationarity[j], position.sides[j]))
        {
            position.sides[j] = Side::Free;
            continue;
        }
        point.bound_multipliers[j] = evaluation.stationarity[j];
    }
    update_multipliers(problem, point, evaluation);
    double residual = largest(residual_of(problem, point, evaluation));
    return Candidate{std::move(position), std::move(evaluation), residual};
}

Candidate held_start(const KktProblem& problem, const KktPoint& start)
{
    Position position{start, sides_named_by(start.bound_multipliers)};
    hold_at_bounds(problem.spec, position.sides, position.point.weights);
    return settled(problem, std::move(position));
}

// The Newton system of the KKT conditions over the unknowns,
//   [ H  -A' ] [ dw      ]   [ -(grad f - A'lambda) ]
//   [ A   0  ] [ dlambda ] = [ -h                   ],
// H the Hessian of the Lagrangian over the free weights and A the gradients of the kept equalities. Each row
// of A and its right-hand side are scaled to the size of H, which scales dlambda by the inverse.
LinearSystem newton_system(const KktEvaluation& evaluation, const Unknowns& unknowns)
{
    std::size_t n = evaluation.stationarity.size();
    std::size_t free_count = unknowns.free.size();
    std::size_t size = free_count + unknowns.equalities.size();
    LinearSystem system{size, std::vector<double>(size * size, 0.0), std::vector<double>(size),
                        std::vector<double>(unknowns.equalities.size())};
    for (std::size_t a = 0; a < free_count; ++a)
    {
        std::size_t j = unknowns.free[a];
        for (std::size_t b = 0; b < free_count; ++b)
        {
            system.matrix[a * size + b] = evaluation.lagrangian_hessian[j * n + unknowns.free[b]];
        }
        system.rhs[a] = -evaluation.stationarity[j];
    }
    double hessian_scale = largest_magnitude(system.matrix);
    for (std::size_t b = 0; b < unknowns.equalities.size(); ++b)
    {
        std::size_t k = unknowns.equalities[b];
        std::vector<double> row = equality_row(evaluation, k, unknowns.free);
        double scale = hessian_scale / largest_magnitude(row);
        system.scales[b] = scale;
        for (std::size_t a = 0; a < free_count; ++a)
        {
            system.matrix[a * size + free_count + b] = -scale * row[a];
            system.matrix[(free_count + b) * size + a] = scale * row[a];
        }
        system.rhs[free_count + b] = -scale * evaluation.gaps[k];
    }
    return system;
}

// What solution_of does with a matrix that fails the rank test of the QR (linalg::has_full_column_rank).
enum class RankTest
{
    // No solution: the step of the weights.
    Required,
    // The solution all the same, when finite: the correction of the multipliers, which has its own tests.
    Ignored
};

bool all_finite(std::span<const double> values)
{
    return std::ranges::all_of(values, [](double value) { return std::isfinite(value); });
}

// The solution of the Newton system by a Householder QR of its matrix; none when there are no free weights,
// when an entry of the solution is not finite, or, under RankTest::Required, when the matrix is rank
// deficient to working precision: the active gradients are dependent or the Hessian is singular along the
// constraints.
std::optional<std::vector<double>> solution_of(LinearSystem system, const Unknowns& unknowns,
                                               RankTest rank_test)
{
    std::optional<std::vector<double>> solution;
    double scale = largest_magnitude(system.matrix);
    if (unknowns.free.empty() || scale == 0.0 || !std::isfinite(scale))
    {
        return solution;
    }
    linalg::Qr qr = linalg::householder_qr(system.matrix, system.size, system.size);
    if (rank_test == RankTest::Required && !linalg::has_full_column_rank(qr))
    {
        return solution;
    }
    std::vector<double> values = linalg::least_squares(qr, std::move(system.rhs));
    if (!all_finite(values))
    {
        return solution;
    }
    solution.emplace(std::move(values));
    return solution;
}

// A free weight that the step takes out of the box stops at the bound it crossed and is held there.
Position moved_by(const Spec& spec, const Candidate& candidate, const Unknowns& unknowns,
                  std::span<const double> step)
{
    Position position = candidate.position;
    for (std::size_t a = 0; a < unknowns.free.size(); ++a)
    {
        std::size_t j = unknowns.free[a];
        double moved = position.point.weights[j] + step[a];
        if (moved < problems::lower_of(spec, j))
        {
            position.sides[j] = Side::Lower;
        }
        if (moved > problems::upper_of(spec, j))
        {
            position.sides[j] = Side::Upper;
        }
        position.point.weights[j] = moved;
    }
    hold_at_bounds(spec, position.sides, position.point.weights);
    return position;
}

// The step of the weights is that of the Newton system (solution_of); the multipliers are estimated again at
// the new weights, which keeps the step a good one where the gradients of the constraints turn quickly (a
// target volatility just above the minimum). No step where the system has no solution.
std::optional<Candidate> newton_step(const KktProblem& problem, const Candidate& candidate)
{
    std::optional<Candidate> next;
    Unknowns unknowns = unknowns_of(candidate.position.sides, candidate.evaluation);
    std::optional<std::vector<double>> step =
        solution_of(newton_system(candidate.evaluation, unknowns), unknowns, RankTest::Required);
    if (!step.has_value())
    {
        return next;
    }
    next.emplace(settled(problem, moved_by(problem.spec, candidate, unknowns, *step)));
    return next;
}

// The step of the free weights and the change of the multipliers of the kept equalities that the Newton
// system at a point gives. The rank test of the QR does not veto it: at cond(Sigma) = 1e8 with n >= 12 the
// scaled KKT matrix failed it on about half the shapes measured (before the portable generators), while the
// solution corrected the multipliers to 3.8e-10 or better there (docs/guide/validation.md); contraction,
// stationarity and sides decide (corrected_multipliers).
struct NewtonCorrection
{
    std::vector<double> weights;
    std::vector<double> multipliers;
};

std::optional<NewtonCorrection> newton_correction(const KktEvaluation& evaluation, const Unknowns& unknowns)
{
    std::optional<NewtonCorrection> correction;
    LinearSystem system = newton_system(evaluation, unknowns);
    std::vector<double> scales = system.scales;
    std::optional<std::vector<double>> solution =
        solution_of(std::move(system), unknowns, RankTest::Ignored);
    if (!solution.has_value())
    {
        return correction;
    }
    auto split = solution->begin() + static_cast<std::ptrdiff_t>(unknowns.free.size());
    correction.emplace(NewtonCorrection{std::vector<double>(solution->begin(), split),
                                        std::vector<double>(split, solution->end())});
    for (std::size_t b = 0; b < scales.size(); ++b)
    {
        correction->multipliers[b] *= scales[b];
    }
    return correction;
}

// beta_j = (grad f - sum_k lambda_k grad h_k)_j + (H dw)_j of held weight j, the lambda of the point moved by
// the correction: what its row of the Newton system leaves.
double held_multiplier(const KktEvaluation& evaluation, const KktPoint& point, const Unknowns& unknowns,
                       const NewtonCorrection& correction, std::size_t j)
{
    std::size_t n = point.weights.size();
    double beta = evaluation.stationarity[j] + point.bound_multipliers[j];
    for (std::size_t b = 0; b < unknowns.equalities.size(); ++b)
    {
        beta -= correction.multipliers[b] * evaluation.rows[unknowns.equalities[b] * n + j];
    }
    for (std::size_t a = 0; a < unknowns.free.size(); ++a)
    {
        beta += evaluation.lagrangian_hessian[j * n + unknowns.free[a]] * correction.weights[a];
    }
    return beta;
}

KktPoint corrected_point(const KktPoint& point, const KktEvaluation& evaluation, const Unknowns& unknowns,
                         const NewtonCorrection& correction)
{
    KktPoint corrected = point;
    for (std::size_t b = 0; b < unknowns.equalities.size(); ++b)
    {
        corrected.equality_multipliers[unknowns.equalities[b]] += correction.multipliers[b];
    }
    for (std::size_t j = 0; j < point.weights.size(); ++j)
    {
        if (point.bound_multipliers[j] != 0.0)
        {
            corrected.bound_multipliers[j] = held_multiplier(evaluation, point, unknowns, correction, j);
        }
    }
    return corrected;
}

bool keeps_sides(std::span<const double> bound_multipliers, std::span<const Side> sides)
{
    for (std::size_t j = 0; j < sides.size(); ++j)
    {
        if (sides[j] != Side::Free && !has_sign_of(bound_multipliers[j], sides[j]))
        {
            return false;
        }
    }
    return true;
}

// A NaN residual is worse than any number.
bool improves_on(double residual, double current)
{
    return residual < current || (std::isnan(current) && !std::isnan(residual));
}

} // namespace

// The second solve, from the corrected multipliers at the same weights, measures how far the first is from
// its limit. The weights do not move, so the gradients of the equalities and the unknowns are those of the
// first solve, and the evaluation changes only with the multipliers.
MultiplierCorrection corrected_multipliers(const KktProblem& problem, const KktPoint& point,
                                           KktEvaluation evaluation)
{
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    MultiplierCorrection result{point, kNaN, kNaN, false};
    std::vector<Side> sides = sides_named_by(point.bound_multipliers);
    Unknowns unknowns = unknowns_of(sides, evaluation);
    std::optional<NewtonCorrection> first = newton_correction(evaluation, unknowns);
    if (!first.has_value())
    {
        return result;
    }
    result.change = norm_of(first->multipliers);
    KktPoint corrected = corrected_point(point, evaluation, unknowns, *first);
    update_multipliers(problem, corrected, evaluation);
    std::optional<NewtonCorrection> second = newton_correction(evaluation, unknowns);
    if (!second.has_value())
    {
        return result;
    }
    result.second_change = norm_of(second->multipliers);
    bool contracts = result.second_change <= 0.5 * result.change;
    // The factor of the KKT gate, measured, not derived: near sigma_min a correction that contracts can still
    // leave 958 times the noise and multipliers ten times worse than least squares, while away from it every
    // correction measured left at most 0.38 (docs/guide/validation.md, "The Newton polish"). Near sigma_min
    // the factor does not separate good corrections from bad ones; it only refuses the worst.
    bool stationary = residual_of(problem, corrected, evaluation).stationarity <= kKktNoiseFactor;
    if (!contracts || !stationary || !keeps_sides(corrected.bound_multipliers, sides))
    {
        return result;
    }
    result.reported_point = std::move(corrected);
    result.accepted = true;
    return result;
}

PolishResult newton_polish(const KktProblem& problem, const KktPoint& start)
{
    Candidate best = held_start(problem, start);
    int steps = 0;
    for (int step = 1; step <= kMaxSteps && !(best.residual <= kNoiseFloor); ++step)
    {
        std::optional<Candidate> next = newton_step(problem, best);
        if (!next.has_value() || !improves_on(next->residual, best.residual))
        {
            break;
        }
        best = std::move(*next);
        steps = step;
    }
    return PolishResult{std::move(best.position.point), std::move(best.evaluation), best.residual, steps};
}

} // namespace po::riskparity
