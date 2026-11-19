#include "problems/parametric_path.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "linalg/cholesky.hpp"
#include "linalg/qr.hpp"
#include "problems/mean_variance.hpp"
#include "qp/goldfarb_idnani.hpp"

namespace po::problems
{
namespace
{

// The rounding of Q(t) and of the slope, as a multiple of the noise models below; the margins it leaves on the
// limit certificate are measured in docs/guide/validation.md.
constexpr double kNoiseFactor = 16.0;
constexpr double kInfinity = std::numeric_limits<double>::infinity();
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kEpsilon = std::numeric_limits<double>::epsilon();

double dot(std::span<const double> a, std::span<const double> b)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        sum += a[i] * b[i];
    }
    return sum;
}

double largest_magnitude(std::span<const double> values)
{
    double largest = 0.0;
    for (double value : values)
    {
        largest = std::max(largest, std::fabs(value));
    }
    return largest;
}

std::vector<double> difference(std::span<const double> later, std::span<const double> earlier)
{
    std::vector<double> change(later.size());
    for (std::size_t i = 0; i < later.size(); ++i)
    {
        change[i] = later[i] - earlier[i];
    }
    return change;
}

// L'w, so that w'Sigma w = |L'w|^2.
std::vector<double> factor_transposed_times(std::span<const double> factor, std::span<const double> weights)
{
    std::size_t n = weights.size();
    std::vector<double> product(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j <= i; ++j)
        {
            product[j] += factor[i * n + j] * weights[i];
        }
    }
    return product;
}

std::vector<double> excess_of(const Data& data)
{
    double rate = data.risk_free().value_or(0.0);
    std::vector<double> excess(data.mu().begin(), data.mu().end());
    for (double& value : excess)
    {
        value -= rate;
    }
    return excess;
}

// The normal r of the face (a - s for a generated row).
std::vector<double> face_normal(const qp::Problem& problem, const Face& face)
{
    std::size_t n = problem.n;
    std::vector<double> normal(n);
    if (face.kind == Face::Kind::Equality)
    {
        std::copy_n(problem.equality_rows.begin() + static_cast<std::ptrdiff_t>(face.row * n), n, normal.begin());
        return normal;
    }
    if (face.kind == Face::Kind::Inequality)
    {
        std::copy_n(problem.inequality_rows.begin() + static_cast<std::ptrdiff_t>(face.row * n), n, normal.begin());
        return normal;
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        normal[i] = problem.gross_exposure->offset[i] - face.signs[i];
    }
    return normal;
}

// The right-hand side of a face held as an equality (an inequality or a generated row).
double face_rhs(const qp::Problem& problem, const Face& face)
{
    if (face.kind == Face::Kind::Inequality)
    {
        return problem.inequality_rhs[face.row];
    }
    return problem.gross_exposure->rhs;
}

// v = mu~ - amount r with amount the midrange of r_i mu~_i, which minimizes |v|max for a normal of entries -1 and 1.
Shift shift_along(const qp::Problem& problem, Face face, std::span<const double> excess)
{
    std::vector<double> normal = face_normal(problem, face);
    std::vector<double> loads(excess.size());
    for (std::size_t i = 0; i < excess.size(); ++i)
    {
        loads[i] = normal[i] * excess[i];
    }
    auto [low, high] = std::ranges::minmax(loads);
    double amount = 0.5 * (low + high);
    std::vector<double> shifted(excess.begin(), excess.end());
    for (std::size_t i = 0; i < shifted.size(); ++i)
    {
        shifted[i] -= amount * normal[i];
    }
    return Shift{std::move(face), amount, std::move(shifted)};
}

// The budget where the rows hold 1'w = 1, whose shift is always exact; Q(t) itself otherwise.
Shift base_shift(const Translation& translation, std::span<const double> excess)
{
    if (!translation.rows.budget.has_value())
    {
        return Shift{Face{}, 0.0, std::vector<double>(excess.begin(), excess.end())};
    }
    return shift_along(translation.problem, Face{Face::Kind::Equality, *translation.rows.budget, {}}, excess);
}

bool unit_normal(std::span<const double> normal)
{
    return std::ranges::all_of(normal, [](double entry) { return std::fabs(entry) == 1.0; });
}

bool held_as_equality(const Face& face)
{
    return face.kind == Face::Kind::Inequality || face.kind == Face::Kind::Generated;
}

bool same_face(const Face& first, const Face& second)
{
    return first.kind == second.kind && first.row == second.row && first.signs == second.signs;
}

// The position of signs among the k sign vectors of solution (lexicographic order), and whether it is there.
std::pair<std::size_t, bool> sign_position(const qp::Solution& solution, std::span<const double> signs)
{
    std::size_t n = signs.size();
    std::size_t k = solution.gross_exposure_multipliers.size();
    std::span<const double> all(solution.gross_exposure_signs);
    for (std::size_t j = 0; j < k; ++j)
    {
        std::span<const double> current = all.subspan(j * n, n);
        if (std::ranges::equal(current, signs))
        {
            return {j, true};
        }
        if (std::ranges::lexicographical_compare(signs, current))
        {
            return {j, false};
        }
    }
    return {k, false};
}

// Adds multiplier to the face row of solution, entering it in the active set where it is not.
void add_to_face(qp::Solution& solution, const Face& face, double multiplier)
{
    if (face.kind == Face::Kind::Inequality)
    {
        solution.inequality_multipliers[face.row] += multiplier;
        auto place = std::ranges::lower_bound(solution.active_inequalities, face.row);
        if (place == solution.active_inequalities.end() || *place != face.row)
        {
            solution.active_inequalities.insert(place, face.row);
        }
        return;
    }
    auto [position, present] = sign_position(solution, face.signs);
    if (present)
    {
        solution.gross_exposure_multipliers[position] += multiplier;
        return;
    }
    std::size_t n = face.signs.size();
    solution.gross_exposure_signs.insert(solution.gross_exposure_signs.begin() +
                                             static_cast<std::ptrdiff_t>(position * n),
                                         face.signs.begin(), face.signs.end());
    solution.gross_exposure_multipliers.insert(
        solution.gross_exposure_multipliers.begin() + static_cast<std::ptrdiff_t>(position), multiplier);
}

// The multipliers of Q(t) from those of the shifted Q(t): Sigma w - t mu~ = N u_v - t amount r, so the face
// multiplier is u_v - t amount; a face held as an equality leaves the equality rows and joins the active rows.
void to_unshifted(qp::Solution& solution, const Shift& shift, double t)
{
    if (shift.face.kind == Face::Kind::Equality)
    {
        solution.equality_multipliers[shift.face.row] -= t * shift.amount;
        return;
    }
    if (!held_as_equality(shift.face))
    {
        return;
    }
    double multiplier = solution.equality_multipliers.back() - t * shift.amount;
    solution.equality_multipliers.pop_back();
    add_to_face(solution, shift.face, multiplier);
}

// The face as an equality row: a generated face is appended to the equality rows, and an explicit inequality face
// moves there, so that the solver never holds the same row twice.
void hold_face(qp::Problem& problem, const Face& face)
{
    std::vector<double> normal = face_normal(problem, face);
    double rhs = face_rhs(problem, face);
    if (face.kind == Face::Kind::Inequality)
    {
        auto row = problem.inequality_rows.begin() + static_cast<std::ptrdiff_t>(face.row * problem.n);
        problem.inequality_rows.erase(row, row + static_cast<std::ptrdiff_t>(problem.n));
        problem.inequality_rhs.erase(problem.inequality_rhs.begin() + static_cast<std::ptrdiff_t>(face.row));
    }
    problem.equality_rows.insert(problem.equality_rows.end(), normal.begin(), normal.end());
    problem.equality_rhs.push_back(rhs);
}

// Undoes hold_face on the problem, and on the solution the removal of an explicit inequality face (its multiplier
// is set by to_unshifted).
void release_face(qp::Problem& problem, const Face& face, qp::Solution& solution)
{
    std::size_t n = problem.n;
    auto normal = problem.equality_rows.end() - static_cast<std::ptrdiff_t>(n);
    if (face.kind == Face::Kind::Inequality)
    {
        problem.inequality_rows.insert(problem.inequality_rows.begin() + static_cast<std::ptrdiff_t>(face.row * n),
                                       normal, problem.equality_rows.end());
        problem.inequality_rhs.insert(problem.inequality_rhs.begin() + static_cast<std::ptrdiff_t>(face.row),
                                      problem.equality_rhs.back());
        double filler = solution.status == qp::Status::Optimal ? 0.0 : kNaN;
        solution.inequality_multipliers.insert(
            solution.inequality_multipliers.begin() + static_cast<std::ptrdiff_t>(face.row), filler);
        for (std::size_t& row : solution.active_inequalities)
        {
            row += row >= face.row ? 1 : 0;
        }
    }
    problem.equality_rows.resize(problem.equality_rows.size() - n);
    problem.equality_rhs.pop_back();
}

PathPoint solve_with(Path& path, const Shift& shift, double t)
{
    qp::Problem& problem = path.translation.problem;
    for (std::size_t i = 0; i < problem.n; ++i)
    {
        problem.linear[i] = -t * shift.excess[i];
    }
    if (held_as_equality(shift.face))
    {
        hold_face(problem, shift.face);
    }
    PathPoint point{t, qp::solve(problem, path.basis), kNaN, shift};
    if (held_as_equality(shift.face))
    {
        release_face(problem, shift.face, point.solution);
    }
    path.iterations += point.solution.iterations;
    if (point.solution.status != qp::Status::Optimal)
    {
        return point;
    }
    to_unshifted(point.solution, shift, t);
    std::vector<double> root = factor_transposed_times(path.factor, point.solution.x);
    point.volatility = std::sqrt(dot(root, root));
    return point;
}

// The face multiplier of Q(t) of a point solved with its face held as an equality.
double face_multiplier(const PathPoint& point)
{
    const Face& face = point.shift.face;
    if (face.kind == Face::Kind::Inequality)
    {
        return point.solution.inequality_multipliers[face.row];
    }
    return point.solution.gross_exposure_multipliers[sign_position(point.solution, face.signs).first];
}

// The optimum restricted to the face solves Q(t) when the face multiplier of Q(t) is nonnegative: then the KKT
// conditions of Q(t) hold with the face among the active rows.
bool certified(const PathPoint& point)
{
    return point.solution.status == qp::Status::Optimal && face_multiplier(point) >= 0.0;
}

// The active inequality face of largest multiplier: an explicit row or a generated row whose normal has every
// entry -1 or 1.
std::optional<Face> loaded_face(const Path& path, const PathPoint& point)
{
    const qp::Problem& problem = path.translation.problem;
    const qp::Solution& solution = point.solution;
    std::optional<Face> face;
    double largest = 0.0;
    for (std::size_t row : solution.active_inequalities)
    {
        Face candidate{Face::Kind::Inequality, row, {}};
        if (solution.inequality_multipliers[row] > largest && unit_normal(face_normal(problem, candidate)))
        {
            largest = solution.inequality_multipliers[row];
            face = std::move(candidate);
        }
    }
    std::size_t n = problem.n;
    for (std::size_t j = 0; j < solution.gross_exposure_multipliers.size(); ++j)
    {
        auto signs = solution.gross_exposure_signs.begin() + static_cast<std::ptrdiff_t>(j * n);
        Face candidate{Face::Kind::Generated, 0, std::vector<double>(signs, signs + static_cast<std::ptrdiff_t>(n))};
        if (solution.gross_exposure_multipliers[j] > largest && unit_normal(face_normal(problem, candidate)))
        {
            largest = solution.gross_exposure_multipliers[j];
            face = std::move(candidate);
        }
    }
    return face;
}

// A weight of Q(t) is known to n eps (max(1, |w|max) + t |v|max): the noise of the stationarity equation
// Sigma w - t v = N u of the shifted Q(t) solved, which does not carry cond(Sigma).
double weight_noise(const PathPoint& earlier, const PathPoint& later)
{
    double n = static_cast<double>(later.solution.x.size());
    double scale =
        std::max({1.0, largest_magnitude(earlier.solution.x), largest_magnitude(later.solution.x)});
    return kNoiseFactor * n * kEpsilon * (scale + later.t * largest_magnitude(later.shift.excess));
}

// b_v from N b_v = -v by Householder least squares is known to n eps |v|max.
double slope_noise(const PathPoint& point)
{
    const std::vector<double>& excess = point.shift.excess;
    return kNoiseFactor * static_cast<double>(excess.size()) * kEpsilon * largest_magnitude(excess);
}

std::size_t generated_count(const PathPoint& point)
{
    return point.solution.gross_exposure_multipliers.size();
}

std::size_t active_count(const qp::Problem& problem, const PathPoint& point)
{
    return problem.equality_rhs.size() + point.solution.active_inequalities.size() + generated_count(point);
}

// The solver keeps the active normals independent, the generated ones among them, so the equality rows and the
// active rows have rank equal to their count, and a count of n fixes w. (MAX_RETURN has one equality row at most,
// the budget.)
bool is_vertex(const Path& path, const PathPoint& point)
{
    const qp::Problem& problem = path.translation.problem;
    return active_count(problem, point) == problem.n;
}

// The active generated rows a - s_j (k x n, row-major), in the order of the sign vectors of the solution.
std::vector<double> generated_rows(const qp::Problem& problem, const PathPoint& point)
{
    std::vector<double> rows = point.solution.gross_exposure_signs;
    std::size_t n = problem.n;
    for (std::size_t k = 0; k < generated_count(point); ++k)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            rows[k * n + i] = problem.gross_exposure->offset[i] - rows[k * n + i];
        }
    }
    return rows;
}

// The rows (row-major, n columns) as the columns first, first + 1, ... of normals (n x m, row-major).
void place_columns(std::span<const double> rows, std::size_t n, std::size_t first, std::size_t m,
                   std::vector<double>& normals)
{
    for (std::size_t k = 0; k < rows.size() / n; ++k)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            normals[i * m + first + k] = rows[k * n + i];
        }
    }
}

// N (n x m, row-major), the active normals as columns: the equality rows first, then the active rows in ascending
// order, then the active generated rows in the order of their sign vectors.
std::vector<double> active_normals(const qp::Problem& problem, const PathPoint& point)
{
    std::size_t n = problem.n;
    std::size_t equalities = problem.equality_rhs.size();
    std::size_t inequalities = point.solution.active_inequalities.size();
    std::size_t m = active_count(problem, point);
    std::vector<double> normals(n * m);
    place_columns(problem.equality_rows, n, 0, m, normals);
    for (std::size_t k = 0; k < inequalities; ++k)
    {
        std::size_t row = point.solution.active_inequalities[k];
        place_columns(std::span<const double>(problem.inequality_rows).subspan(row * n, n), n, equalities + k, m,
                      normals);
    }
    place_columns(generated_rows(problem, point), n, equalities + inequalities, m, normals);
    return normals;
}

double row_tolerance(std::span<const double> row, double rhs, double x_scale)
{
    return qp::feasibility_tolerance(rhs, std::sqrt(dot(row, row)), x_scale);
}

// a'x - ||x||_1 >= b, decided as the QP solver decides it: by the rule of its generated row at s = sign(x)
// (s_i = 1 where x_i = 0), the most violated row of the family.
bool family_holds(const qp::GrossExposure& family, std::span<const double> x)
{
    std::vector<double> normal = family.offset;
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        normal[i] -= x[i] < 0.0 ? -1.0 : 1.0;
    }
    double slack = dot(normal, x) - family.rhs;
    return slack >= -row_tolerance(normal, family.rhs, largest_magnitude(x));
}

bool rows_hold(const qp::Problem& problem, std::span<const double> x)
{
    std::size_t n = problem.n;
    double x_scale = largest_magnitude(x);
    for (std::size_t k = 0; k < problem.equality_rhs.size(); ++k)
    {
        std::span<const double> row = std::span<const double>(problem.equality_rows).subspan(k * n, n);
        double rhs = problem.equality_rhs[k];
        if (std::fabs(dot(row, x) - rhs) > row_tolerance(row, rhs, x_scale))
        {
            return false;
        }
    }
    for (std::size_t k = 0; k < problem.inequality_rhs.size(); ++k)
    {
        std::span<const double> row = std::span<const double>(problem.inequality_rows).subspan(k * n, n);
        double rhs = problem.inequality_rhs[k];
        if (dot(row, x) - rhs < -row_tolerance(row, rhs, x_scale))
        {
            return false;
        }
    }
    return !problem.gross_exposure.has_value() || family_holds(*problem.gross_exposure, x);
}

// The rows (row-major, n columns) transposed times their multipliers: their term in the stationarity equation.
std::vector<double> transposed_times(std::span<const double> rows, std::span<const double> multipliers,
                                     std::size_t n)
{
    std::vector<double> product(n, 0.0);
    for (std::size_t k = 0; k < multipliers.size(); ++k)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            product[i] += rows[k * n + i] * multipliers[k];
        }
    }
    return product;
}

// The multipliers of Q(t) from those of the shifted Q(t), or their slope b from b_v: amount (t times the shift amount,
// or the shift amount for a slope) less at the face of point, the face being among its active rows.
Slope unshifted(Slope slope, const PathPoint& point, double amount)
{
    const Face& face = point.shift.face;
    if (face.kind == Face::Kind::Equality)
    {
        slope.equality[face.row] -= amount;
        return slope;
    }
    if (face.kind == Face::Kind::Inequality)
    {
        slope.inequality[face.row] -= amount;
        return slope;
    }
    if (face.kind == Face::Kind::Generated)
    {
        slope.gross_exposure[sign_position(point.solution, face.signs).first] -= amount;
    }
    return slope;
}

// |v + E'b_E + C'b_C + sum_j b_j (a - s_j)|max relative to the largest of its four terms, within the rounding of
// b, in the frame of the point's shift (v and b_v).
bool slope_stationary(const Path& path, const PathPoint& point, const Slope& slope)
{
    const qp::Problem& problem = path.translation.problem;
    const std::vector<double>& excess = point.shift.excess;
    std::size_t n = problem.n;
    std::vector<double> equality = transposed_times(problem.equality_rows, slope.equality, n);
    std::vector<double> inequality = transposed_times(problem.inequality_rows, slope.inequality, n);
    std::vector<double> generated = transposed_times(generated_rows(problem, point), slope.gross_exposure, n);
    std::vector<double> residual(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        residual[i] = excess[i] + equality[i] + inequality[i] + generated[i];
    }
    double scale = std::max({largest_magnitude(excess), largest_magnitude(equality), largest_magnitude(inequality),
                             largest_magnitude(generated)});
    return largest_magnitude(residual) <= kNoiseFactor * static_cast<double>(n) * kEpsilon * scale;
}

// The t where a multiplier of value multiplier at point.t and slope rate reaches zero; infinity unless it
// decreases, a slope below -slope_noise.
double zero_crossing(const PathPoint& point, double multiplier, double rate)
{
    if (!(rate < -slope_noise(point)))
    {
        return kInfinity;
    }
    return point.t - multiplier / rate;
}

// Q(t) solved on the face of the path, when there is one and the solution is certified.
std::optional<PathPoint> on_path_face(Path& path, double t)
{
    std::optional<PathPoint> point;
    if (path.face.face.kind == Face::Kind::None)
    {
        return point;
    }
    point = solve_with(path, path.face, t);
    if (!certified(*point))
    {
        point.reset();
    }
    return point;
}

// Without a budget, Q(t) solved on the most loaded active inequality face of point (solved with the base shift),
// when it is not the face of the path, already tried, and the solution is certified; that face becomes the path's.
std::optional<PathPoint> on_loaded_face(Path& path, const PathPoint& point)
{
    std::optional<PathPoint> shifted;
    if (path.base.face.kind != Face::Kind::None || point.solution.status != qp::Status::Optimal)
    {
        return shifted;
    }
    std::optional<Face> face = loaded_face(path, point);
    if (!face.has_value() || same_face(*face, path.face.face))
    {
        return shifted;
    }
    Shift shift = shift_along(path.translation.problem, std::move(*face), path.excess);
    shifted = solve_with(path, shift, point.t);
    if (!certified(*shifted))
    {
        shifted.reset();
        return shifted;
    }
    path.face = std::move(shift);
    return shifted;
}

// An inequality multiplier is nonnegative. The combinations taken here are nonnegative in exact arithmetic; a
// negative value is rounding: a slope above -slope_noise (not a decrease), that slope times the distance
// along it, or the rounding of u(point) at the exit itself. It is set to 0.
double nonnegative(double multiplier)
{
    return std::max(multiplier, 0.0);
}

// The multipliers u of N u = target over the active normals of point, by least squares, placed by row: the
// equality rows, every inequality row (0 where inactive) and the active generated rows in the order of the point's
// sign vectors.
Slope active_multipliers(const qp::Problem& problem, const PathPoint& point, std::vector<double> target)
{
    std::size_t equalities = problem.equality_rhs.size();
    std::size_t inequalities = point.solution.active_inequalities.size();
    std::vector<double> b = linalg::least_squares(
        linalg::householder_qr(active_normals(problem, point), problem.n, active_count(problem, point)),
        std::move(target));
    auto generated = b.begin() + static_cast<std::ptrdiff_t>(equalities + inequalities);
    Slope slope{std::vector<double>(b.begin(), b.begin() + static_cast<std::ptrdiff_t>(equalities)),
                std::vector<double>(problem.inequality_rhs.size(), 0.0), std::vector<double>(generated, b.end())};
    for (std::size_t k = 0; k < inequalities; ++k)
    {
        slope.inequality[point.solution.active_inequalities[k]] = b[equalities + k];
    }
    return slope;
}

double sign_of(double value)
{
    if (value > 0.0)
    {
        return 1.0;
    }
    if (value < 0.0)
    {
        return -1.0;
    }
    return 0.0;
}

// The entry the active sign vectors share, 0 where two of them disagree: there the face holds w_i = 0 (the zeros D).
// Empty when no generated row is active.
std::vector<double> common_signs(const PathPoint& point)
{
    std::size_t n = point.solution.x.size();
    std::span<const double> all(point.solution.gross_exposure_signs);
    std::vector<double> signs;
    if (generated_count(point) == 0)
    {
        return signs;
    }
    signs.assign(all.begin(), all.begin() + static_cast<std::ptrdiff_t>(n));
    for (std::size_t j = 1; j < generated_count(point); ++j)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            if (all[j * n + i] != signs[i])
            {
                signs[i] = 0.0;
            }
        }
    }
    return signs;
}

// The family's part of the face of point: sign(w_i) on the support, 0 on D and where w_i = 0; empty when no
// generated row is active. It does not depend on which sign vectors hold the face.
std::vector<double> support_signs(const PathPoint& point)
{
    std::vector<double> signs = common_signs(point);
    for (std::size_t i = 0; i < signs.size(); ++i)
    {
        if (signs[i] != 0.0)
        {
            signs[i] = sign_of(point.solution.x[i]);
        }
    }
    return signs;
}

// L w, so that Sigma w = L (L'w).
std::vector<double> factor_times(std::span<const double> factor, std::span<const double> vector)
{
    std::size_t n = vector.size();
    std::vector<double> product(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j <= i; ++j)
        {
            product[i] += factor[i * n + j] * vector[j];
        }
    }
    return product;
}

// |Sigma|max, the largest diagonal entry |L_i|^2.
double largest_variance(std::span<const double> factor, std::size_t n)
{
    double largest = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        std::span<const double> row = factor.subspan(i * n, i + 1);
        largest = std::max(largest, dot(row, row));
    }
    return largest;
}

// The rows of the face of point over its free weights F as the columns of normals (F.size() x rhs.size(),
// row-major): the equality rows, the active rows that are not bounds and, with the family active, its row a - s for
// the sign s every active sign vector has on F, with rhs less their terms in the fixed weights. A weight off F is
// fixed exactly: at 0 on D, at its bound under an active bound row (one entry -1 or 1). normals has full column
// rank: a row that depends over F on the rows before it is left out. Rows independent over all n weights can
// depend over F (the budget 1 and the generated row -s with s = 1 on F differ only on D); every row is active at
// the point, so the rows kept imply the ones left out on the face. The multipliers stay those of every active
// row over all n weights.
struct FaceRows
{
    std::vector<std::size_t> free;
    // The fixed weights, 0 on F.
    std::vector<double> fixed;
    std::vector<double> normals;
    std::vector<double> rhs;
};

// The weight a row of one entry -1 or 1 bounds.
std::optional<std::size_t> bounded_weight(std::span<const double> row)
{
    auto nonzero = [](double entry) { return entry != 0.0; };
    if (std::ranges::count_if(row, nonzero) != 1)
    {
        return std::nullopt;
    }
    auto weight = static_cast<std::size_t>(std::ranges::find_if(row, nonzero) - row.begin());
    if (std::fabs(row[weight]) != 1.0)
    {
        return std::nullopt;
    }
    return weight;
}

// The row over F as the next column, its rhs less its term in the fixed weights.
void append_row(FaceRows& face, std::span<const double> row, double rhs, std::vector<double>& columns)
{
    for (std::size_t coordinate : face.free)
    {
        columns.push_back(row[coordinate]);
    }
    face.rhs.push_back(rhs - dot(row, face.fixed));
}

std::vector<double> transposed(std::span<const double> a, std::size_t rows, std::size_t cols)
{
    std::vector<double> result(rows * cols);
    for (std::size_t i = 0; i < rows; ++i)
    {
        for (std::size_t j = 0; j < cols; ++j)
        {
            result[j * rows + i] = a[i * cols + j];
        }
    }
    return result;
}

// The active rows of point that are not bounds; each bound fixes its weight in face.
std::vector<std::size_t> fix_bounds(const qp::Problem& problem, const PathPoint& point, FaceRows& face,
                                    std::vector<bool>& off_face)
{
    std::size_t n = problem.n;
    std::vector<std::size_t> rows;
    for (std::size_t row : point.solution.active_inequalities)
    {
        std::span<const double> normal = std::span<const double>(problem.inequality_rows).subspan(row * n, n);
        std::optional<std::size_t> weight = bounded_weight(normal);
        if (!weight.has_value())
        {
            rows.push_back(row);
            continue;
        }
        face.fixed[*weight] = problem.inequality_rhs[row] * normal[*weight];
        off_face[*weight] = true;
    }
    return rows;
}

// rows (count x free, row-major, a row per column) have full column rank as a free x count matrix.
bool full_column_rank(std::span<const double> rows, std::size_t count, std::size_t free)
{
    if (count > free)
    {
        return false;
    }
    return linalg::has_full_column_rank(linalg::householder_qr(transposed(rows, count, free), free, count));
}

// normals from rows (face.rhs.size() x |F|, row-major), each row kept when it is independent of the rows kept
// before it; the rhs of a row left out goes with it.
void keep_independent_rows(FaceRows& face, std::span<const double> rows)
{
    std::size_t p = face.free.size();
    std::vector<double> kept;
    std::vector<double> rhs;
    for (std::size_t k = 0; k < face.rhs.size(); ++k)
    {
        std::span<const double> row = rows.subspan(k * p, p);
        kept.insert(kept.end(), row.begin(), row.end());
        if (!full_column_rank(kept, rhs.size() + 1, p))
        {
            kept.resize(kept.size() - p);
            continue;
        }
        rhs.push_back(face.rhs[k]);
    }
    face.normals = transposed(kept, rhs.size(), p);
    face.rhs = std::move(rhs);
}

FaceRows face_rows(const qp::Problem& problem, const PathPoint& point, std::span<const double> signs)
{
    std::size_t n = problem.n;
    FaceRows face;
    face.fixed.assign(n, 0.0);
    std::vector<bool> off_face(n, false);
    std::vector<std::size_t> active = fix_bounds(problem, point, face, off_face);
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!off_face[i] && (signs.empty() || signs[i] != 0.0))
        {
            face.free.push_back(i);
        }
    }
    std::vector<double> columns;
    for (std::size_t k = 0; k < problem.equality_rhs.size(); ++k)
    {
        append_row(face, std::span<const double>(problem.equality_rows).subspan(k * n, n), problem.equality_rhs[k],
                   columns);
    }
    for (std::size_t row : active)
    {
        append_row(face, std::span<const double>(problem.inequality_rows).subspan(row * n, n),
                   problem.inequality_rhs[row], columns);
    }
    if (!signs.empty())
    {
        append_row(face, difference(problem.gross_exposure->offset, signs), problem.gross_exposure->rhs, columns);
    }
    keep_independent_rows(face, columns);
    return face;
}

// The face in the coordinates y = Q'w_F of the QR N_F = Q R of its rows: N_F'w_F = rhs fixes y_1 = R'^-1 rhs (the
// first m entries) and leaves y_2 free, so the last r = |F| - m columns of Q are a basis Z of the face. factor is
// L'E_F Q (n x |F|, row-major, E_F the columns F of the identity), and hessian the QR of its last r columns L'E_F Z,
// whose R_r'R_r is the reduced Hessian Z'Sigma_FF Z.
struct ReducedFace
{
    FaceRows rows;
    linalg::Qr normals;
    std::vector<double> factor;
    linalg::Qr hessian;
};

// L'E_F Q: row j is Q' applied to column j of E_F'L, the entries L_ij for i in F.
std::vector<double> factor_on_face(std::span<const double> factor, std::size_t n, const FaceRows& rows,
                                   const linalg::Qr& normals)
{
    std::size_t p = rows.free.size();
    std::vector<double> product(n * p);
    for (std::size_t j = 0; j < n; ++j)
    {
        std::span<double> column = std::span<double>(product).subspan(j * p, p);
        for (std::size_t a = 0; a < p; ++a)
        {
            column[a] = factor[rows.free[a] * n + j];
        }
        linalg::apply_transposed(normals, column);
    }
    return product;
}

// The columns first, first + 1, ... of a (rows x cols, row-major).
std::vector<double> columns_from(std::span<const double> a, std::size_t rows, std::size_t cols, std::size_t first)
{
    std::size_t count = cols - first;
    std::vector<double> result(rows * count);
    for (std::size_t i = 0; i < rows; ++i)
    {
        std::copy_n(a.begin() + static_cast<std::ptrdiff_t>(i * cols + first), count,
                    result.begin() + static_cast<std::ptrdiff_t>(i * count));
    }
    return result;
}

// The face reduced to its null space (its rows have full column rank, by face_rows): nullopt when the reduced
// Hessian has a pivot below its noise.
std::optional<ReducedFace> reduced_face(const Path& path, FaceRows rows)
{
    const qp::Problem& problem = path.translation.problem;
    std::optional<ReducedFace> face(std::in_place);
    face->rows = std::move(rows);
    std::size_t p = face->rows.free.size();
    std::size_t m = face->rows.rhs.size();
    face->normals = linalg::householder_qr(face->rows.normals, p, m);
    face->factor = factor_on_face(path.factor, problem.n, face->rows, face->normals);
    face->hessian = linalg::householder_qr(columns_from(face->factor, problem.n, p, m), problem.n, p - m);
    if (!linalg::has_full_column_rank(face->hessian))
    {
        face.reset();
    }
    return face;
}

// w(t) = weights + (t - point.t) direction on the face, over all n weights, w_D = 0.
struct FaceLine
{
    std::vector<double> weights;
    std::vector<double> direction;
};

// base with Q [y_1; y_2] on F.
std::vector<double> from_face(const ReducedFace& face, std::span<const double> fixed, std::span<const double> free,
                              std::vector<double> base)
{
    std::vector<double> y(fixed.begin(), fixed.end());
    y.insert(y.end(), free.begin(), free.end());
    linalg::apply(face.normals, y);
    for (std::size_t a = 0; a < y.size(); ++a)
    {
        base[face.rows.free[a]] = y[a];
    }
    return base;
}

// The stationarity of the shifted Q(t) on the face, Z'(Sigma w - t v) = 0, solved in y: with w_p the fixed weights
// and Q [y_1; 0] on F, and L'E_F Z = Q_r R_r, y_2(t) = R_r^-1 (t R_r'^-1 Z'v - (Q_r'L'w_p)_r), so
// dy_2 = R_r^-1 R_r'^-1 Z'v.
FaceLine face_line(const Path& path, const ReducedFace& face, const PathPoint& point)
{
    std::size_t n = point.solution.x.size();
    std::size_t p = face.rows.free.size();
    std::size_t m = face.rows.rhs.size();
    std::vector<double> fixed = face.rows.rhs;
    linalg::solve_upper_transposed(face.normals, fixed);
    std::vector<double> particular = factor_transposed_times(path.factor, face.rows.fixed);
    for (std::size_t j = 0; j < n; ++j)
    {
        particular[j] += dot(std::span<const double>(face.factor).subspan(j * p, m), fixed);
    }
    linalg::apply_transposed(face.hessian, particular);
    std::vector<double> drift(p);
    for (std::size_t a = 0; a < p; ++a)
    {
        drift[a] = point.shift.excess[face.rows.free[a]];
    }
    linalg::apply_transposed(face.normals, drift);
    drift.erase(drift.begin(), drift.begin() + static_cast<std::ptrdiff_t>(m));
    linalg::solve_upper_transposed(face.hessian, drift);
    std::vector<double> anchored(p - m);
    for (std::size_t b = 0; b < anchored.size(); ++b)
    {
        anchored[b] = point.t * drift[b] - particular[b];
    }
    linalg::solve_upper(face.hessian, anchored);
    linalg::solve_upper(face.hessian, drift);
    std::vector<double> none(m, 0.0);
    return FaceLine{from_face(face, fixed, anchored, face.rows.fixed),
                    from_face(face, none, drift, std::vector<double>(n, 0.0))};
}

// The face in the range space of Sigma = L L', when every weight is free. In y = L'w the face solve minimizes
// 0.5 |y|^2 - t (L^-1 v)'y over G'y = rhs for G = L^-1 N = Q R, so y = Q [R'^-1 rhs; (Q't L^-1 v)_2],
// dy = Q [0; (Q'L^-1 v)_2] and w = L'^-1 y. nullopt when L^-1 N has a pivot below its noise.
std::optional<linalg::Qr> range_normals(const Path& path, const FaceRows& rows)
{
    std::size_t n = rows.free.size();
    std::size_t m = rows.rhs.size();
    std::vector<double> solved(n * m);
    std::vector<double> column(n);
    for (std::size_t k = 0; k < m; ++k)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            column[i] = rows.normals[i * m + k];
        }
        linalg::solve_lower(path.factor, n, column);
        for (std::size_t i = 0; i < n; ++i)
        {
            solved[i * m + k] = column[i];
        }
    }
    std::optional<linalg::Qr> normals(std::in_place, linalg::householder_qr(solved, n, m));
    if (!linalg::has_full_column_rank(*normals))
    {
        normals.reset();
    }
    return normals;
}

// L'^-1 Q [head; (Q'L^-1 right)_2].
std::vector<double> from_range(const Path& path, const linalg::Qr& normals, std::span<const double> head,
                               std::vector<double> right)
{
    std::size_t n = right.size();
    linalg::solve_lower(path.factor, n, right);
    linalg::apply_transposed(normals, right);
    std::ranges::copy(head, right.begin());
    linalg::apply(normals, right);
    linalg::solve_lower_transposed(path.factor, n, right);
    return right;
}

FaceLine range_line(const Path& path, const FaceRows& rows, const linalg::Qr& normals, const PathPoint& point)
{
    std::vector<double> head = rows.rhs;
    linalg::solve_upper_transposed(normals, head);
    std::vector<double> anchored = point.shift.excess;
    for (double& entry : anchored)
    {
        entry *= point.t;
    }
    std::vector<double> none(head.size(), 0.0);
    return FaceLine{from_range(path, normals, head, std::move(anchored)),
                    from_range(path, normals, none, point.shift.excess)};
}

// The branch of fewer flops (leading terms, a multiply and an add counted apart) on a face whose n weights are all
// free, with m rows and r = n - m: the null space applies Q' to the n columns of L, 4 n^2 m, and factors L'Z
// (n x r), 2 n r^2; the range space makes m + 4 triangular solves with L, (m + 4) n^2. The QR of the n x m rows
// (N or L^-1 N) is common to both. A face with a fixed weight stays in the null space: its range space would need
// the factor of Sigma_FF, n p^2 / 3 to n p^2 flops to form from the rows F of L and p^3 / 3 to factor, the order of
// the null space's 2 n r^2 itself, and a face the certificate refuses would pay for both branches.
bool range_space_cheaper(std::size_t n, std::size_t m)
{
    double assets = static_cast<double>(n);
    double rows = static_cast<double>(m);
    double dimension = assets - rows;
    double null_space = 4.0 * assets * assets * rows + 2.0 * assets * dimension * dimension;
    double range_space = (rows + 4.0) * assets * assets;
    return range_space < null_space;
}

// The h with |L'(weights + h direction)| = s: |L'w|^2 = c + 2 b h + a h^2 reaches s^2 at
// h = (s^2 - c) / (b + sqrt(b^2 + a (s^2 - c))), without the cancellation of -b + sqrt(...) since b >= 0 where the
// volatility is nondecreasing. nullopt when there is no such root next to h = 0.
std::optional<double> crossing_step(const Path& path, const FaceLine& line)
{
    std::vector<double> start = factor_transposed_times(path.factor, line.weights);
    std::vector<double> direction = factor_transposed_times(path.factor, line.direction);
    double volatility = std::sqrt(dot(start, start));
    double curvature = dot(direction, direction);
    double slope = dot(start, direction);
    double gap = (path.cap - volatility) * (path.cap + volatility);
    double denominator = slope + std::sqrt(slope * slope + curvature * gap);
    if (!(denominator > 0.0))
    {
        return std::nullopt;
    }
    return gap / denominator;
}

// A point on the face line of a point at t, with its multipliers u_v of N u_v = Sigma w - t v in the frame of the
// point's shift.
struct Candidate
{
    double t = 0.0;
    std::vector<double> weights;
    Slope multipliers;
};

Candidate candidate_at(const Path& path, const PathPoint& point, double t, std::vector<double> weights)
{
    std::vector<double> gradient = factor_times(path.factor, factor_transposed_times(path.factor, weights));
    for (std::size_t i = 0; i < gradient.size(); ++i)
    {
        gradient[i] -= t * point.shift.excess[i];
    }
    Slope multipliers = active_multipliers(path.translation.problem, point, std::move(gradient));
    return Candidate{t, std::move(weights), std::move(multipliers)};
}

// The rounding of the multipliers of the shifted Q(t*): n eps of the terms of Sigma w - t* v.
double dual_tolerance(const Path& path, const PathPoint& point, const Candidate& candidate)
{
    double n = static_cast<double>(candidate.weights.size());
    double scale = candidate.t * largest_magnitude(point.shift.excess) +
                   largest_variance(path.factor, candidate.weights.size()) * largest_magnitude(candidate.weights);
    return kNoiseFactor * n * kEpsilon * scale;
}

// Every row and the family hold to the feasibility tolerance of the QP solver, and each weight of F keeps the sign
// of the face (a weight past zero is on another face).
bool primal_holds(const Path& path, std::span<const double> signs, std::span<const double> weights)
{
    const qp::Problem& problem = path.translation.problem;
    if (!std::ranges::all_of(weights, [](double weight) { return std::isfinite(weight); }) ||
        !rows_hold(problem, weights))
    {
        return false;
    }
    if (signs.empty())
    {
        return true;
    }
    double row_norm = std::sqrt(static_cast<double>(problem.n));
    double tolerance = qp::feasibility_tolerance(problem.gross_exposure->rhs, row_norm, largest_magnitude(weights));
    for (std::size_t i = 0; i < weights.size(); ++i)
    {
        if (signs[i] * weights[i] < -tolerance)
        {
            return false;
        }
    }
    return true;
}

// The active inequality rows have nonnegative multipliers in Q(t*): u_v less t* amount at the face of the shift.
bool explicit_duals_hold(const PathPoint& point, const Candidate& candidate, double tolerance)
{
    Slope multipliers = unshifted(candidate.multipliers, point, candidate.t * point.shift.amount);
    return std::ranges::all_of(point.solution.active_inequalities,
                               [&](std::size_t row) { return multipliers.inequality[row] >= -tolerance; });
}

// rho = t* v - Sigma w + E'y_v + C'z_v: what stationarity leaves to the family, sum_j zeta_v,j (s_j - a).
std::vector<double> family_residual(const Path& path, const PathPoint& point, const Candidate& candidate)
{
    const qp::Problem& problem = path.translation.problem;
    std::size_t n = problem.n;
    std::vector<double> residual = factor_times(path.factor, factor_transposed_times(path.factor, candidate.weights));
    std::vector<double> equality = transposed_times(problem.equality_rows, candidate.multipliers.equality, n);
    std::vector<double> inequality = transposed_times(problem.inequality_rows, candidate.multipliers.inequality, n);
    for (std::size_t i = 0; i < n; ++i)
    {
        residual[i] = candidate.t * point.shift.excess[i] - residual[i] + equality[i] + inequality[i];
    }
    return residual;
}

// The family's term in the stationarity of Q(t*) is zeta (o - g) for its offset o, with zeta = sum_j zeta_j >= 0 and g
// in the subdifferential of ||w||_1: g_i = sigma_i on the support, |g_i| <= 1 on D. In the frame of the shift, with
// zeta_v the sum of the shifted multipliers and q = rho + zeta_v o, zeta = zeta_v - t* amount and
// zeta g = q - t* amount s_f when the shift is along a generated face s_f (amount 0 for any other face), so:
// sigma_i q_i = zeta_v on the support and -zeta_v + 2 t* amount <= s_f,i q_i <= zeta_v on D. These are quantities
// of the order of t*|v| + |Sigma||w|, except 2 t* amount, and hold for every set of sign vectors of the face.
struct FamilyLoad
{
    double total = 0.0;
    double shift = 0.0;
    std::span<const double> face_signs;
};

bool subgradient_holds(const FamilyLoad& load, double sign, double face_sign, double q, double tolerance)
{
    if (sign != 0.0)
    {
        return std::fabs(sign * q - load.total) <= tolerance;
    }
    double along = face_sign * q;
    return along <= load.total + tolerance && along >= -load.total + 2.0 * load.shift - tolerance;
}

FamilyLoad family_load(const PathPoint& point, const Candidate& candidate)
{
    FamilyLoad load;
    for (double zeta : candidate.multipliers.gross_exposure)
    {
        load.total += zeta;
    }
    if (point.shift.face.kind == Face::Kind::Generated)
    {
        load.shift = candidate.t * point.shift.amount;
        load.face_signs = point.shift.face.signs;
    }
    return load;
}

bool family_dual_holds(const Path& path, const PathPoint& point, std::span<const double> signs,
                       const Candidate& candidate, double tolerance)
{
    std::vector<double> residual = family_residual(path, point, candidate);
    if (signs.empty())
    {
        return std::ranges::all_of(residual, [&](double entry) { return std::fabs(entry) <= tolerance; });
    }
    FamilyLoad load = family_load(point, candidate);
    if (load.total - load.shift < -tolerance)
    {
        return false;
    }
    std::span<const double> offset(path.translation.problem.gross_exposure->offset);
    for (std::size_t i = 0; i < residual.size(); ++i)
    {
        double face_sign = 1.0;
        if (!load.face_signs.empty())
        {
            face_sign = load.face_signs[i];
        }
        if (!subgradient_holds(load, signs[i], face_sign, residual[i] + load.total * offset[i], tolerance))
        {
            return false;
        }
    }
    return true;
}

// The KKT conditions of Q(t*) at the candidate, which suffice since Q(t*) is convex: stationarity and
// complementarity hold on the face by construction, so primal feasibility and the signs of the multipliers remain.
bool crossing_certified(const Path& path, const PathPoint& point, std::span<const double> signs,
                        const Candidate& candidate)
{
    double tolerance = dual_tolerance(path, point, candidate);
    return primal_holds(path, signs, candidate.weights) && explicit_duals_hold(point, candidate, tolerance) &&
           family_dual_holds(path, point, signs, candidate, tolerance);
}

// The candidate as a solution of Q(t*) in the convention of PathPoint::solution: the active set of point, its
// multipliers less t* amount at the face of the shift; an inequality multiplier within the tolerance below 0 is 0.
qp::Solution crossing_solution(const PathPoint& point, Candidate candidate)
{
    Slope multipliers = unshifted(std::move(candidate.multipliers), point, candidate.t * point.shift.amount);
    qp::Solution solution = point.solution;
    solution.x = std::move(candidate.weights);
    solution.equality_multipliers = std::move(multipliers.equality);
    solution.inequality_multipliers = std::move(multipliers.inequality);
    for (double& multiplier : solution.inequality_multipliers)
    {
        multiplier = nonnegative(multiplier);
    }
    solution.gross_exposure_multipliers = std::move(multipliers.gross_exposure);
    return solution;
}

// The crossing on the face line of point, when the KKT conditions of Q(t*) certify it.
std::optional<Crossing> certified_crossing(const Path& path, const PathPoint& point,
                                           std::span<const double> signs, FaceLine line)
{
    std::optional<double> step = crossing_step(path, line);
    if (!step.has_value() || !std::isfinite(point.t + *step) || !(point.t + *step > 0.0))
    {
        return std::nullopt;
    }
    for (std::size_t i = 0; i < line.weights.size(); ++i)
    {
        line.weights[i] += *step * line.direction[i];
    }
    Candidate candidate = candidate_at(path, point, point.t + *step, std::move(line.weights));
    if (!crossing_certified(path, point, signs, candidate))
    {
        return std::nullopt;
    }
    return Crossing{candidate.t, crossing_solution(point, std::move(candidate))};
}

std::optional<Crossing> null_space_crossing(const Path& path, const PathPoint& point,
                                            std::span<const double> signs, FaceRows rows)
{
    std::optional<ReducedFace> face = reduced_face(path, std::move(rows));
    if (!face.has_value())
    {
        return std::nullopt;
    }
    return certified_crossing(path, point, signs, face_line(path, *face, point));
}

// The weights of the range space carry eps cond(Sigma) where those of the null space carry eps cond(Z'Sigma Z), and
// its stationarity residual is eps |L| |g|, with g = t L^-1 v the unconstrained optimum in y, where the null space's
// is eps |Sigma| |w|. Neither error is taken on trust: the root puts vol = s at the computed weights by the
// cancellation-free formula, the multipliers come from the same least squares over N, and the certificate measures
// that residual against the dual tolerance. It passes where |g| is of the order of |y| = |L'w|, and refuses a face
// that holds the portfolio far from t Sigma^-1 v (the start of the Goldfarb-Idnani solve in
// docs/research/2026-11-18-max-return-path-and-portable-generators.md, section 1); there Cheaper takes the null space.
// nullopt as well when a weight is fixed.
std::optional<Crossing> range_space_crossing(const Path& path, const PathPoint& point,
                                             std::span<const double> signs, const FaceRows& rows)
{
    if (rows.free.size() < path.translation.problem.n)
    {
        return std::nullopt;
    }
    std::optional<linalg::Qr> normals = range_normals(path, rows);
    if (!normals.has_value())
    {
        return std::nullopt;
    }
    return certified_crossing(path, point, signs, range_line(path, rows, *normals, point));
}

std::optional<Crossing> range_or_null_space_crossing(const Path& path, const PathPoint& point,
                                                     std::span<const double> signs, FaceRows rows)
{
    std::optional<Crossing> crossing = range_space_crossing(path, point, signs, rows);
    if (!crossing.has_value())
    {
        crossing = null_space_crossing(path, point, signs, std::move(rows));
    }
    return crossing;
}

} // namespace

Path make_path(const Spec& spec, const Data& data)
{
    Translation translation = mean_variance_problem(spec, data, MeanVariance{1.0, 0.0});
    std::vector<double> excess = excess_of(data);
    Shift base = base_shift(translation, excess);
    return Path{std::move(translation), std::move(excess), data.covariance_factor(), *spec.max_volatility,
                std::move(base), Shift{}, 0, qp::BasisCache{}};
}

// Without a budget, the common mean of mu~ falls on the inequality face that holds the portfolio once it is
// active (a cash bound, the leverage): the face of the path first, then the most loaded face of Q(t) itself.
// Every shifted solution is certified, and Q(t) unshifted is the answer otherwise.
PathPoint solve_at(Path& path, double t)
{
    std::optional<PathPoint> point = on_path_face(path, t);
    if (point.has_value())
    {
        return std::move(*point);
    }
    point = solve_with(path, path.base, t);
    std::optional<PathPoint> shifted = on_loaded_face(path, *point);
    if (shifted.has_value())
    {
        return std::move(*shifted);
    }
    return std::move(*point);
}

bool same_segment(const PathPoint& earlier, const PathPoint& later)
{
    return earlier.solution.active_inequalities == later.solution.active_inequalities &&
           support_signs(earlier) == support_signs(later);
}

bool weights_moved(const PathPoint& earlier, const PathPoint& later)
{
    double change = largest_magnitude(difference(later.solution.x, earlier.solution.x));
    return change > weight_noise(earlier, later);
}

bool weights_constant(const Path& path, const PathPoint& earlier, const PathPoint& later)
{
    if (!same_segment(earlier, later))
    {
        return false;
    }
    return is_vertex(path, later) || !weights_moved(earlier, later);
}

Slope multiplier_slope(const Path& path, const PathPoint& point)
{
    std::vector<double> target(point.shift.excess.size());
    for (std::size_t i = 0; i < target.size(); ++i)
    {
        target[i] = -point.shift.excess[i];
    }
    return active_multipliers(path.translation.problem, point, std::move(target));
}

double vertex_exit(const PathPoint& point, const Slope& shifted)
{
    Slope slope = unshifted(shifted, point, point.shift.amount);
    const qp::Solution& solution = point.solution;
    double exit = kInfinity;
    for (std::size_t row : solution.active_inequalities)
    {
        exit = std::min(exit, zero_crossing(point, solution.inequality_multipliers[row], slope.inequality[row]));
    }
    for (std::size_t k = 0; k < generated_count(point); ++k)
    {
        exit = std::min(exit, zero_crossing(point, solution.gross_exposure_multipliers[k], slope.gross_exposure[k]));
    }
    return exit;
}

bool limit_holds(const Path& path, const PathPoint& point, const Slope& slope)
{
    return rows_hold(path.translation.problem, point.solution.x) && slope_stationary(path, point, slope);
}

qp::Solution moved_along(const PathPoint& point, const Slope& shifted, double t)
{
    Slope slope = unshifted(shifted, point, point.shift.amount);
    qp::Solution solution = point.solution;
    double distance = t - point.t;
    for (std::size_t k = 0; k < solution.equality_multipliers.size(); ++k)
    {
        solution.equality_multipliers[k] += distance * slope.equality[k];
    }
    for (std::size_t k = 0; k < solution.inequality_multipliers.size(); ++k)
    {
        solution.inequality_multipliers[k] = nonnegative(solution.inequality_multipliers[k] +
                                                         distance * slope.inequality[k]);
    }
    for (std::size_t k = 0; k < solution.gross_exposure_multipliers.size(); ++k)
    {
        solution.gross_exposure_multipliers[k] = nonnegative(solution.gross_exposure_multipliers[k] +
                                                             distance * slope.gross_exposure[k]);
    }
    return solution;
}

qp::Solution slope_solution(const PathPoint& point, const Slope& shifted)
{
    Slope slope = unshifted(shifted, point, point.shift.amount);
    qp::Solution solution = point.solution;
    solution.equality_multipliers = slope.equality;
    for (std::size_t k = 0; k < solution.inequality_multipliers.size(); ++k)
    {
        solution.inequality_multipliers[k] = nonnegative(slope.inequality[k]);
    }
    for (std::size_t k = 0; k < solution.gross_exposure_multipliers.size(); ++k)
    {
        solution.gross_exposure_multipliers[k] = nonnegative(slope.gross_exposure[k]);
    }
    return solution;
}

// A face with as many rows as free weights has no direction.
std::optional<Crossing> tangent_crossing(const Path& path, const PathPoint& point, FaceSolve solve)
{
    if (point.solution.status != qp::Status::Optimal)
    {
        return std::nullopt;
    }
    std::vector<double> signs = common_signs(point);
    FaceRows rows = face_rows(path.translation.problem, point, signs);
    std::size_t n = path.translation.problem.n;
    std::size_t p = rows.free.size();
    std::size_t m = rows.rhs.size();
    if (m >= p)
    {
        return std::nullopt;
    }
    if (solve == FaceSolve::RangeSpace)
    {
        return range_space_crossing(path, point, signs, rows);
    }
    if (solve == FaceSolve::NullSpace || p < n || !range_space_cheaper(n, m))
    {
        return null_space_crossing(path, point, signs, std::move(rows));
    }
    return range_or_null_space_crossing(path, point, signs, std::move(rows));
}
} // namespace po::problems
