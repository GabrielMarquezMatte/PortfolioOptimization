#include "problems/parametric_path.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <utility>
#include <vector>

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
    PathPoint point{t, qp::solve(problem), kNaN, shift};
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

// The slope b of the multipliers of Q(t) from that of the shifted Q(t): b_v less the shift amount at the face, the
// face of point being among its active rows.
Slope unshifted(Slope slope, const PathPoint& point)
{
    double amount = -point.shift.amount;
    const Face& face = point.shift.face;
    if (face.kind == Face::Kind::Equality)
    {
        slope.equality[face.row] += amount;
        return slope;
    }
    if (face.kind == Face::Kind::Inequality)
    {
        slope.inequality[face.row] += amount;
        return slope;
    }
    if (face.kind == Face::Kind::Generated)
    {
        slope.gross_exposure[sign_position(point.solution, face.signs).first] += amount;
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

// (1 - fraction) lower + fraction upper, entry by entry.
void interpolate(std::span<const double> lower, std::span<const double> upper, double fraction,
                 std::span<double> result)
{
    for (std::size_t k = 0; k < result.size(); ++k)
    {
        result[k] = (1.0 - fraction) * lower[k] + fraction * upper[k];
    }
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

} // namespace

Path make_path(const Spec& spec, const Data& data)
{
    Translation translation = mean_variance_problem(spec, data, MeanVariance{1.0, 0.0});
    std::vector<double> excess = excess_of(data);
    Shift base = base_shift(translation, excess);
    return Path{std::move(translation), std::move(excess), data.covariance_factor(), *spec.max_volatility,
                std::move(base), Shift{}, 0};
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
           earlier.solution.gross_exposure_signs == later.solution.gross_exposure_signs;
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
    const qp::Problem& problem = path.translation.problem;
    std::size_t equalities = problem.equality_rhs.size();
    std::size_t inequalities = point.solution.active_inequalities.size();
    std::vector<double> target(point.shift.excess.size());
    for (std::size_t i = 0; i < target.size(); ++i)
    {
        target[i] = -point.shift.excess[i];
    }
    std::vector<double> b = linalg::least_squares(
        linalg::householder_qr(active_normals(problem, point), problem.n, active_count(problem, point)), target);
    auto generated = b.begin() + static_cast<std::ptrdiff_t>(equalities + inequalities);
    Slope slope{std::vector<double>(b.begin(), b.begin() + static_cast<std::ptrdiff_t>(equalities)),
                std::vector<double>(problem.inequality_rhs.size(), 0.0), std::vector<double>(generated, b.end())};
    for (std::size_t k = 0; k < inequalities; ++k)
    {
        slope.inequality[point.solution.active_inequalities[k]] = b[equalities + k];
    }
    return slope;
}

double vertex_exit(const PathPoint& point, const Slope& shifted)
{
    Slope slope = unshifted(shifted, point);
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
    Slope slope = unshifted(shifted, point);
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
    Slope slope = unshifted(shifted, point);
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

// |L'w|^2 = c + 2 b tau + a tau^2 reaches s^2 at the root in [0, 1], written without the cancellation of
// -b + sqrt(...) since b >= 0 where the volatility is nondecreasing.
double crossing_fraction(const Path& path, const Bracket& bracket)
{
    std::vector<double> start = factor_transposed_times(path.factor, bracket.lower.solution.x);
    std::vector<double> direction =
        factor_transposed_times(path.factor, difference(bracket.upper.solution.x, bracket.lower.solution.x));
    double curvature = dot(direction, direction);
    double slope = dot(start, direction);
    double gap = (path.cap - bracket.lower.volatility) * (path.cap + bracket.lower.volatility);
    return std::clamp(gap / (slope + std::sqrt(slope * slope + curvature * gap)), 0.0, 1.0);
}

qp::Solution on_segment(const Bracket& bracket, double fraction)
{
    const qp::Solution& lower = bracket.lower.solution;
    const qp::Solution& upper = bracket.upper.solution;
    qp::Solution solution = lower;
    interpolate(lower.x, upper.x, fraction, solution.x);
    interpolate(lower.equality_multipliers, upper.equality_multipliers, fraction, solution.equality_multipliers);
    interpolate(lower.inequality_multipliers, upper.inequality_multipliers, fraction,
                solution.inequality_multipliers);
    interpolate(lower.gross_exposure_multipliers, upper.gross_exposure_multipliers, fraction,
                solution.gross_exposure_multipliers);
    return solution;
}

} // namespace po::problems
