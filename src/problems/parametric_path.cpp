#include "problems/parametric_path.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
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

// The rounding of Q(t), as a multiple of the noise models below. On equal-mean paths, where w(t) is constant,
// the measured weight change stays below 0.3 times weight_noise / kNoiseFactor up to a condition number of
// 1e8 (docs/guide/validation.md).
constexpr double kNoiseFactor = 16.0;
constexpr double kInfinity = std::numeric_limits<double>::infinity();
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

double free_scale_of(const Data& data, std::vector<double> excess)
{
    std::size_t n = data.asset_count();
    linalg::solve_lower(data.covariance_factor(), n, excess);
    linalg::solve_lower_transposed(data.covariance_factor(), n, excess);
    return largest_magnitude(excess);
}

// The solver starts from the unconstrained minimizer t Sigma^-1 mu~ and every later step carries its
// rounding: a weight is known to n eps (max(1, |w|max) + t |Sigma^-1 mu~|max), the noise of the stationarity
// equation.
double weight_noise(const Path& path, const PathPoint& earlier, const PathPoint& later)
{
    double n = static_cast<double>(later.solution.x.size());
    double scale =
        std::max({1.0, largest_magnitude(earlier.solution.x), largest_magnitude(later.solution.x)});
    return kNoiseFactor * n * kEpsilon * (scale + later.t * path.free_scale);
}

// b from N b = -mu~ by Householder least squares is known to n eps |mu~|max.
double slope_noise(const Path& path)
{
    return kNoiseFactor * static_cast<double>(path.excess.size()) * kEpsilon * largest_magnitude(path.excess);
}

// The solver keeps the active normals independent, so the equality rows and the active rows have rank equal
// to their count, and a count of n fixes w. (MAX_RETURN has one equality row at most, the budget.)
bool is_vertex(const Path& path, const PathPoint& point)
{
    const qp::Problem& problem = path.translation.problem;
    return problem.equality_rhs.size() + point.solution.active_inequalities.size() == problem.n;
}

// Normal k of the active set: the equality rows first, then the active rows in ascending order.
std::span<const double> active_normal(const qp::Problem& problem, const PathPoint& point, std::size_t k)
{
    std::size_t n = problem.n;
    std::size_t equalities = problem.equality_rhs.size();
    if (k < equalities)
    {
        return std::span<const double>(problem.equality_rows).subspan(k * n, n);
    }
    std::size_t row = point.solution.active_inequalities[k - equalities];
    return std::span<const double>(problem.inequality_rows).subspan(row * n, n);
}

// N (n x m, row-major), the active normals as columns.
std::vector<double> active_normals(const qp::Problem& problem, const PathPoint& point)
{
    std::size_t n = problem.n;
    std::size_t m = problem.equality_rhs.size() + point.solution.active_inequalities.size();
    std::vector<double> normals(n * m);
    for (std::size_t k = 0; k < m; ++k)
    {
        std::span<const double> normal = active_normal(problem, point, k);
        for (std::size_t i = 0; i < n; ++i)
        {
            normals[i * m + k] = normal[i];
        }
    }
    return normals;
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
    std::vector<double> excess = excess_of(data);
    double free_scale = free_scale_of(data, excess);
    return Path{mean_variance_problem(spec, data, MeanVariance{1.0, 0.0}), std::move(excess),
                data.covariance_factor(), *spec.max_volatility, free_scale, 0};
}

PathPoint solve_at(Path& path, double t)
{
    qp::Problem& problem = path.translation.problem;
    for (std::size_t i = 0; i < problem.n; ++i)
    {
        problem.linear[i] = -t * path.excess[i];
    }
    PathPoint point{t, qp::solve(problem), std::numeric_limits<double>::quiet_NaN()};
    path.iterations += point.solution.iterations;
    if (point.solution.status == qp::Status::Optimal)
    {
        std::vector<double> root = factor_transposed_times(path.factor, point.solution.x);
        point.volatility = std::sqrt(dot(root, root));
    }
    return point;
}

bool same_segment(const PathPoint& earlier, const PathPoint& later)
{
    return earlier.solution.active_inequalities == later.solution.active_inequalities;
}

bool weights_moved(const Path& path, const PathPoint& earlier, const PathPoint& later)
{
    double change = largest_magnitude(difference(later.solution.x, earlier.solution.x));
    return change > weight_noise(path, earlier, later);
}

bool weights_constant(const Path& path, const PathPoint& earlier, const PathPoint& later)
{
    if (!same_segment(earlier, later))
    {
        return false;
    }
    return is_vertex(path, later) || !weights_moved(path, earlier, later);
}

Slope multiplier_slope(const Path& path, const PathPoint& point)
{
    const qp::Problem& problem = path.translation.problem;
    std::size_t equalities = problem.equality_rhs.size();
    std::size_t m = equalities + point.solution.active_inequalities.size();
    std::vector<double> target(path.excess.size());
    for (std::size_t i = 0; i < target.size(); ++i)
    {
        target[i] = -path.excess[i];
    }
    std::vector<double> b =
        linalg::least_squares(linalg::householder_qr(active_normals(problem, point), problem.n, m), target);
    Slope slope{std::vector<double>(b.begin(), b.begin() + static_cast<std::ptrdiff_t>(equalities)),
                std::vector<double>(problem.inequality_rhs.size(), 0.0)};
    for (std::size_t k = equalities; k < m; ++k)
    {
        slope.inequality[point.solution.active_inequalities[k - equalities]] = b[k];
    }
    return slope;
}

double vertex_exit(const Path& path, const PathPoint& point, const Slope& slope)
{
    double noise = slope_noise(path);
    double exit = kInfinity;
    for (std::size_t row : point.solution.active_inequalities)
    {
        double rate = slope.inequality[row];
        if (!(rate < -noise))
        {
            continue;
        }
        exit = std::min(exit, point.t - point.solution.inequality_multipliers[row] / rate);
    }
    return exit;
}

qp::Solution moved_along(const PathPoint& point, const Slope& slope, double t)
{
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
    return solution;
}

qp::Solution slope_solution(const PathPoint& point, const Slope& slope)
{
    qp::Solution solution = point.solution;
    solution.equality_multipliers = slope.equality;
    for (std::size_t k = 0; k < solution.inequality_multipliers.size(); ++k)
    {
        solution.inequality_multipliers[k] = nonnegative(slope.inequality[k]);
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
    for (std::size_t i = 0; i < solution.x.size(); ++i)
    {
        solution.x[i] = (1.0 - fraction) * lower.x[i] + fraction * upper.x[i];
    }
    for (std::size_t k = 0; k < solution.equality_multipliers.size(); ++k)
    {
        solution.equality_multipliers[k] =
            (1.0 - fraction) * lower.equality_multipliers[k] + fraction * upper.equality_multipliers[k];
    }
    for (std::size_t k = 0; k < solution.inequality_multipliers.size(); ++k)
    {
        solution.inequality_multipliers[k] =
            (1.0 - fraction) * lower.inequality_multipliers[k] + fraction * upper.inequality_multipliers[k];
    }
    return solution;
}

} // namespace po::problems
