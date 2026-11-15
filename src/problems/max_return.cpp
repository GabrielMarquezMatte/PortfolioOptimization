#include "problems/max_return.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include "problems/constraints.hpp"
#include "problems/mean_variance.hpp"
#include "qp/goldfarb_idnani.hpp"

namespace po::problems
{
namespace
{

constexpr double kVolatilityTolerance = 1e-12;
constexpr double kSameWeightsTolerance = 1e-14;
constexpr int kMaxDoublings = 200;
constexpr int kMaxBisections = 200;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInfinity = std::numeric_limits<double>::infinity();
constexpr double kEpsilon = std::numeric_limits<double>::epsilon();

// Q(t): G = Sigma with the factor L of Data, and g = -t (mu - rf 1); only g depends on t. cap is s.
struct Path
{
    Translation translation;
    std::vector<double> excess;
    std::span<const double> factor;
    double cap = 0.0;
    int iterations = 0;
};

struct PathPoint
{
    double t = 0.0;
    qp::Solution solution;
    double volatility = kNaN;
};

// vol(w(lower.t)) < s <= vol(w(upper.t)).
struct Bracket
{
    PathPoint lower;
    PathPoint upper;
};

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

Path make_path(const Spec& spec, const Data& data)
{
    return Path{mean_variance_problem(spec, data, MeanVariance{1.0, 0.0}), excess_of(data),
                data.covariance_factor(), *spec.max_volatility, 0};
}

PathPoint solve_at(Path& path, double t)
{
    qp::Problem& problem = path.translation.problem;
    for (std::size_t i = 0; i < problem.n; ++i)
    {
        problem.linear[i] = -t * path.excess[i];
    }
    PathPoint point{t, qp::solve(problem), kNaN};
    path.iterations += point.solution.iterations;
    if (point.solution.status == qp::Status::Optimal)
    {
        std::vector<double> root = factor_transposed_times(path.factor, point.solution.x);
        point.volatility = std::sqrt(dot(root, root));
    }
    return point;
}

Result stopped(const Data& data, const Path& path, qp::Status status)
{
    qp::Solution solution;
    solution.status = status;
    solution.iterations = path.iterations;
    return to_result(data, path.translation.rows, solution);
}

Result answer(const Data& data, const Path& path, const qp::Solution& solution, double multiplier_scale,
              double volatility_multiplier)
{
    Result result = to_result(data, path.translation.rows, solution, multiplier_scale);
    result.volatility_multiplier = volatility_multiplier;
    result.iterations = path.iterations;
    return result;
}

// The weights of earlier with the multipliers a u(earlier) + b u(later). An inequality multiplier is
// nonnegative, and decreases() keeps every combination taken here at or above zero up to rounding, which is
// dropped.
qp::Solution combined(const PathPoint& earlier, const PathPoint& later, double a, double b)
{
    qp::Solution solution = earlier.solution;
    for (std::size_t k = 0; k < solution.equality_multipliers.size(); ++k)
    {
        solution.equality_multipliers[k] =
            a * earlier.solution.equality_multipliers[k] + b * later.solution.equality_multipliers[k];
    }
    for (std::size_t k = 0; k < solution.inequality_multipliers.size(); ++k)
    {
        double value =
            a * earlier.solution.inequality_multipliers[k] + b * later.solution.inequality_multipliers[k];
        solution.inequality_multipliers[k] = std::max(value, 0.0);
    }
    return solution;
}

// Where w(t) is constant, the QP multipliers are u(t) = a + b t and Sigma w - t (mu - rf 1) = N u(t) for
// every t, so -(mu - rf 1) = N b: the user multipliers are the slope b, and the volatility multiplier is 0.
Result at_limit(const Data& data, const Path& path, const PathPoint& earlier, const PathPoint& later)
{
    return answer(data, path, combined(earlier, later, -1.0, 1.0), later.t - earlier.t, 0.0);
}

// w(t) stays at w(earlier) up to exit, where an active multiplier reaches zero. There
// Sigma w - exit (mu - rf 1) = N u(exit), so the user multipliers are u(exit) / exit and the volatility
// multiplier is -s / exit, the right derivative of the optimum in s.
Result at_exit(const Data& data, const Path& path, const PathPoint& earlier, const PathPoint& later,
               double exit)
{
    double fraction = (exit - earlier.t) / (later.t - earlier.t);
    return answer(data, path, combined(earlier, later, 1.0 - fraction, fraction), exit, -path.cap / exit);
}

// The volatility grows as soon as t > 0: at s = sigma_min the frontier is vertical, the volatility multiplier
// is -infinity and the other multipliers have no finite value.
Result vertical(const Data& data, const Path& path, const PathPoint& start)
{
    Result result = to_result(data, path.translation.rows, start.solution);
    result.budget_multiplier = kNaN;
    result.target_return_multiplier = kNaN;
    std::ranges::fill(result.bound_multipliers, kNaN);
    result.volatility_multiplier = -kInfinity;
    result.iterations = path.iterations;
    return result;
}

double multiplier_magnitude(const PathPoint& point)
{
    return std::max(largest_magnitude(point.solution.equality_multipliers),
                    largest_magnitude(point.solution.inequality_multipliers));
}

// A multiplier whose slope is zero moves by rounding, which grows with t: only a larger fall is a decrease.
bool decreases(const PathPoint& earlier, const PathPoint& later, std::size_t row)
{
    double scale = std::max(multiplier_magnitude(earlier), multiplier_magnitude(later));
    double fall = earlier.solution.inequality_multipliers[row] - later.solution.inequality_multipliers[row];
    return fall > kSameWeightsTolerance * scale;
}

// With u(t) affine through both points, the t where the first decreasing active multiplier reaches zero;
// infinity when none decreases, so that w stays optimal for every larger t.
double vertex_exit(const PathPoint& earlier, const PathPoint& later)
{
    double exit = kInfinity;
    for (std::size_t row : later.solution.active_inequalities)
    {
        if (!decreases(earlier, later, row))
        {
            continue;
        }
        double start = earlier.solution.inequality_multipliers[row];
        double slope = (later.solution.inequality_multipliers[row] - start) / (later.t - earlier.t);
        exit = std::min(exit, earlier.t - start / slope);
    }
    return exit;
}

bool same_weights(std::span<const double> earlier, std::span<const double> later)
{
    double scale = std::max(largest_magnitude(earlier), largest_magnitude(later));
    return largest_magnitude(difference(later, earlier)) <= kSameWeightsTolerance * scale;
}

bool stays(const PathPoint& earlier, const PathPoint& later)
{
    return earlier.solution.active_inequalities == later.solution.active_inequalities &&
           same_weights(earlier.solution.x, later.solution.x);
}

// With w = w_a + tau (w_b - w_a), |L'w|^2 = c + 2 b tau + a tau^2 reaches s^2 at the root in [0, 1], written
// without the cancellation of -b + sqrt(...) since b >= 0 where the volatility is nondecreasing.
double crossing_t(const Path& path, const Bracket& bracket)
{
    std::vector<double> start = factor_transposed_times(path.factor, bracket.lower.solution.x);
    std::vector<double> direction =
        factor_transposed_times(path.factor, difference(bracket.upper.solution.x, bracket.lower.solution.x));
    double curvature = dot(direction, direction);
    double slope = dot(start, direction);
    double gap = (path.cap - bracket.lower.volatility) * (path.cap + bracket.lower.volatility);
    double fraction = gap / (slope + std::sqrt(slope * slope + curvature * gap));
    return bracket.lower.t + fraction * (bracket.upper.t - bracket.lower.t);
}

bool reaches_cap(const Path& path, const PathPoint& point)
{
    return std::fabs(point.volatility - path.cap) <= kVolatilityTolerance * path.cap;
}

void narrow(const Path& path, Bracket& bracket, PathPoint point)
{
    if (point.volatility >= path.cap)
    {
        bracket.upper = std::move(point);
        return;
    }
    bracket.lower = std::move(point);
}

// On one segment w(t) is affine and the root of vol(w(t)) = s is exact; a root that misses s by more than the
// tolerance narrows the bracket like any other point.
Result bisect(const Data& data, Path& path, Bracket bracket)
{
    for (int step = 0; step < kMaxBisections; ++step)
    {
        if (bracket.lower.solution.active_inequalities == bracket.upper.solution.active_inequalities)
        {
            PathPoint crossing = solve_at(path, crossing_t(path, bracket));
            if (crossing.solution.status != qp::Status::Optimal)
            {
                return stopped(data, path, crossing.solution.status);
            }
            if (reaches_cap(path, crossing))
            {
                return answer(data, path, crossing.solution, crossing.t, -path.cap / crossing.t);
            }
            narrow(path, bracket, std::move(crossing));
        }
        PathPoint middle = solve_at(path, 0.5 * (bracket.lower.t + bracket.upper.t));
        if (middle.solution.status != qp::Status::Optimal)
        {
            return stopped(data, path, middle.solution.status);
        }
        narrow(path, bracket, std::move(middle));
    }
    return stopped(data, path, qp::Status::MaxIterations);
}

Result climb(const Data& data, Path& path, PathPoint lower, double first_t)
{
    double t = first_t;
    for (int doubling = 0; doubling < kMaxDoublings; ++doubling)
    {
        PathPoint upper = solve_at(path, t);
        if (upper.solution.status != qp::Status::Optimal)
        {
            return stopped(data, path, upper.solution.status);
        }
        if (upper.volatility >= path.cap)
        {
            return bisect(data, path, Bracket{std::move(lower), std::move(upper)});
        }
        if (stays(lower, upper) && vertex_exit(lower, upper) == kInfinity)
        {
            return at_limit(data, path, lower, upper);
        }
        lower = std::move(upper);
        t *= 2.0;
    }
    return stopped(data, path, qp::Status::MaxIterations);
}

// t_b = sup{t : w(t) = w(0)} decides the multipliers at s = sigma_min: infinite, at a vertex exit, or zero.
// moved is the smallest t seen with w(t) != w(0), which bounds t_b.
Result leave_minimum(const Data& data, const Path& path, const PathPoint& start, const PathPoint& point,
                     double moved)
{
    double exit = vertex_exit(start, point);
    if (exit == kInfinity && moved == kInfinity)
    {
        return at_limit(data, path, start, point);
    }
    if (exit <= moved)
    {
        return at_exit(data, path, start, point, exit);
    }
    return vertical(data, path, start);
}

// Halves t from first_t until w(t) = w(0); when it never does, t_b = 0.
Result at_minimum(const Data& data, Path& path, const PathPoint& start, double first_t)
{
    double moved = kInfinity;
    double t = first_t;
    for (int step = 0; step < kMaxBisections; ++step)
    {
        PathPoint point = solve_at(path, t);
        if (point.solution.status != qp::Status::Optimal)
        {
            return stopped(data, path, point.solution.status);
        }
        if (stays(start, point))
        {
            return leave_minimum(data, path, start, point, moved);
        }
        moved = t;
        t *= 0.5;
    }
    return vertical(data, path, start);
}

// t has the units of a variance over a return. With a zero minimum volatility (all cash), s^2 sets the scale.
double first_step(std::span<const double> excess, double minimum, double cap)
{
    double variance = minimum * minimum;
    if (variance == 0.0)
    {
        variance = cap * cap;
    }
    return variance / std::max(largest_magnitude(excess), kEpsilon);
}

} // namespace

Result solve_max_return(const Spec& spec, const Data& data)
{
    Path path = make_path(spec, data);
    PathPoint start = solve_at(path, 0.0);
    if (start.solution.status != qp::Status::Optimal)
    {
        return stopped(data, path, start.solution.status);
    }
    double minimum = start.volatility;
    if (path.cap < minimum * (1.0 - kVolatilityTolerance))
    {
        return stopped(data, path, qp::Status::Infeasible);
    }
    double first_t = first_step(path.excess, minimum, path.cap);
    if (std::fabs(path.cap - minimum) <= kVolatilityTolerance * path.cap)
    {
        return at_minimum(data, path, start, first_t);
    }
    return climb(data, path, std::move(start), first_t);
}

} // namespace po::problems
