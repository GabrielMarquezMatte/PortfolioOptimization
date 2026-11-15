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

constexpr double kAtMinimumTolerance = 1e-12;
constexpr double kSameWeightsTolerance = 1e-14;
constexpr int kMaxDoublings = 200;
constexpr int kMaxBisections = 200;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kEpsilon = std::numeric_limits<double>::epsilon();

// Q(t): G = Sigma with the factor L of Data, and g = -t (mu - rf 1); only g depends on t.
struct Path
{
    Translation translation;
    std::vector<double> excess;
    std::span<const double> factor;
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
                data.covariance_factor(), 0};
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

// Where w(t) is constant, the QP multipliers are u(t) = a + b t and Sigma w - t (mu - rf 1) = N u(t) for
// every t, so -(mu - rf 1) = N b: the user multipliers are the slope b, and the volatility multiplier is 0.
Result at_limit(const Data& data, const Path& path, const PathPoint& earlier, const PathPoint& later)
{
    qp::Solution slope = earlier.solution;
    slope.equality_multipliers =
        difference(later.solution.equality_multipliers, earlier.solution.equality_multipliers);
    slope.inequality_multipliers =
        difference(later.solution.inequality_multipliers, earlier.solution.inequality_multipliers);
    return answer(data, path, slope, later.t - earlier.t, 0.0);
}

bool same_weights(std::span<const double> earlier, std::span<const double> later)
{
    double scale = std::max(largest_magnitude(earlier), largest_magnitude(later));
    return largest_magnitude(difference(later, earlier)) <= kSameWeightsTolerance * scale;
}

// A multiplier that decreases reaches zero at a larger t, where w(t) moves again.
bool multipliers_grow(const PathPoint& earlier, const PathPoint& later)
{
    const std::vector<double>& before = earlier.solution.inequality_multipliers;
    const std::vector<double>& after = later.solution.inequality_multipliers;
    return std::ranges::all_of(later.solution.active_inequalities,
                               [&](std::size_t row) { return after[row] >= before[row]; });
}

bool is_last_segment(const PathPoint& earlier, const PathPoint& later)
{
    return earlier.solution.active_inequalities == later.solution.active_inequalities &&
           same_weights(earlier.solution.x, later.solution.x) && multipliers_grow(earlier, later);
}

// With w = w_a + tau (w_b - w_a), |L'w|^2 = c + 2 b tau + a tau^2 reaches s^2 at the root in [0, 1], written
// without the cancellation of -b + sqrt(...) since b >= 0 where the volatility is nondecreasing.
double crossing_fraction(std::span<const double> factor, const Bracket& bracket, double cap)
{
    std::vector<double> start = factor_transposed_times(factor, bracket.lower.solution.x);
    std::vector<double> direction =
        factor_transposed_times(factor, difference(bracket.upper.solution.x, bracket.lower.solution.x));
    double curvature = dot(direction, direction);
    double slope = dot(start, direction);
    double gap = (cap - bracket.lower.volatility) * (cap + bracket.lower.volatility);
    return gap / (slope + std::sqrt(slope * slope + curvature * gap));
}

Result cross(const Data& data, Path& path, const Bracket& bracket, double cap)
{
    double width = bracket.upper.t - bracket.lower.t;
    double t = bracket.lower.t + crossing_fraction(path.factor, bracket, cap) * width;
    PathPoint point = solve_at(path, t);
    if (point.solution.status != qp::Status::Optimal)
    {
        return stopped(data, path, point.solution.status);
    }
    return answer(data, path, point.solution, t, -cap / t);
}

bool same_active_set(const Bracket& bracket)
{
    return bracket.lower.solution.active_inequalities == bracket.upper.solution.active_inequalities;
}

Result bisect(const Data& data, Path& path, Bracket bracket, double cap)
{
    for (int step = 0; !same_active_set(bracket); ++step)
    {
        if (step == kMaxBisections)
        {
            return stopped(data, path, qp::Status::MaxIterations);
        }
        PathPoint middle = solve_at(path, 0.5 * (bracket.lower.t + bracket.upper.t));
        if (middle.solution.status != qp::Status::Optimal)
        {
            return stopped(data, path, middle.solution.status);
        }
        if (middle.volatility >= cap)
        {
            bracket.upper = std::move(middle);
            continue;
        }
        bracket.lower = std::move(middle);
    }
    return cross(data, path, bracket, cap);
}

Result climb(const Data& data, Path& path, PathPoint lower, double first_t, double cap)
{
    double t = first_t;
    for (int doubling = 0; doubling < kMaxDoublings; ++doubling)
    {
        PathPoint upper = solve_at(path, t);
        if (upper.solution.status != qp::Status::Optimal)
        {
            return stopped(data, path, upper.solution.status);
        }
        if (upper.volatility >= cap)
        {
            return bisect(data, path, Bracket{std::move(lower), std::move(upper)}, cap);
        }
        if (is_last_segment(lower, upper))
        {
            return at_limit(data, path, lower, upper);
        }
        lower = std::move(upper);
        t *= 2.0;
    }
    return stopped(data, path, qp::Status::MaxIterations);
}

Result at_minimum(const Data& data, Path& path, const PathPoint& start, double first_t)
{
    PathPoint next = solve_at(path, first_t);
    if (next.solution.status != qp::Status::Optimal)
    {
        return stopped(data, path, next.solution.status);
    }
    return at_limit(data, path, start, next);
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
    double cap = *spec.max_volatility;
    double minimum = start.volatility;
    if (cap < minimum * (1.0 - kAtMinimumTolerance))
    {
        return stopped(data, path, qp::Status::Infeasible);
    }
    double first_t = first_step(path.excess, minimum, cap);
    if (std::fabs(cap - minimum) <= kAtMinimumTolerance * cap)
    {
        return at_minimum(data, path, start, first_t);
    }
    return climb(data, path, std::move(start), first_t, cap);
}

} // namespace po::problems
