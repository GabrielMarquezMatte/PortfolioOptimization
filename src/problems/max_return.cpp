#include "problems/max_return.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include "problems/parametric_path.hpp"
#include "qp/problem.hpp"

namespace po::problems
{
namespace
{

constexpr double kVolatilityTolerance = 1e-12;
constexpr int kMaxDoublings = 200;
constexpr int kMaxBisections = 200;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInfinity = std::numeric_limits<double>::infinity();
constexpr double kEpsilon = std::numeric_limits<double>::epsilon();

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
