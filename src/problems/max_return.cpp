#include "problems/max_return.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>

#include "problems/parametric_path.hpp"
#include "qp/problem.hpp"

namespace po::problems
{
namespace
{

// s within this of sigma_min, relative, is the case s = sigma_min.
constexpr double kMinimumTolerance = 1e-12;
// Q(t*) confirms the crossing when it lies on the bracket's segment and its volatility is within this
// fraction of the bracket's volatility range from s: an error of the affine model is a sizable fraction of
// the range, while the rounding of vol(Q(t*)) is about 1e-15 cond(Sigma) s (docs/guide/validation.md).
constexpr double kCrossingSanity = 1e-6;
// Across a breakpoint at adjacent doubles, the end is the answer when its volatility is within this of s,
// relative: 140 times the largest measured rounding of vol(Q(t)) (at cond(Sigma) = 1e8), and 160 times below
// the smallest deviation measured where Q(t) is lost in its rounding (docs/guide/validation.md).
constexpr double kAdjacentEnd = 1e-5;
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

// Feasibility is decided once, at t = 0: a QP that fails at t > 0 has broken down along the path.
Result broke_down(const Data& data, const Path& path)
{
    return stopped(data, path, qp::Status::MaxIterations);
}

// The path reports the multipliers of Q(t) whatever the shift a point was solved with: over multiplier_scale they
// are the user's.
Result answer(const Data& data, const Path& path, const qp::Solution& solution, double multiplier_scale,
              double volatility_multiplier)
{
    Result result = to_result(data, path.translation.rows, solution, multiplier_scale);
    result.volatility_multiplier = volatility_multiplier;
    result.iterations = path.iterations;
    return result;
}

// Where w(t) is constant, Sigma w - t v = N u(t) for every t, so -v = N b with b = du/dt: the user
// multipliers are b and the volatility multiplier is 0.
Result at_limit(const Data& data, const Path& path, const PathPoint& point, const Slope& slope)
{
    return answer(data, path, slope_solution(point, slope), 1.0, 0.0);
}

// w(t) stays at w(point) up to exit, where an active multiplier reaches zero. There
// Sigma w - exit v = N u(exit), so the user multipliers are u(exit) / exit and the volatility multiplier is
// -s / exit, the right derivative of the optimum in s.
Result at_exit(const Data& data, const Path& path, const PathPoint& point, const Slope& slope, double exit)
{
    return answer(data, path, moved_along(point, slope, exit), exit, -path.cap / exit);
}

// The volatility grows as soon as t > 0: at s = sigma_min the frontier is vertical, the volatility multiplier
// is -infinity and the other multipliers have no finite value.
Result vertical(const Data& data, const Path& path, const PathPoint& start)
{
    Result result = to_result(data, path.translation.rows, start.solution);
    result.budget_multiplier = kNaN;
    result.target_return_multiplier = kNaN;
    result.leverage_multiplier = kNaN;
    std::ranges::fill(result.bound_multipliers, kNaN);
    result.volatility_multiplier = -kInfinity;
    result.iterations = path.iterations;
    return result;
}

// The point of the segment where vol = s: there Sigma w - t* v = N u, so the user multipliers are u / t*
// and the volatility multiplier is -s / t*.
Result at_crossing(const Data& data, const Path& path, const Bracket& bracket, double fraction)
{
    double t = bracket.lower.t + fraction * (bracket.upper.t - bracket.lower.t);
    return answer(data, path, on_segment(bracket, fraction), t, -path.cap / t);
}

bool confirms(const Path& path, const Bracket& bracket, const PathPoint& crossing)
{
    double range = bracket.upper.volatility - bracket.lower.volatility;
    return same_segment(bracket.lower, crossing) &&
           std::fabs(crossing.volatility - path.cap) <= kCrossingSanity * range;
}

// t_a and t_b are adjacent doubles. On one segment the root is still exact; across a breakpoint, the lower
// end, whose volatility is below s by rounding only, with the volatility multiplier -vol / t that keeps it
// stationary (the upper end when the lower is t = 0). An end farther from s is Q(t) lost in its own rounding
// (vol(t) is continuous): the path has broken down.
Result at_adjacent(const Data& data, const Path& path, const Bracket& bracket)
{
    if (same_segment(bracket.lower, bracket.upper))
    {
        return at_crossing(data, path, bracket, crossing_fraction(path, bracket));
    }
    const PathPoint* end = &bracket.lower;
    if (end->t == 0.0)
    {
        end = &bracket.upper;
    }
    if (std::fabs(end->volatility - path.cap) > kAdjacentEnd * path.cap)
    {
        return broke_down(data, path);
    }
    return answer(data, path, end->solution, end->t, -end->volatility / end->t);
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

// On one segment w(t) is affine and the root of vol(w(t)) = s is exact; Q(t*) confirms that t* lies on that
// segment. A crossing that does not narrows the bracket like any other point.
Result bisect(const Data& data, Path& path, Bracket bracket)
{
    for (int step = 0; step < kMaxBisections; ++step)
    {
        if (same_segment(bracket.lower, bracket.upper))
        {
            double fraction = crossing_fraction(path, bracket);
            double t = bracket.lower.t + fraction * (bracket.upper.t - bracket.lower.t);
            PathPoint crossing = solve_at(path, t);
            if (crossing.solution.status != qp::Status::Optimal)
            {
                return broke_down(data, path);
            }
            if (confirms(path, bracket, crossing))
            {
                return at_crossing(data, path, bracket, fraction);
            }
            narrow(path, bracket, std::move(crossing));
        }
        double middle = 0.5 * (bracket.lower.t + bracket.upper.t);
        if (middle <= bracket.lower.t || middle >= bracket.upper.t)
        {
            return at_adjacent(data, path, bracket);
        }
        PathPoint point = solve_at(path, middle);
        if (point.solution.status != qp::Status::Optimal)
        {
            return broke_down(data, path);
        }
        narrow(path, bracket, std::move(point));
    }
    return broke_down(data, path);
}

// w(t) stays at w(point) for every larger t, certified at point by the KKT conditions of the limit: no active
// multiplier decreasing and limit_holds.
bool stays_beyond(const Path& path, const PathPoint& point, Slope& slope)
{
    slope = multiplier_slope(path, point);
    return vertex_exit(point, slope) == kInfinity && limit_holds(path, point, slope);
}

Result climb(const Data& data, Path& path, PathPoint lower, double first_t)
{
    double t = first_t;
    for (int doubling = 0; doubling < kMaxDoublings; ++doubling)
    {
        PathPoint upper = solve_at(path, t);
        if (upper.solution.status != qp::Status::Optimal)
        {
            return broke_down(data, path);
        }
        if (upper.volatility >= path.cap)
        {
            return bisect(data, path, Bracket{std::move(lower), std::move(upper)});
        }
        Slope slope;
        if (stays_beyond(path, upper, slope))
        {
            return at_limit(data, path, upper, slope);
        }
        lower = std::move(upper);
        t *= 2.0;
    }
    return broke_down(data, path);
}

// t_b = sup{t : w(t) = w(0)} decides the multipliers at s = sigma_min: zero when t_b is infinite, -s / t_b at
// a vertex exit, and -infinity when t_b = 0. moved is the smallest t seen on the segment of w(0) where the
// weights moved beyond their rounding, which bounds t_b: -infinity needs that evidence and no exit before it;
// otherwise the multipliers decide.
Result leave_minimum(const Data& data, const Path& path, const PathPoint& start, const PathPoint& point,
                     double moved)
{
    Slope slope = multiplier_slope(path, point);
    double exit = vertex_exit(point, slope);
    if (moved == kInfinity && exit == kInfinity && !limit_holds(path, point, slope))
    {
        return broke_down(data, path);
    }
    if (moved == kInfinity && exit == kInfinity)
    {
        return at_limit(data, path, point, slope);
    }
    if (exit <= moved)
    {
        return at_exit(data, path, point, slope, exit);
    }
    return vertical(data, path, start);
}

// Halves t from first_t until w(t) = w(0).
Result at_minimum(const Data& data, Path& path, const PathPoint& start, double first_t)
{
    double moved = kInfinity;
    double t = first_t;
    for (int step = 0; step < kMaxBisections; ++step)
    {
        PathPoint point = solve_at(path, t);
        if (point.solution.status != qp::Status::Optimal)
        {
            return broke_down(data, path);
        }
        if (weights_constant(path, start, point))
        {
            return leave_minimum(data, path, start, point, moved);
        }
        if (same_segment(start, point) && weights_moved(start, point))
        {
            moved = t;
        }
        t *= 0.5;
    }
    if (moved == kInfinity)
    {
        return broke_down(data, path);
    }
    return vertical(data, path, start);
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

// t has the units of a variance over a return: sigma_min^2, or, with a zero minimum volatility (all cash),
// the smaller of s^2 and the largest asset variance. The rounding of Q(t) grows with t, so a large cap must
// not start the path far beyond its last breakpoint.
double first_step(const Data& data, std::span<const double> excess, double minimum, double cap)
{
    double variance = minimum * minimum;
    if (variance == 0.0)
    {
        std::size_t n = data.asset_count();
        double largest = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            largest = std::max(largest, data.covariance()[i * n + i]);
        }
        variance = std::min(cap * cap, largest);
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
    if (path.cap < minimum * (1.0 - kMinimumTolerance))
    {
        return stopped(data, path, qp::Status::Infeasible);
    }
    double first_t = first_step(data, path.base.excess, minimum, path.cap);
    if (std::fabs(path.cap - minimum) <= kMinimumTolerance * path.cap)
    {
        return at_minimum(data, path, start, first_t);
    }
    return climb(data, path, std::move(start), first_t);
}

} // namespace po::problems
