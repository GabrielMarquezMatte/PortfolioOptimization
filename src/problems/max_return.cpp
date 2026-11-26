#include "problems/max_return.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <utility>

#include "problems/constraints.hpp"
#include "problems/parametric_path.hpp"
#include "qp/problem.hpp"

namespace po::problems
{
namespace
{

constexpr int kMaxDoublings = 200;
constexpr int kMaxBisections = 200;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kInfinity = std::numeric_limits<double>::infinity();
constexpr double kEpsilon = std::numeric_limits<double>::epsilon();
// The rounding of a computed volatility |L'w|, as a multiple of n eps of it.
constexpr double kVolatilityNoise = 16.0;

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

// The point of the face where vol = s: there Sigma w - t* v = N u, so the user multipliers are u / t* and the
// volatility multiplier is -s / t*.
Result at_tangent(const Data& data, const Path& path, const Crossing& crossing)
{
    return answer(data, path, crossing.solution, crossing.t, -path.cap / crossing.t);
}

// vol(lower) < s <= vol(upper) with the two within the rounding of a volatility: vol = s to it at either end. The
// weights cannot tell: far along the path their rounding, t |v| n eps, covers points of different volatility.
bool volatility_at_cap(const Path& path, const Bracket& bracket)
{
    double n = static_cast<double>(bracket.upper.solution.x.size());
    return bracket.upper.volatility - bracket.lower.volatility <= kVolatilityNoise * n * kEpsilon * path.cap;
}

// t_a and t_b are adjacent doubles, with a breakpoint between them or Q(t) lost in its rounding: the crossing on the
// face of either end that certifies, the lower first (the upper first when the lower is t = 0). Where neither does,
// vol = s at both ends to its rounding and w(t) is constant on the face of the first (fixed by its rows, or
// N b = -v: no face line to cross on), the first end is the crossing, with the multipliers of at_exit there.
// Otherwise a breakdown.
Result from_either_end(const Data& data, const Path& path, const Bracket& bracket)
{
    const PathPoint* first = &bracket.lower;
    const PathPoint* second = &bracket.upper;
    if (first->t == 0.0)
    {
        std::swap(first, second);
    }
    for (const PathPoint* end : {first, second})
    {
        std::optional<Crossing> crossing = tangent_crossing(path, *end);
        if (crossing.has_value())
        {
            return at_tangent(data, path, *crossing);
        }
    }
    Slope slope = multiplier_slope(path, *first);
    if (volatility_at_cap(path, bracket) &&
        (fixed_by_face(path, *first) || limit_holds(path, *first, slope)))
    {
        return at_exit(data, path, *first, slope, first->t);
    }
    return broke_down(data, path);
}

// Whether the point became the lower end.
bool narrow(const Path& path, Bracket& bracket, PathPoint point)
{
    if (point.volatility >= path.cap)
    {
        bracket.upper = std::move(point);
        return false;
    }
    bracket.lower = std::move(point);
    return true;
}

// The Goldfarb-Idnani points decide the bracket. When both ends are on one face, the tangent of that face at the
// lower end gives the crossing; it is tried once per lower end, and one that does not certify leaves the bisection
// to go on.
Result bisect(const Data& data, Path& path, Bracket bracket)
{
    bool untried = true;
    for (int step = 0; step < kMaxBisections; ++step)
    {
        if (untried && same_segment(bracket.lower, bracket.upper))
        {
            std::optional<Crossing> crossing = tangent_crossing(path, bracket.lower);
            if (crossing.has_value())
            {
                return at_tangent(data, path, *crossing);
            }
            untried = false;
        }
        double middle = 0.5 * (bracket.lower.t + bracket.upper.t);
        if (middle <= bracket.lower.t || middle >= bracket.upper.t)
        {
            return from_either_end(data, path, bracket);
        }
        PathPoint point = solve_at(path, middle);
        if (point.solution.status != qp::Status::Optimal)
        {
            return broke_down(data, path);
        }
        if (narrow(path, bracket, std::move(point)))
        {
            untried = true;
        }
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
    if (path.cap < minimum * (1.0 - kMinimumVolatilityTolerance))
    {
        return stopped(data, path, qp::Status::Infeasible);
    }
    double first_t = first_step(data, path.base.excess, minimum, path.cap);
    if (std::fabs(path.cap - minimum) <= kMinimumVolatilityTolerance * path.cap)
    {
        return at_minimum(data, path, start, first_t);
    }
    return climb(data, path, std::move(start), first_t);
}

} // namespace po::problems
