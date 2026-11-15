#pragma once

#include <limits>
#include <span>
#include <vector>

#include "core/data.hpp"
#include "core/spec.hpp"
#include "problems/constraints.hpp"
#include "qp/problem.hpp"

namespace po::problems
{

// Q(t): minimize 0.5 w'Sigma w - t (mu - rf 1)'w over the constraint rows of the spec, with G = Sigma and the
// factor L of Data; only g = -t (mu - rf 1) depends on t. Its optimum w(t) is piecewise affine in t with a
// nondecreasing volatility; cap is the volatility cap s of MAX_RETURN.
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
    double volatility = std::numeric_limits<double>::quiet_NaN();
};

// vol(w(lower.t)) < s <= vol(w(upper.t)).
struct Bracket
{
    PathPoint lower;
    PathPoint upper;
};

Path make_path(const Spec& spec, const Data& data);

// Solves Q(t) and adds its iterations to the path; the volatility is NaN unless the solution is Optimal.
PathPoint solve_at(Path& path, double t);

// The weights of earlier with the multipliers a u(earlier) + b u(later). An inequality multiplier is
// nonnegative, and decreases() keeps every combination taken here at or above zero up to rounding, which is
// dropped.
qp::Solution combined(const PathPoint& earlier, const PathPoint& later, double a, double b);

// With u(t) affine through both points, the t where the first decreasing active multiplier reaches zero;
// infinity when none decreases, so that w stays optimal for every larger t.
double vertex_exit(const PathPoint& earlier, const PathPoint& later);

// Same active set and the same weights at both points.
bool stays(const PathPoint& earlier, const PathPoint& later);

// With w = w_a + tau (w_b - w_a), the t where |L'w| reaches s on the segment of the bracket.
double crossing_t(const Path& path, const Bracket& bracket);

// The largest magnitude of the entries, 0 for none.
double largest_magnitude(std::span<const double> values);

} // namespace po::problems
