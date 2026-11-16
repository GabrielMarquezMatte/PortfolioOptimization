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

// Q(t): minimize 0.5 w'Sigma w - t v'w over the constraint rows of the spec, v = mu~ - shift 1 with
// mu~ = mu - rf 1, G = Sigma and the factor L of Data; only g = -t v depends on t. Where the rows hold 1'w = 1,
// shift is the midrange of mu~ and -t v'w = -t mu~'w + t shift: the same optimum w(t), whose QP budget
// multiplier is that of mu~ plus t shift. Elsewhere shift = 0. w(t) is piecewise affine in t with a
// nondecreasing volatility. On a segment the active set is fixed and w and the QP multipliers u are affine in
// t, with Sigma w - t v = N u for the active normals N. cap is the volatility cap s of MAX_RETURN.
struct Path
{
    Translation translation;
    // v
    std::vector<double> excess;
    std::span<const double> factor;
    double cap = 0.0;
    double shift = 0.0;
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

// Where w(t) is constant, du/dt = b with N b = -v (the user multipliers of the inactive cap).
struct Slope
{
    std::vector<double> equality;
    std::vector<double> inequality;
};

Path make_path(const Spec& spec, const Data& data);

// Solves Q(t) and adds its iterations to the path; the volatility is NaN unless the solution is Optimal.
PathPoint solve_at(Path& path, double t);

bool same_segment(const PathPoint& earlier, const PathPoint& later);

// The weights differ by more than the rounding of Q(t) can explain (see weight_noise in the source).
bool weights_moved(const Path& path, const PathPoint& earlier, const PathPoint& later);

// w(t) is the same at both points: the same active set, and either a vertex (the equalities and the active
// rows have rank n, so they fix w) or weights that have not moved.
bool weights_constant(const Path& path, const PathPoint& earlier, const PathPoint& later);

// N b = -v over the active normals of point, by least squares. The solver keeps the active normals linearly
// independent, so N has full column rank.
Slope multiplier_slope(const Path& path, const PathPoint& point);

// Where w(t) is constant from point on: the t where the first decreasing active multiplier reaches zero, or
// infinity when none decreases and w stays optimal for every larger t.
double vertex_exit(const Path& path, const PathPoint& point, const Slope& slope);

// The KKT conditions of max v'w over the rows (the user's, with a volatility multiplier of 0) at point, with
// the multipliers b of slope: every row holds to the feasibility tolerance of the QP solver, and N b = -v to the
// rounding of b. Together with b >= 0 on the active rows (vertex_exit infinite), they certify that w(point)
// maximizes mu~'w over the rows and solves Q(t) for every larger t.
bool limit_holds(const Path& path, const PathPoint& point, const Slope& slope);

// The weights of point with the multipliers u(point) + (t - point.t) b, along a stretch where w is constant.
qp::Solution moved_along(const PathPoint& point, const Slope& slope, double t);

// The weights of point with the multipliers b.
qp::Solution slope_solution(const PathPoint& point, const Slope& slope);

// The fraction tau in [0, 1] where |L'w| reaches s, with w = w_a + tau (w_b - w_a) on the bracket's segment.
double crossing_fraction(const Path& path, const Bracket& bracket);

// The point of the bracket's segment at fraction: weights and multipliers interpolated, an exact solution of
// Q(t_a + fraction (t_b - t_a)) since the KKT conditions on one active set are linear in (w, u, t).
qp::Solution on_segment(const Bracket& bracket, double fraction);

} // namespace po::problems
