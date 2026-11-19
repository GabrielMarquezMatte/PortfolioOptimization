#pragma once

#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <vector>

#include "core/data.hpp"
#include "core/spec.hpp"
#include "problems/constraints.hpp"
#include "qp/goldfarb_idnani.hpp"
#include "qp/problem.hpp"

namespace po::problems
{

// A row r'w >= rho (or = rho) whose normal has every entry -1 or 1: the budget, a cash bound, or a generated row
// s'w <= L of the leverage. It holds the whole portfolio, so it can carry the common mean of mu~.
struct Face
{
    enum class Kind
    {
        None,
        Equality,
        Inequality,
        Generated
    };
    Kind kind = Kind::None;
    // The equality or inequality row.
    std::size_t row = 0;
    // The sign vector s of a generated row a - s.
    std::vector<double> signs;
};

// Q(t) solved with v = mu~ - amount r for the face r'w >= rho, held as an equality, with amount the midrange of
// the entries r_i mu~_i: -t v'w = -t mu~'w + t amount rho on the face, so the optimum of Q(t) restricted to the face
// is the same, and its face multiplier in Q(t) is that of v less t amount. The solver starts from t Sigma^-1 v, and
// a common mean far above the spread of mu~ would leave its rounding in every weight. An equality face (the budget)
// always holds, so the shifted solution is that of Q(t). An inequality face is certified a posteriori: the
// restricted optimum solves Q(t) exactly when its face multiplier in Q(t) is nonnegative (the KKT conditions of
// Q(t) with the face among the active rows). Kind None is Q(t) itself.
struct Shift
{
    Face face;
    double amount = 0.0;
    // v
    std::vector<double> excess;
};

// Q(t): minimize 0.5 w'Sigma w - t mu~'w over the constraint rows of the spec and its gross exposure family, with
// mu~ = mu - rf 1, G = Sigma and the factor L of Data; only g = -t v depends on t, for the v of the shift each
// point is solved with. w(t) is piecewise affine in t with a nondecreasing volatility. On a segment the active set
// is fixed and w and the QP multipliers u of Q(t) are affine in t, with Sigma w - t mu~ = N u for the active
// normals N: the active explicit rows and the active generated rows a - s_j of the family, ordinary columns of N.
// The segment is the face, not its sign vectors: with the family active, the s_j agree with sign(w) on the support
// and w_i = 0 where two of them disagree (the zeros D), and the same face is held by several sets of s_j, between
// which the solver alternates under rounding. A weight that crosses zero under the leverage changes the face, and
// with it the segment. base is the shift of the budget where the rows hold 1'w = 1 (Kind None otherwise); without a
// budget, face is the last inequality face certified, tried first at the next point. cap is the volatility cap s of
// MAX_RETURN.
struct Path
{
    Translation translation;
    // mu~
    std::vector<double> excess;
    std::span<const double> factor;
    double cap = 0.0;
    Shift base;
    Shift face;
    int iterations = 0;
    // Every point is solved with the factor of Data, so J0 is built once.
    qp::BasisCache basis;
};

// solution holds the multipliers of Q(t) whatever the shift it was solved with; shift is that shift.
struct PathPoint
{
    double t = 0.0;
    qp::Solution solution;
    double volatility = std::numeric_limits<double>::quiet_NaN();
    Shift shift;
};

// vol(w(lower.t)) < s <= vol(w(upper.t)).
struct Bracket
{
    PathPoint lower;
    PathPoint upper;
};

// Where w(t) is constant, du/dt = b with N b = -mu~ (the user multipliers of the inactive cap), held as b_v with
// N b_v = -v in the frame of the point's shift: b is b_v less the shift amount at the face, a large common mean
// added to small entries, so b_v is kept and the functions below map it. gross_exposure is that of the active
// generated rows, in the order of the point's sign vectors.
struct Slope
{
    std::vector<double> equality;
    std::vector<double> inequality;
    std::vector<double> gross_exposure;
};

Path make_path(const Spec& spec, const Data& data);

// Solves Q(t) with the face of the path, then with its base shift and the face of the most loaded active
// inequality face, and adds the iterations to the path; the volatility is NaN unless the solution is Optimal.
PathPoint solve_at(Path& path, double t);

// The same face: the same active explicit rows, the family active at both or at neither and, where it is, the same
// support with the same signs and the same zeros (the zeros D where the sign vectors disagree, and the weights that
// are 0). The sign vectors themselves are not compared.
bool same_segment(const PathPoint& earlier, const PathPoint& later);

// The weights differ by more than the rounding of Q(t) can explain (see weight_noise in the source).
bool weights_moved(const PathPoint& earlier, const PathPoint& later);

// w(t) is the same at both points: the same face, and either a vertex (the equalities and the active
// rows have rank n, so they fix w) or weights that have not moved.
bool weights_constant(const Path& path, const PathPoint& earlier, const PathPoint& later);

// N b_v = -v over the active normals of point, for the v of its shift, by least squares. The solver keeps the
// active normals linearly independent, so N has full column rank.
Slope multiplier_slope(const Path& path, const PathPoint& point);

// Where w(t) is constant from point on: the t where the first decreasing active multiplier reaches zero, or
// infinity when none decreases and w stays optimal for every larger t.
double vertex_exit(const PathPoint& point, const Slope& slope);

// The KKT conditions of max mu~'w over the rows and the family (the user's, with a volatility multiplier of 0) at
// point, with the multipliers b of slope: every row and the family hold to the feasibility tolerance of the QP
// solver, and N b_v = -v, the generated normals among N, to the rounding of b_v. Together with b >= 0 on the active rows (vertex_exit infinite), they certify that w(point)
// maximizes mu~'w over the rows and solves Q(t) for every larger t.
bool limit_holds(const Path& path, const PathPoint& point, const Slope& slope);

// The weights of point with the multipliers u(point) + (t - point.t) b, along a stretch where w is constant.
qp::Solution moved_along(const PathPoint& point, const Slope& slope, double t);

// The weights of point with the multipliers b.
qp::Solution slope_solution(const PathPoint& point, const Slope& slope);

// A solution of Q(t) at the crossing vol(w(t)) = cap, multipliers in the convention of PathPoint::solution.
struct Crossing
{
    double t = 0.0;
    qp::Solution solution;
};

// From the face and shift of point: the exact tangent dw on that face (null space of the active normals over the
// free coordinates, QR of L'Z), anchored at a face solve at point.t; the root of ||L'(w_a + h dw)|| = cap by the
// cancellation-free formula; certified by the KKT conditions of Q(t*) at one point, Q(t*) being convex
// (docs/research/2026-11-18-max-return-path-and-portable-generators.md, section 2). nullopt when the certificate
// fails or the reduced Hessian has a pivot below its noise.
std::optional<Crossing> tangent_crossing(const Path& path, const PathPoint& point);

} // namespace po::problems
