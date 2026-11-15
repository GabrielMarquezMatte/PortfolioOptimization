#include "problems/parametric_path.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

#include "problems/mean_variance.hpp"
#include "qp/goldfarb_idnani.hpp"

namespace po::problems
{
namespace
{

constexpr double kSameWeightsTolerance = 1e-14;
constexpr double kInfinity = std::numeric_limits<double>::infinity();

double dot(std::span<const double> a, std::span<const double> b)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        sum += a[i] * b[i];
    }
    return sum;
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

bool same_weights(std::span<const double> earlier, std::span<const double> later)
{
    double scale = std::max(largest_magnitude(earlier), largest_magnitude(later));
    return largest_magnitude(difference(later, earlier)) <= kSameWeightsTolerance * scale;
}

} // namespace

double largest_magnitude(std::span<const double> values)
{
    double largest = 0.0;
    for (double value : values)
    {
        largest = std::max(largest, std::fabs(value));
    }
    return largest;
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
    PathPoint point{t, qp::solve(problem), std::numeric_limits<double>::quiet_NaN()};
    path.iterations += point.solution.iterations;
    if (point.solution.status == qp::Status::Optimal)
    {
        std::vector<double> root = factor_transposed_times(path.factor, point.solution.x);
        point.volatility = std::sqrt(dot(root, root));
    }
    return point;
}

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

bool stays(const PathPoint& earlier, const PathPoint& later)
{
    return earlier.solution.active_inequalities == later.solution.active_inequalities &&
           same_weights(earlier.solution.x, later.solution.x);
}

// |L'w|^2 = c + 2 b tau + a tau^2 reaches s^2 at the root in [0, 1], written without the cancellation of
// -b + sqrt(...) since b >= 0 where the volatility is nondecreasing.
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

} // namespace po::problems
