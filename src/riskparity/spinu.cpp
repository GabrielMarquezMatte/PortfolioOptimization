#include "riskparity/spinu.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <utility>

#include "core/portable_math.hpp"
#include "linalg/cholesky.hpp"

namespace po::riskparity
{
namespace
{

constexpr int kMaxIterations = 200;
constexpr double kStepFactor = 16.0;
constexpr double kFullStepRegion = 0.25;
constexpr double kStallRatio = 0.5;
constexpr double kStallCeiling = 1e-8;
constexpr double kBoundaryFraction = 0.99;
constexpr double kSufficientDecrease = 1e-4;
constexpr int kMaxHalvings = 60;

struct Problem
{
    std::span<const double> covariance;
    std::span<const double> budget;
    std::size_t n;
};

// F(y) = 1/2 y'Sigma y - b' log y at y, with Sigma y kept to build the gradient.
struct Point
{
    std::vector<double> y;
    std::vector<double> sigma_y;
    double value;
};

struct NewtonStep
{
    std::vector<double> direction;
    double slope;
};

double dot(std::span<const double> a, std::span<const double> b)
{
    double total = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        total += a[i] * b[i];
    }
    return total;
}

std::vector<double> covariance_times(const Problem& problem, std::span<const double> y)
{
    std::vector<double> product(problem.n);
    for (std::size_t i = 0; i < problem.n; ++i)
    {
        product[i] = dot(problem.covariance.subspan(i * problem.n, problem.n), y);
    }
    return product;
}

Point evaluate(const Problem& problem, std::vector<double> y)
{
    std::vector<double> sigma_y = covariance_times(problem, y);
    double value = 0.5 * dot(y, sigma_y);
    for (std::size_t i = 0; i < problem.n; ++i)
    {
        value -= problem.budget[i] * portable::log(y[i]);
    }
    return Point{std::move(y), std::move(sigma_y), value};
}

// y'Sigma y = sum b = 1 at the minimizer, and y_i proportional to sqrt(b_i) / sigma_i is its value for Sigma
// diagonal: the start is that vector scaled to y'Sigma y = 1.
Point starting_point(const Problem& problem)
{
    std::vector<double> y(problem.n);
    for (std::size_t i = 0; i < problem.n; ++i)
    {
        y[i] = std::sqrt(problem.budget[i] / problem.covariance[i * problem.n + i]);
    }
    double scale = 1.0 / std::sqrt(dot(y, covariance_times(problem, y)));
    for (double& entry : y)
    {
        entry *= scale;
    }
    return evaluate(problem, std::move(y));
}

std::vector<double> gradient_of(const Problem& problem, const Point& point)
{
    std::vector<double> gradient(problem.n);
    for (std::size_t i = 0; i < problem.n; ++i)
    {
        gradient[i] = point.sigma_y[i] - problem.budget[i] / point.y[i];
    }
    return gradient;
}

std::vector<double> hessian_of(const Problem& problem, const Point& point)
{
    std::vector<double> hessian(problem.covariance.begin(), problem.covariance.end());
    for (std::size_t i = 0; i < problem.n; ++i)
    {
        hessian[i * problem.n + i] += problem.budget[i] / point.y[i] / point.y[i];
    }
    return hessian;
}

std::optional<NewtonStep> newton_step(const Problem& problem, const Point& point)
{
    std::vector<double> gradient = gradient_of(problem, point);
    std::optional<std::vector<double>> factor = linalg::cholesky(hessian_of(problem, point), problem.n);
    std::optional<NewtonStep> step(std::in_place);
    if (!factor.has_value())
    {
        step.reset();
        return step;
    }
    step->direction.assign(problem.n, 0.0);
    for (std::size_t i = 0; i < problem.n; ++i)
    {
        step->direction[i] = -gradient[i];
    }
    linalg::solve_lower(*factor, problem.n, step->direction);
    linalg::solve_lower_transposed(*factor, problem.n, step->direction);
    step->slope = dot(gradient, step->direction);
    return step;
}

// The largest step that keeps y positive, shortened to a fraction of the way to the boundary and capped at
// the full Newton step.
double boundary_step(const Point& point, std::span<const double> direction)
{
    double longest = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < direction.size(); ++i)
    {
        if (direction[i] < 0.0)
        {
            longest = std::fmin(longest, -point.y[i] / direction[i]);
        }
    }
    return std::fmin(1.0, kBoundaryFraction * longest);
}

std::vector<double> moved(const Point& point, std::span<const double> direction, double step)
{
    std::vector<double> y(point.y.size());
    for (std::size_t i = 0; i < y.size(); ++i)
    {
        y[i] = point.y[i] + step * direction[i];
    }
    return y;
}

// Armijo backtracking; nullopt when no step of the sequence decreases F enough. The caller then takes the
// step to the boundary fraction anyway (advance), so this is not a stop.
std::optional<Point> armijo(const Problem& problem, const Point& point, const NewtonStep& step)
{
    double length = boundary_step(point, step.direction);
    for (int halving = 0; halving < kMaxHalvings; ++halving)
    {
        Point candidate = evaluate(problem, moved(point, step.direction, length));
        double sufficient = point.value + kSufficientDecrease * length * step.slope;
        if (candidate.value < point.value && candidate.value <= sufficient)
        {
            return candidate;
        }
        length *= 0.5;
    }
    return std::nullopt;
}

// max_i |d_i| / y_i: the relative change of every weight at once, so that a step of size delta leaves each
// weight with a relative error of order delta^2.
double relative_size(const Point& point, std::span<const double> direction)
{
    double largest = 0.0;
    for (std::size_t i = 0; i < direction.size(); ++i)
    {
        largest = std::fmax(largest, std::fabs(direction[i] / point.y[i]));
    }
    return largest;
}

// Close to the minimizer the full Newton step is taken without testing F, which no longer changes at the
// rounding of its value while the weights still do; far from it the step is damped by Armijo.
Point advance(const Problem& problem, const Point& point, const NewtonStep& step, double size)
{
    if (size <= kFullStepRegion)
    {
        return evaluate(problem, moved(point, step.direction, 1.0));
    }
    std::optional<Point> next = armijo(problem, point, step);
    if (next.has_value())
    {
        return std::move(*next);
    }
    return evaluate(problem, moved(point, step.direction, boundary_step(point, step.direction)));
}

std::vector<double> normalized(std::span<const double> y)
{
    double total = 0.0;
    for (double entry : y)
    {
        total += entry;
    }
    std::vector<double> weights(y.size());
    for (std::size_t i = 0; i < y.size(); ++i)
    {
        weights[i] = y[i] / total;
    }
    return weights;
}

} // namespace

SpinuResult spinu(const Data& data, std::span<const double> budget)
{
    Problem problem{data.covariance(), budget, data.asset_count()};
    Point point = starting_point(problem);
    int iterations = 0;
    bool converged = false;
    double previous = std::numeric_limits<double>::infinity();
    while (iterations < kMaxIterations)
    {
        std::optional<NewtonStep> step = newton_step(problem, point);
        if (!step.has_value())
        {
            break;
        }
        double size = relative_size(point, step->direction);
        point = advance(problem, point, *step, size);
        ++iterations;
        if (size <= kStepFactor * std::numeric_limits<double>::epsilon())
        {
            converged = true;
            break;
        }
        if (size >= kStallRatio * previous && size <= kStallCeiling)
        {
            converged = true;
            break;
        }
        previous = size;
    }
    return SpinuResult{normalized(point.y), iterations, converged};
}

} // namespace po::riskparity
