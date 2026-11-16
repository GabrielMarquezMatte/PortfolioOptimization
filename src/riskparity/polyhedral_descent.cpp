#include "riskparity/polyhedral_descent.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

#include "problems/constraints.hpp"
#include "problems/mean_variance.hpp"
#include "qp/goldfarb_idnani.hpp"
#include "qp/problem.hpp"
#include "riskparity/objective.hpp"

namespace po::riskparity
{
namespace
{

constexpr double kSufficientDecrease = 1e-4;
constexpr int kMaxHalvings = 60;

// L = f + L3 at weights. Without a volatility term, volatility_gradient is empty and gradient is grad f.
struct Iterate
{
    std::vector<double> weights;
    Linearization linearization;
    std::vector<double> gradient;
    double volatility = 0.0;
    // grad h3 = Sigma w / vol(w).
    std::vector<double> volatility_gradient;
    // dL3 / dh3 = -multiplier + penalty h3.
    double volatility_slope = 0.0;
};

// The model as a QP in the step d = x - w, whose Hessian, linear term and right-hand sides are replaced at
// every step: the rows of the polyhedron in w (rhs kept here) shifted to the iterate. tau is the weight of
// the proximal term.
struct Descent
{
    const DescentProblem& problem;
    problems::Translation translation;
    std::vector<double> equality_rhs;
    std::vector<double> inequality_rhs;
    double tau;
};

// N'u = E'y + C'z of the QP solution, and the largest |u_k| max_j |N_kj| of its terms.
struct NormalTerm
{
    std::vector<double> total;
    double scale = 0.0;
};

// The QP solved at an iterate, whose x is the step d to the minimizer of the model, with N'u and d'Qd.
struct Step
{
    qp::Solution solution;
    NormalTerm normal;
    double curvature = 0.0;
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

double largest_magnitude(std::span<const double> values)
{
    double largest = 0.0;
    for (double value : values)
    {
        largest = std::fmax(largest, std::fabs(value));
    }
    return largest;
}

std::vector<double> matrix_times(std::span<const double> matrix, std::span<const double> x)
{
    std::size_t n = x.size();
    std::vector<double> product(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        product[i] = dot(matrix.subspan(i * n, n), x);
    }
    return product;
}

Descent prepared(const DescentProblem& problem, const DescentOptions& options)
{
    std::size_t n = problem.data.asset_count();
    double trace = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        trace += problem.data.covariance()[i * n + i];
    }
    problems::MeanVariance rows_only{2.0, 0.0};
    Descent descent{problem, problems::mean_variance_problem(problem.spec, problem.data, rows_only), {}, {},
                    options.proximal * trace / (2.0 * static_cast<double>(n))};
    descent.equality_rhs = descent.translation.problem.equality_rhs;
    descent.inequality_rhs = descent.translation.problem.inequality_rhs;
    return descent;
}

void add_volatility_term(const DescentProblem& problem, Iterate& iterate)
{
    if (!problem.volatility.has_value())
    {
        return;
    }
    const VolatilityTerm& term = *problem.volatility;
    iterate.volatility_gradient = matrix_times(problem.data.covariance(), iterate.weights);
    iterate.volatility = std::sqrt(dot(iterate.weights, iterate.volatility_gradient));
    for (double& entry : iterate.volatility_gradient)
    {
        entry /= iterate.volatility;
    }
    iterate.volatility_slope = -term.multiplier + term.penalty * (iterate.volatility - term.target);
    for (std::size_t i = 0; i < iterate.gradient.size(); ++i)
    {
        iterate.gradient[i] += iterate.volatility_slope * iterate.volatility_gradient[i];
    }
}

Iterate evaluate(const DescentProblem& problem, std::vector<double> weights)
{
    Iterate iterate;
    iterate.linearization = linearize(problem.data, problem.budget, weights, false);
    iterate.gradient = iterate.linearization.gradient;
    iterate.weights = std::move(weights);
    add_volatility_term(problem, iterate);
    return iterate;
}

// rhs - N w for every row: the polyhedron in the step d.
void shift_rows(std::span<const double> rows, std::span<const double> rhs, std::span<const double> weights,
                std::vector<double>& shifted)
{
    std::size_t n = weights.size();
    for (std::size_t k = 0; k < rhs.size(); ++k)
    {
        shifted[k] = rhs[k] - dot(rows.subspan(k * n, n), weights);
    }
}

// min 1/2 d'Qd + grad L'd over the shifted rows, Q = 2 J'J + tau I (+ penalty grad h3 grad h3'): the QP
// min 1/2 x'Qx + (grad L - Q w)'x over the polyhedron in x = w + d, posed in d so that d is not the rounding
// of the difference of two points of the size of w. The solver factors Q, which tau makes positive definite:
// J w = 0, so J'J is singular.
void set_model(Descent& descent, const Iterate& iterate)
{
    qp::Problem& model = descent.translation.problem;
    std::size_t n = model.n;
    model.hessian = gauss_newton(iterate.linearization.jacobian, n);
    model.hessian_factor.clear();
    for (std::size_t i = 0; i < n; ++i)
    {
        model.hessian[i * n + i] += descent.tau;
    }
    if (descent.problem.volatility.has_value())
    {
        double penalty = descent.problem.volatility->penalty;
        std::span<const double> gradient = iterate.volatility_gradient;
        for (std::size_t k = 0; k < n * n; ++k)
        {
            model.hessian[k] += penalty * gradient[k / n] * gradient[k % n];
        }
    }
    model.linear = iterate.gradient;
    shift_rows(model.equality_rows, descent.equality_rhs, iterate.weights, model.equality_rhs);
    shift_rows(model.inequality_rows, descent.inequality_rhs, iterate.weights, model.inequality_rhs);
}

// Adds sum_k u_k row_k to normal.
void add_rows(std::span<const double> rows, std::span<const double> multipliers, NormalTerm& normal)
{
    std::size_t n = normal.total.size();
    for (std::size_t k = 0; k < multipliers.size(); ++k)
    {
        std::span<const double> row = rows.subspan(k * n, n);
        for (std::size_t j = 0; j < n; ++j)
        {
            normal.total[j] += multipliers[k] * row[j];
        }
        normal.scale = std::fmax(normal.scale, std::fabs(multipliers[k]) * largest_magnitude(row));
    }
}

Step step_at(Descent& descent, const Iterate& iterate)
{
    set_model(descent, iterate);
    const qp::Problem& model = descent.translation.problem;
    Step step{qp::solve(model), NormalTerm{std::vector<double>(model.n, 0.0)}, 0.0};
    add_rows(model.equality_rows, step.solution.equality_multipliers, step.normal);
    add_rows(model.inequality_rows, step.solution.inequality_multipliers, step.normal);
    step.curvature = dot(step.solution.x, matrix_times(model.hessian, step.solution.x));
    return step;
}

// The largest |grad f_j| a gap of at most 1 in every entry can produce, 2 max_j sum_i |J_ij|: the scale of
// the stationarity at f = 0, where every term of grad f vanishes with the gap.
double unit_gap_scale(std::span<const double> jacobian, std::size_t n)
{
    double largest = 0.0;
    for (std::size_t j = 0; j < n; ++j)
    {
        double column = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            column += std::fabs(jacobian[i * n + j]);
        }
        largest = std::fmax(largest, 2.0 * column);
    }
    return largest;
}

// |grad L - N'u| relative to the largest term of grad f, of the volatility term and of N'u.
double stationarity_of(const Iterate& iterate, const Step& step, double tolerance)
{
    std::size_t n = iterate.weights.size();
    std::vector<double> residual(n);
    for (std::size_t j = 0; j < n; ++j)
    {
        residual[j] = iterate.gradient[j] - step.normal.total[j];
    }
    double scale = std::fmax(step.normal.scale, largest_magnitude(iterate.linearization.gradient));
    double volatility_scale = largest_magnitude(iterate.volatility_gradient);
    scale = std::fmax(scale, std::fabs(iterate.volatility_slope) * volatility_scale);
    if (std::sqrt(iterate.linearization.value) <= tolerance)
    {
        scale = std::fmax(scale, unit_gap_scale(iterate.linearization.jacobian, n));
    }
    double largest = largest_magnitude(residual);
    if (scale == 0.0)
    {
        return largest;
    }
    return largest / scale;
}

// The slack at weights of every inequality row with a nonzero multiplier, relative to max(|rhs|, |w|max). The
// shifted right-hand side of the model is minus that slack.
double complementarity_of(const Descent& descent, const Iterate& iterate, const qp::Solution& solution)
{
    std::span<const double> shifted = descent.translation.problem.inequality_rhs;
    double weight_scale = largest_magnitude(iterate.weights);
    double largest = 0.0;
    for (std::size_t k = 0; k < shifted.size(); ++k)
    {
        if (solution.inequality_multipliers[k] == 0.0)
        {
            continue;
        }
        double row_scale = std::fmax(std::fabs(descent.inequality_rhs[k]), weight_scale);
        largest = std::fmax(largest, std::fabs(shifted[k]) / row_scale);
    }
    return largest;
}

// L(w + s) - L(w), from s: dh3 = dv / (vol + vol') and dL3 = dh3 (dL3/dh3 + penalty dh3 / 2).
double merit_change(const DescentProblem& problem, const Iterate& iterate, std::span<const double> step)
{
    Change change = value_change(problem.data, problem.budget, iterate.weights, step);
    if (!problem.volatility.has_value())
    {
        return change.value;
    }
    double updated = std::sqrt(iterate.volatility * iterate.volatility + change.variance);
    double volatility_change = change.variance / (iterate.volatility + updated);
    double penalty = problem.volatility->penalty;
    return change.value + volatility_change * (iterate.volatility_slope + 0.5 * penalty * volatility_change);
}

std::vector<double> scaled(std::span<const double> direction, double length)
{
    std::vector<double> step(direction.begin(), direction.end());
    for (double& entry : step)
    {
        entry *= length;
    }
    return step;
}

std::vector<double> moved(std::span<const double> weights, std::span<const double> step)
{
    std::vector<double> point(weights.begin(), weights.end());
    for (std::size_t i = 0; i < point.size(); ++i)
    {
        point[i] += step[i];
    }
    return point;
}

// Armijo on the Lagrangian l(w) = L(w) - u'(N w - rhs) of the step's multipliers, from the full step down.
// Its slope along d is (grad L - N'u)'d = -d'Qd by the stationarity of the QP, and on the polyhedron
// L(w + s) - L(w) = l(w + s) - l(w) + u'N s <= l(w + s) - l(w), since E s = 0 and z_k C_k s <= 0 for the
// rows active at w + d, so every accepted step decreases L. The test on l rather than on L is what keeps it
// deciding near a stationary point: the rounding of the step across the rows moves L by u'N s, which there is
// far larger than the decrease, and l does not see it. nullopt when d'Qd is not positive or no halving
// decreases l enough.
std::optional<Iterate> armijo(const DescentProblem& problem, const Iterate& iterate, const Step& step)
{
    std::optional<Iterate> next;
    double slope = -step.curvature;
    if (!(slope < 0.0))
    {
        return next;
    }
    double length = 1.0;
    for (int halving = 0; halving < kMaxHalvings; ++halving)
    {
        std::vector<double> trial = scaled(step.solution.x, length);
        double change = merit_change(problem, iterate, trial) - dot(step.normal.total, trial);
        if (change <= kSufficientDecrease * length * slope)
        {
            next.emplace(evaluate(problem, moved(iterate.weights, trial)));
            return next;
        }
        length *= 0.5;
    }
    return next;
}

void record(const Descent& descent, const Iterate& iterate, const Step& step, double tolerance,
            DescentResult& result)
{
    result.equality_multipliers = step.solution.equality_multipliers;
    result.bound_multipliers =
        problems::to_result(descent.problem.data, descent.translation.rows, step.solution).bound_multipliers;
    result.stationarity = std::numeric_limits<double>::quiet_NaN();
    if (step.solution.status != qp::Status::Optimal)
    {
        return;
    }
    double complementarity = complementarity_of(descent, iterate, step.solution);
    result.stationarity = std::fmax(stationarity_of(iterate, step, tolerance), complementarity);
}

} // namespace

DescentResult descend(const DescentProblem& problem, std::span<const double> start,
                      const DescentOptions& options)
{
    Descent descent = prepared(problem, options);
    Iterate iterate = evaluate(problem, std::vector<double>(start.begin(), start.end()));
    DescentResult result{{}, {}, {}, std::numeric_limits<double>::quiet_NaN(), 0, false};
    while (true)
    {
        Step step = step_at(descent, iterate);
        record(descent, iterate, step, options.tolerance, result);
        if (std::isnan(result.stationarity))
        {
            break;
        }
        if (result.stationarity <= options.tolerance)
        {
            result.converged = true;
            break;
        }
        if (result.iterations >= options.max_iterations)
        {
            break;
        }
        std::optional<Iterate> next = armijo(problem, iterate, step);
        if (!next.has_value())
        {
            break;
        }
        iterate = std::move(*next);
        ++result.iterations;
    }
    result.weights = std::move(iterate.weights);
    return result;
}

} // namespace po::riskparity
