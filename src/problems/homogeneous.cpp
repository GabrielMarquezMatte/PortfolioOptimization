#include "problems/homogeneous.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include "linalg/cholesky.hpp"
#include "problems/constraints.hpp"
#include "qp/goldfarb_idnani.hpp"
#include "qp/problem.hpp"

namespace po::problems
{
namespace
{

// kappa at most this multiple of its rounding (kappa_noise) is the closure point kappa = 0: 400 times the
// largest error of kappa measured at the boundary, and 1.6e5 times below the smallest kappa of an OPTIMAL result
// of the battery, both in units of kappa_noise (docs/guide/validation.md).
constexpr double kKappaNoiseFactor = 16.0;
constexpr double kEpsilon = std::numeric_limits<double>::epsilon();

// The QP rows whose multipliers map back to the user's: the bound rows and the row 1'y >= 0.
struct Homogeneous
{
    qp::Problem problem;
    std::vector<BoundRow> bounds;
    std::size_t kappa_row = 0;
};

double bound_of(const std::vector<double>& per_weight, double scalar, std::size_t weight)
{
    if (per_weight.empty())
    {
        return scalar;
    }
    return per_weight[weight];
}

double lower_of(const Spec& spec, std::size_t weight)
{
    return bound_of(spec.asset_lower, spec.lower, weight);
}

double upper_of(const Spec& spec, std::size_t weight)
{
    return bound_of(spec.asset_upper, spec.upper, weight);
}

std::size_t add_row(qp::Problem& problem, std::span<const double> row)
{
    problem.inequality_rows.insert(problem.inequality_rows.end(), row.begin(), row.end());
    problem.inequality_rhs.push_back(0.0);
    return problem.inequality_rhs.size() - 1;
}

// scale 1'y + sign y_weight >= 0: the lower row y_i - lo_i 1'y (scale -lo_i, sign 1) or the upper row
// hi_i 1'y - y_i (scale hi_i, sign -1).
std::vector<double> bound_row(std::size_t n, std::size_t weight, double scale, double sign)
{
    std::vector<double> row(n, scale);
    row[weight] += sign;
    return row;
}

void add_bound_rows(const Spec& spec, Homogeneous& homogeneous)
{
    std::size_t n = homogeneous.problem.n;
    for (std::size_t i = 0; i < n; ++i)
    {
        double lower = lower_of(spec, i);
        double upper = upper_of(spec, i);
        if (std::isfinite(lower))
        {
            std::size_t row = add_row(homogeneous.problem, bound_row(n, i, -lower, 1.0));
            homogeneous.bounds.push_back(BoundRow{i, row, true});
        }
        if (std::isfinite(upper))
        {
            std::size_t row = add_row(homogeneous.problem, bound_row(n, i, upper, -1.0));
            homogeneous.bounds.push_back(BoundRow{i, row, false});
        }
    }
}

// G = 2 Sigma with the factor sqrt(2) L of Data, so that positive definiteness stays the one decision of
// Data::create.
void add_curvature(const Data& data, qp::Problem& problem)
{
    problem.hessian.assign(data.covariance().begin(), data.covariance().end());
    for (double& entry : problem.hessian)
    {
        entry *= 2.0;
    }
    problem.hessian_factor.assign(data.covariance_factor().begin(), data.covariance_factor().end());
    for (double& entry : problem.hessian_factor)
    {
        entry *= std::sqrt(2.0);
    }
}

Homogeneous homogeneous_problem(const Spec& spec, const Data& data, std::span<const double> v)
{
    Homogeneous homogeneous;
    qp::Problem& problem = homogeneous.problem;
    problem.n = data.asset_count();
    add_curvature(data, problem);
    problem.linear.assign(problem.n, 0.0);
    problem.equality_rows.assign(v.begin(), v.end());
    problem.equality_rhs.push_back(1.0);
    add_bound_rows(spec, homogeneous);
    homogeneous.kappa_row = add_row(problem, std::vector<double>(problem.n, 1.0));
    if (spec.leverage.has_value())
    {
        problem.gross_exposure = qp::GrossExposure{std::vector<double>(problem.n, *spec.leverage), 0.0};
    }
    return homogeneous;
}

Result without_point(const Data& data, qp::Status status, int iterations)
{
    qp::Solution solution;
    solution.status = status;
    solution.iterations = iterations;
    return to_result(data, ConstraintRows{}, solution);
}

// The supremum is not attained: no point, so NaN weights and multipliers.
Result unbounded(const Data& data, int iterations)
{
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    Result result;
    result.status = SolveStatus::Unbounded;
    result.weights.assign(data.weight_count(), kNaN);
    result.bound_multipliers.assign(data.weight_count(), kNaN);
    result.residual = kNaN;
    result.iterations = iterations;
    return result;
}

double sum_of(std::span<const double> values)
{
    double total = 0.0;
    for (double value : values)
    {
        total += value;
    }
    return total;
}

// sqrt(y'Sigma y) = |L'y|.
double volatility_of(const Data& data, std::span<const double> y)
{
    std::size_t n = data.asset_count();
    std::span<const double> factor = data.covariance_factor();
    double variance = 0.0;
    for (std::size_t j = 0; j < n; ++j)
    {
        double entry = 0.0;
        for (std::size_t i = j; i < n; ++i)
        {
            entry += factor[i * n + j] * y[i];
        }
        variance += entry * entry;
    }
    return std::sqrt(variance);
}

// |Sigma| |y|, componentwise.
std::vector<double> magnitude_product(const Data& data, std::span<const double> y)
{
    std::size_t n = data.asset_count();
    std::span<const double> covariance = data.covariance();
    std::vector<double> product(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            product[i] += std::fabs(covariance[i * n + j] * y[j]);
        }
    }
    return product;
}

double magnitude_dot(std::span<const double> a, std::span<const double> b)
{
    double total = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        total += std::fabs(a[i]) * std::fabs(b[i]);
    }
    return total;
}

// The rounding of kappa = 1'y. Where only v'y = 1 is active, Sigma y = sigma_y^2 v and kappa = sigma_y^2 z'v with
// z = Sigma^-1 1 (two triangular solves on the factor of Data). A backward error |dSigma| <= n eps |Sigma| of the
// solve and the rounding |dv| <= eps m of v (m = |mu| + |rf| or sigma, the magnitudes v is computed from) move it
// to first order by |z|'|dSigma||y| + sigma_y^2 |z|'|dv|, so kappa is known to
// n eps (|z|'|Sigma||y| + sigma_y^2 |z|'m).
double kappa_noise(const Data& data, std::span<const double> y, std::span<const double> magnitude)
{
    std::size_t n = data.asset_count();
    std::vector<double> z(n, 1.0);
    linalg::solve_lower(data.covariance_factor(), n, z);
    linalg::solve_lower_transposed(data.covariance_factor(), n, z);
    double sigma = volatility_of(data, y);
    double curvature = magnitude_dot(z, magnitude_product(data, y));
    double rounding = sigma * sigma * magnitude_dot(z, magnitude);
    return static_cast<double>(n) * kEpsilon * (curvature + rounding);
}

std::vector<double> user_weights(const Data& data, std::span<const double> y, double kappa)
{
    std::vector<double> weights(y.begin(), y.end());
    for (double& weight : weights)
    {
        weight /= kappa;
    }
    if (data.risk_free().has_value())
    {
        weights.push_back(0.0);
    }
    return weights;
}

// With c = kappa / (2 sigma_y^3), 2 Sigma y = 2 sigma_y^2 v + C'z + eta 1 + sum_j zeta_j (L 1 - s_j) turns
// grad f(w) = -kappa v / sigma_y + kappa Sigma y / sigma_y^3 into c (C'z + eta 1 + sum_j zeta_j (L 1 - s_j))
// (docs/guide/math.md): beta_i += c z_i at a lower row, beta_j -= c z_j at an upper row,
// lambda_lev = -c sum zeta, and the budget term
// c (eta - sum_lower z_i lo_i + sum_upper z_j hi_j + L sum zeta), which with cash is beta_0 = -that.
Result at_point(const Spec& spec, const Data& data, const Homogeneous& homogeneous,
                const qp::Solution& solution)
{
    double kappa = sum_of(solution.x);
    double sigma = volatility_of(data, solution.x);
    double scale = kappa / (2.0 * sigma * sigma * sigma);
    std::span<const double> z = solution.inequality_multipliers;
    Result result;
    result.status = SolveStatus::Optimal;
    result.weights = user_weights(data, solution.x, kappa);
    result.bound_multipliers.assign(data.weight_count(), 0.0);
    double budget = z[homogeneous.kappa_row];
    for (const BoundRow& bound : homogeneous.bounds)
    {
        if (bound.lower)
        {
            result.bound_multipliers[bound.weight] += scale * z[bound.row];
            budget -= z[bound.row] * lower_of(spec, bound.weight);
            continue;
        }
        result.bound_multipliers[bound.weight] -= scale * z[bound.row];
        budget += z[bound.row] * upper_of(spec, bound.weight);
    }
    if (spec.leverage.has_value())
    {
        double zeta = sum_of(solution.gross_exposure_multipliers);
        budget += *spec.leverage * zeta;
        result.leverage_multiplier = -scale * zeta;
    }
    result.budget_multiplier = scale * budget;
    if (data.risk_free().has_value())
    {
        result.bound_multipliers.back() = -result.budget_multiplier;
        result.budget_multiplier = std::numeric_limits<double>::quiet_NaN();
    }
    result.iterations = solution.iterations;
    return result;
}

// Whether kappa = 1'y is bounded away from 0 on the homogeneous cone. At kappa = 0 the rows read y_i >= 0 (finite
// lower), y_i <= 0 (finite upper), 1'y = 0 and, with a leverage, ||y||_1 <= 0, and v'y = 1 needs y != 0 there. So
// the cone is bounded exactly when a leverage is set or no pair i != j has upper_i = +inf and lower_j = -inf (a
// direction with 1'd = 0 needs some d_i > 0 and some d_j < 0, j != i); a Fixed box is finite on both sides. Then
// UNBOUNDED is impossible.
bool cone_bounded(const Spec& spec, std::size_t n)
{
    if (spec.leverage.has_value())
    {
        return true;
    }
    std::vector<std::size_t> rising;
    std::vector<std::size_t> falling;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!std::isfinite(upper_of(spec, i)))
        {
            rising.push_back(i);
        }
        if (!std::isfinite(lower_of(spec, i)))
        {
            falling.push_back(i);
        }
    }
    if (rising.empty() || falling.empty())
    {
        return true;
    }
    return rising.size() == 1 && falling.size() == 1 && rising.front() == falling.front();
}

// The closure point kappa = 0: the row 1'y >= 0 active, or kappa <= 0. On an unbounded cone also a kappa within
// kKappaNoiseFactor times its rounding (kappa_noise) with the row inactive, held by the solver within its
// feasibility tolerance: its sign, and with it w = y / kappa, is not decided by the data.
bool at_closure(const Spec& spec, const Data& data, const qp::Solution& solution, bool kappa_active,
                std::span<const double> magnitude)
{
    double kappa = sum_of(solution.x);
    if (kappa_active || !(kappa > 0.0))
    {
        return true;
    }
    if (cone_bounded(spec, data.asset_count()))
    {
        return false;
    }
    return !(kappa > kKappaNoiseFactor * kappa_noise(data, solution.x, magnitude));
}

bool cash_can_be_zero(const Spec& spec, const Data& data)
{
    std::size_t n = data.asset_count();
    return !data.risk_free().has_value() || (lower_of(spec, n) <= 0.0 && upper_of(spec, n) >= 0.0);
}

// The supremum at kappa = 0 (at_closure) is UNBOUNDED.
Result solve_homogeneous(const Spec& spec, const Data& data, std::span<const double> v,
                         std::span<const double> magnitude)
{
    if (!cash_can_be_zero(spec, data))
    {
        return without_point(data, qp::Status::Infeasible, 0);
    }
    Homogeneous homogeneous = homogeneous_problem(spec, data, v);
    qp::Solution solution = qp::solve(homogeneous.problem);
    if (solution.status != qp::Status::Optimal)
    {
        return without_point(data, solution.status, solution.iterations);
    }
    bool kappa_active = std::ranges::binary_search(solution.active_inequalities, homogeneous.kappa_row);
    if (at_closure(spec, data, solution, kappa_active, magnitude))
    {
        return unbounded(data, solution.iterations);
    }
    return at_point(spec, data, homogeneous, solution);
}

std::vector<double> excess_return(const Data& data)
{
    double rate = data.risk_free().value_or(0.0);
    std::vector<double> excess(data.mu().begin(), data.mu().end());
    for (double& value : excess)
    {
        value -= rate;
    }
    return excess;
}

// |mu| + |rf|, the magnitudes whose rounding mu - rf 1 carries.
std::vector<double> excess_magnitude(const Data& data)
{
    double rate = std::fabs(data.risk_free().value_or(0.0));
    std::vector<double> magnitude(data.mu().begin(), data.mu().end());
    for (double& value : magnitude)
    {
        value = std::fabs(value) + rate;
    }
    return magnitude;
}

std::vector<double> asset_volatilities(const Data& data)
{
    std::size_t n = data.asset_count();
    std::vector<double> sigma(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        sigma[i] = std::sqrt(data.covariance()[i * n + i]);
    }
    return sigma;
}

} // namespace

Result solve_max_sharpe(const Spec& spec, const Data& data)
{
    return solve_homogeneous(spec, data, excess_return(data), excess_magnitude(data));
}

Result solve_max_diversification(const Spec& spec, const Data& data)
{
    std::vector<double> sigma = asset_volatilities(data);
    return solve_homogeneous(spec, data, sigma, sigma);
}

} // namespace po::problems
