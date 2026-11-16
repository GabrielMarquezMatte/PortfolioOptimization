#include "problems/homogeneous.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include "problems/constraints.hpp"
#include "qp/goldfarb_idnani.hpp"
#include "qp/problem.hpp"

namespace po::problems
{
namespace
{

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

Result unbounded(const Data& data, int iterations)
{
    Result result = without_point(data, qp::Status::Infeasible, iterations);
    result.status = SolveStatus::Unbounded;
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

bool cash_can_be_zero(const Spec& spec, const Data& data)
{
    std::size_t n = data.asset_count();
    return !data.risk_free().has_value() || (lower_of(spec, n) <= 0.0 && upper_of(spec, n) >= 0.0);
}

// The row 1'y >= 0 in the active set is the supremum at kappa = 0. A kappa that is not positive with the row
// inactive is the same closure point, held by the solver within its feasibility tolerance.
Result solve_homogeneous(const Spec& spec, const Data& data, std::span<const double> v)
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
    if (kappa_active || !(sum_of(solution.x) > 0.0))
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
    return solve_homogeneous(spec, data, excess_return(data));
}

Result solve_max_diversification(const Spec& spec, const Data& data)
{
    return solve_homogeneous(spec, data, asset_volatilities(data));
}

} // namespace po::problems
