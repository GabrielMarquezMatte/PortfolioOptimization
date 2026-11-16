#include "problems/constraints.hpp"

#include <cmath>
#include <limits>
#include <span>
#include <utility>

namespace po::problems
{
namespace
{

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

std::size_t add_equality(qp::Problem& problem, std::span<const double> row, double rhs)
{
    problem.equality_rows.insert(problem.equality_rows.end(), row.begin(), row.end());
    problem.equality_rhs.push_back(rhs);
    return problem.equality_rhs.size() - 1;
}

std::size_t add_inequality(qp::Problem& problem, std::span<const double> row, double rhs)
{
    problem.inequality_rows.insert(problem.inequality_rows.end(), row.begin(), row.end());
    problem.inequality_rhs.push_back(rhs);
    return problem.inequality_rhs.size() - 1;
}

std::vector<double> unit_row(std::size_t n, std::size_t index, double value)
{
    std::vector<double> row(n, 0.0);
    row[index] = value;
    return row;
}

void add_equalities(const Spec& spec, const Data& data, qp::Problem& problem, ConstraintRows& rows)
{
    std::size_t n = data.asset_count();
    if (!data.risk_free().has_value())
    {
        rows.budget = add_equality(problem, std::vector<double>(n, 1.0), 1.0);
    }
    if (!spec.target_return.has_value())
    {
        return;
    }
    double rate = data.risk_free().value_or(0.0);
    std::vector<double> excess(data.mu().begin(), data.mu().end());
    for (double& value : excess)
    {
        value -= rate;
    }
    rows.target_return = add_equality(problem, excess, *spec.target_return - rate);
}

void add_asset_bounds(const Spec& spec, std::size_t n, qp::Problem& problem, ConstraintRows& rows)
{
    for (std::size_t i = 0; i < n; ++i)
    {
        double lower = lower_of(spec, i);
        double upper = upper_of(spec, i);
        if (std::isfinite(lower))
        {
            rows.bounds.push_back(BoundRow{i, add_inequality(problem, unit_row(n, i, 1.0), lower), true});
        }
        if (std::isfinite(upper))
        {
            rows.bounds.push_back(BoundRow{i, add_inequality(problem, unit_row(n, i, -1.0), -upper), false});
        }
    }
}

void add_cash_bounds(const Spec& spec, std::size_t n, qp::Problem& problem, ConstraintRows& rows)
{
    double lower = lower_of(spec, n);
    double upper = upper_of(spec, n);
    if (std::isfinite(lower))
    {
        rows.bounds.push_back(BoundRow{n, add_inequality(problem, std::vector<double>(n, -1.0), lower - 1.0), true});
    }
    if (std::isfinite(upper))
    {
        rows.bounds.push_back(BoundRow{n, add_inequality(problem, std::vector<double>(n, 1.0), 1.0 - upper), false});
    }
}

SolveStatus status_of(qp::Status status)
{
    switch (status)
    {
    case qp::Status::Optimal:
        return SolveStatus::Optimal;
    case qp::Status::Infeasible:
        return SolveStatus::Infeasible;
    case qp::Status::MaxIterations:
        return SolveStatus::MaxIterations;
    }
    std::unreachable();
}

double multiplier_at(std::optional<std::size_t> row, std::span<const double> multipliers, double scale)
{
    if (!row.has_value())
    {
        return kNaN;
    }
    return multipliers[*row] / scale;
}

std::vector<double> user_weights(const Data& data, std::span<const double> x)
{
    std::vector<double> weights(x.begin(), x.end());
    if (!data.risk_free().has_value())
    {
        return weights;
    }
    double cash = 1.0;
    for (double weight : x)
    {
        cash -= weight;
    }
    weights.push_back(cash);
    return weights;
}

std::vector<double> bound_multipliers(const Data& data, const ConstraintRows& rows, std::span<const double> z,
                                      double scale)
{
    std::vector<double> multipliers(data.weight_count(), 0.0);
    for (const BoundRow& bound : rows.bounds)
    {
        if (bound.lower)
        {
            multipliers[bound.weight] += z[bound.row] / scale;
            continue;
        }
        multipliers[bound.weight] -= z[bound.row] / scale;
    }
    return multipliers;
}

Result without_point(const Data& data, const qp::Solution& solution)
{
    Result result;
    result.status = status_of(solution.status);
    result.weights.assign(data.weight_count(), kNaN);
    result.bound_multipliers.assign(data.weight_count(), kNaN);
    result.residual = kNaN;
    result.iterations = solution.iterations;
    return result;
}

// The family row (-s_j)'x >= -L enters the QP stationarity as -sum_j zeta_j s_j = lambda_lev g, with
// g = sum_j zeta_j s_j / sum_j zeta_j in the subdifferential of ||w||_1: lambda_lev = -sum_j zeta_j, 0 when no
// generated row is active.
double leverage_multiplier(const ConstraintRows& rows, const qp::Solution& solution, double scale)
{
    if (!rows.gross_exposure)
    {
        return kNaN;
    }
    double multiplier = 0.0;
    for (double zeta : solution.gross_exposure_multipliers)
    {
        multiplier -= zeta;
    }
    return multiplier / scale;
}

Result at_point(const Data& data, const ConstraintRows& rows, const qp::Solution& solution, double scale)
{
    Result result;
    result.status = SolveStatus::Optimal;
    result.weights = user_weights(data, solution.x);
    result.budget_multiplier = multiplier_at(rows.budget, solution.equality_multipliers, scale);
    result.target_return_multiplier = multiplier_at(rows.target_return, solution.equality_multipliers, scale);
    result.bound_multipliers = bound_multipliers(data, rows, solution.inequality_multipliers, scale);
    result.leverage_multiplier = leverage_multiplier(rows, solution, scale);
    result.iterations = solution.iterations;
    return result;
}

} // namespace

double lower_of(const Spec& spec, std::size_t weight)
{
    if (spec.asset_lower.empty())
    {
        return spec.lower;
    }
    return spec.asset_lower[weight];
}

double upper_of(const Spec& spec, std::size_t weight)
{
    if (spec.asset_upper.empty())
    {
        return spec.upper;
    }
    return spec.asset_upper[weight];
}

ConstraintRows add_constraint_rows(const Spec& spec, const Data& data, qp::Problem& problem)
{
    ConstraintRows rows;
    add_equalities(spec, data, problem, rows);
    add_asset_bounds(spec, data.asset_count(), problem, rows);
    if (data.risk_free().has_value())
    {
        add_cash_bounds(spec, data.asset_count(), problem, rows);
    }
    if (spec.leverage.has_value())
    {
        // sum_i |w_i| <= L over the risky weights, the QP variables: 0'x - ||x||_1 >= -L.
        std::vector<double> offset(data.asset_count(), 0.0);
        problem.gross_exposure = qp::GrossExposure{std::move(offset), -*spec.leverage};
        rows.gross_exposure = true;
    }
    return rows;
}

Result to_result(const Data& data, const ConstraintRows& rows, const qp::Solution& solution,
                 double multiplier_scale)
{
    if (solution.status != qp::Status::Optimal)
    {
        return without_point(data, solution);
    }
    return at_point(data, rows, solution, multiplier_scale);
}

} // namespace po::problems
