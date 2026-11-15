#include "problems/tracking_error.hpp"

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "linalg/qr.hpp"
#include "qp/goldfarb_idnani.hpp"

namespace po::problems
{
namespace
{

std::vector<double> excess_returns(const Data& data)
{
    std::vector<double> excess(data.returns().begin(), data.returns().end());
    double rate = data.risk_free().value_or(0.0);
    for (double& value : excess)
    {
        value -= rate;
    }
    return excess;
}

std::vector<double> normal_matrix(const std::vector<double>& excess, std::size_t periods, std::size_t n)
{
    std::vector<double> normal(n * n, 0.0);
    for (std::size_t t = 0; t < periods; ++t)
    {
        const double* row = excess.data() + t * n;
        for (std::size_t i = 0; i < n; ++i)
        {
            for (std::size_t j = i; j < n; ++j)
            {
                normal[i * n + j] += row[i] * row[j];
            }
        }
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = i + 1; j < n; ++j)
        {
            normal[j * n + i] = normal[i * n + j];
        }
    }
    return normal;
}

std::vector<double> linear_term(const Data& data, const std::vector<double>& excess)
{
    std::size_t n = data.asset_count();
    double rate = data.risk_free().value_or(0.0);
    std::vector<double> linear(n, 0.0);
    for (std::size_t t = 0; t < data.observation_count(); ++t)
    {
        double target = data.benchmark()[t] - rate;
        for (std::size_t i = 0; i < n; ++i)
        {
            linear[i] -= 2.0 * excess[t * n + i] * target;
        }
    }
    return linear;
}

Failure rank_failure()
{
    return Failure{FailureKind::InvalidArgument,
                   "returns minus the risk-free rate do not have full column rank numerically"};
}

// G = 2 R~'R~ with the factor sqrt(2) R_qr', since R~'R~ = R_qr'R_qr.
Translation assemble(const Spec& spec, const Data& data, std::vector<double> normal,
                     std::vector<double> factor, const std::vector<double>& excess)
{
    Translation translation;
    qp::Problem& problem = translation.problem;
    problem.n = data.asset_count();
    problem.hessian = std::move(normal);
    for (double& entry : problem.hessian)
    {
        entry *= 2.0;
    }
    problem.hessian_factor = std::move(factor);
    double factor_scale = std::sqrt(2.0);
    for (double& entry : problem.hessian_factor)
    {
        entry *= factor_scale;
    }
    problem.linear = linear_term(data, excess);
    translation.rows = add_constraint_rows(spec, data, problem);
    return translation;
}

} // namespace

std::expected<Translation, Failure> tracking_error_problem(const Spec& spec, const Data& data)
{
    std::size_t periods = data.observation_count();
    std::size_t n = data.asset_count();
    std::vector<double> excess = excess_returns(data);
    if (periods < n)
    {
        return std::unexpected(rank_failure());
    }
    linalg::Qr qr = linalg::householder_qr(excess, periods, n);
    if (!linalg::has_full_column_rank(qr))
    {
        return std::unexpected(rank_failure());
    }
    return assemble(spec, data, normal_matrix(excess, periods, n), linalg::lower_factor(qr), excess);
}

std::expected<Result, Failure> solve_tracking_error(const Spec& spec, const Data& data)
{
    std::expected<Translation, Failure> translation = tracking_error_problem(spec, data);
    if (!translation.has_value())
    {
        return std::unexpected(std::move(translation.error()));
    }
    return to_result(data, translation->rows, qp::solve(translation->problem));
}

} // namespace po::problems
