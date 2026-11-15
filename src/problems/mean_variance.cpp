#include "problems/mean_variance.hpp"

#include <cmath>
#include <cstddef>
#include <utility>

#include "qp/goldfarb_idnani.hpp"

namespace po::problems
{
namespace
{

constexpr double kDefaultRiskAversion = 5.0;
constexpr double kDefaultImportance = 0.5;

void add_curvature(const Data& data, double quadratic, qp::Problem& problem)
{
    problem.hessian.assign(data.covariance().begin(), data.covariance().end());
    for (double& entry : problem.hessian)
    {
        entry *= quadratic;
    }
    double factor_scale = std::sqrt(quadratic);
    problem.hessian_factor.assign(data.covariance_factor().begin(), data.covariance_factor().end());
    for (double& entry : problem.hessian_factor)
    {
        entry *= factor_scale;
    }
}

// Without a linear term mu is not read: MIN_VARIANCE needs no mu.
void add_linear(const Data& data, double linear, qp::Problem& problem)
{
    problem.linear.assign(problem.n, 0.0);
    if (linear == 0.0)
    {
        return;
    }
    double rate = data.risk_free().value_or(0.0);
    for (std::size_t i = 0; i < problem.n; ++i)
    {
        problem.linear[i] = -linear * (data.mu()[i] - rate);
    }
}

} // namespace

MeanVariance mean_variance_of(const Spec& spec, const Data& data)
{
    switch (spec.objective)
    {
    case Objective::MinVariance:
        return MeanVariance{2.0, 0.0};
    case Objective::Utility:
        return MeanVariance{spec.risk_aversion.value_or(kDefaultRiskAversion), 1.0};
    case Objective::Importance:
    {
        double lambda = spec.importance.value_or(kDefaultImportance);
        return MeanVariance{2.0 * (1.0 - lambda), lambda};
    }
    case Objective::Kelly:
        return MeanVariance{1.0 / (1.0 + data.risk_free().value_or(0.0)), 1.0};
    default:
        std::unreachable();
    }
}

// G = quadratic Sigma and g = -linear (mu - rf 1), so that 0.5 w'Gw + g'w is the objective f up to a constant and
// the QP multipliers are the user's. The factor of G is sqrt(quadratic) L from the factor L of Sigma, so that
// positive definiteness stays the one decision of Data::create.
Translation mean_variance_problem(const Spec& spec, const Data& data, MeanVariance scales)
{
    Translation translation;
    qp::Problem& problem = translation.problem;
    problem.n = data.asset_count();
    add_curvature(data, scales.quadratic, problem);
    add_linear(data, scales.linear, problem);
    translation.rows = add_constraint_rows(spec, data, problem);
    return translation;
}

Result solve_mean_variance(const Spec& spec, const Data& data)
{
    Translation translation = mean_variance_problem(spec, data, mean_variance_of(spec, data));
    return to_result(data, translation.rows, qp::solve(translation.problem));
}

} // namespace po::problems
