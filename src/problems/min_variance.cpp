#include "problems/min_variance.hpp"

#include <numbers>

#include "qp/goldfarb_idnani.hpp"

namespace po::problems
{

// G = 2 Sigma, so that 0.5 w'Gw is the objective w'Sigma w exactly and the QP multipliers are the user's.
// Its factor is sqrt(2) L from the factor L of Sigma, so that positive definiteness stays the one decision of
// Data::create.
Translation min_variance_problem(const Spec& spec, const Data& data)
{
    Translation translation;
    qp::Problem& problem = translation.problem;
    problem.n = data.asset_count();
    problem.hessian.assign(data.covariance().begin(), data.covariance().end());
    for (double& entry : problem.hessian)
    {
        entry *= 2.0;
    }
    problem.hessian_factor.assign(data.covariance_factor().begin(), data.covariance_factor().end());
    for (double& entry : problem.hessian_factor)
    {
        entry *= std::numbers::sqrt2;
    }
    problem.linear.assign(problem.n, 0.0);
    translation.rows = add_constraint_rows(spec, data, problem);
    return translation;
}

Result solve_min_variance(const Spec& spec, const Data& data)
{
    Translation translation = min_variance_problem(spec, data);
    return to_result(data, translation.rows, qp::solve(translation.problem));
}

} // namespace po::problems
