#include "problems/min_variance.hpp"

#include "qp/goldfarb_idnani.hpp"

namespace po::problems
{

// G = 2 Sigma, so that 0.5 w'Gw is the objective w'Sigma w exactly and the QP multipliers are the user's.
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
