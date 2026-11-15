#pragma once

#include "qp/problem.hpp"

namespace po::qp
{

// Goldfarb and Idnani (1983) dual active-set method. Throws std::invalid_argument when the Hessian is not
// positive definite or a matrix does not match n. When the status is not Optimal, x and the multipliers
// are NaN.
Solution solve(const Problem& problem);

} // namespace po::qp
