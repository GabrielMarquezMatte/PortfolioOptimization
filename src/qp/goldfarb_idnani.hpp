#pragma once

#include "qp/problem.hpp"

namespace po::qp
{

// The relative tolerance of every decision of the solver on a row: feasibility and dependence.
inline constexpr double kRelativeTolerance = 1e-12;

// A row r'x >= rhs holds at x when its slack is at least -feasibility_tolerance(rhs, ||r||, |x|max), with
// feasibility_tolerance = kRelativeTolerance max(1, |rhs|, ||r|| |x|max): relative to the point, never an
// absolute number. The solver decides by it which row is violated, and the generated row of the gross
// exposure family by the same rule.
double feasibility_tolerance(double rhs, double row_norm, double x_scale);

// Goldfarb and Idnani (1983) dual active-set method. Throws std::invalid_argument when a matrix does not
// match n or, without a hessian_factor, when the Hessian is not positive definite (linalg::cholesky). When
// the status is not Optimal, x and the multipliers are NaN.
Solution solve(const Problem& problem);

} // namespace po::qp
