#pragma once

#include <vector>

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

// The basis J0 = L^{-T} that every solve starts from, kept between solves of problems with the same
// hessian_factor (the points of a parametric path): factor holds the bits of the L it was built from.
struct BasisCache
{
    std::vector<double> factor;
    std::vector<double> basis;
};

// solve(problem), starting from the basis of cache when problem.hessian_factor has the bits of cache.factor;
// the basis is built from the factor and kept in cache otherwise. Without a hessian_factor, solve(problem).
Solution solve(const Problem& problem, BasisCache& cache);

} // namespace po::qp
