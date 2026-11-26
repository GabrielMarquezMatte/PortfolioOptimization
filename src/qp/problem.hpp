#pragma once

#include <cstddef>
#include <optional>
#include <vector>

namespace po::qp
{

// a'x - ||x||_1 >= b: the implicit family of the rows (a - s)'x >= b over every s in {-1, 1}^n.
struct GrossExposure
{
    std::vector<double> offset;
    double rhs = 0.0;
};

// min 0.5 x'Gx + g'x  subject to  E x = e (rows of E),  C x >= c (rows of C) and, when gross_exposure is set,
// its family. hessian_factor, when not empty, is the lower Cholesky factor L of G = L L' (n x n, row-major): the
// solver uses it as given, neither factors G nor decides again whether G is positive definite, and never reads
// hessian, which may then be empty (or G, n x n). Empty means the solver factors hessian, which must be G, n x n.
struct Problem
{
    std::size_t n = 0;
    std::vector<double> hessian;
    std::vector<double> hessian_factor;
    std::vector<double> linear;
    std::vector<double> equality_rows;
    std::vector<double> equality_rhs;
    std::vector<double> inequality_rows;
    std::vector<double> inequality_rhs;
    std::optional<GrossExposure> gross_exposure;
};

enum class Status
{
    Optimal,
    Infeasible,
    MaxIterations
};

// Stationarity: G x + g = E' y + C' z + sum_j zeta_j (a - s_j), with z >= 0 and z_i (C_i x - c_i) = 0, over the
// generated rows (a - s_j)'x >= b of the gross exposure family active at the solution, zeta_j >= 0.
struct Solution
{
    Status status = Status::Infeasible;
    std::vector<double> x;
    std::vector<double> equality_multipliers;
    std::vector<double> inequality_multipliers;
    // Inequality rows active at the solution, ascending; empty unless Optimal.
    std::vector<std::size_t> active_inequalities;
    // The active generated rows, one sign vector s_j of n entries -1 or +1 each (k x n, row-major), in
    // lexicographic order of s_j, and their multipliers zeta_j (k); empty unless Optimal. The active set is
    // identified by active_inequalities and these sign vectors.
    std::vector<double> gross_exposure_signs;
    std::vector<double> gross_exposure_multipliers;
    int iterations = 0;
};

} // namespace po::qp
