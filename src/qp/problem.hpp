#pragma once

#include <cstddef>
#include <vector>

namespace po::qp
{

// min 0.5 x'Gx + g'x  subject to  E x = e (rows of E),  C x >= c (rows of C).
struct Problem
{
    std::size_t n = 0;
    std::vector<double> hessian;
    std::vector<double> linear;
    std::vector<double> equality_rows;
    std::vector<double> equality_rhs;
    std::vector<double> inequality_rows;
    std::vector<double> inequality_rhs;
};

enum class Status
{
    Optimal,
    Infeasible,
    MaxIterations
};

// Stationarity: G x + g = E' y + C' z, with z >= 0 and z_i (C_i x - c_i) = 0.
struct Solution
{
    Status status = Status::Infeasible;
    std::vector<double> x;
    std::vector<double> equality_multipliers;
    std::vector<double> inequality_multipliers;
    int iterations = 0;
};

} // namespace po::qp
