#pragma once

#include <cstddef>
#include <span>
#include <vector>

#include "core/data.hpp"

namespace po::riskparity
{

// Sigma x, each entry accurate to its rounding even when its terms cancel.
std::vector<double> covariance_times(const Data& data, std::span<const double> x);

// r(w) = c(w) - b and f = |r|^2, with c_i = w_i (Sigma w)_i / (w'Sigma w) and b the budget. Holds f, its
// gradient 2 J'r, the Jacobian J of c (n x n, row-major) and, when asked for, the exact Hessian of f,
// 2 J'J + 2 sum_i r_i hess(c_i) (n x n, row-major; empty otherwise). w'Sigma w must be positive.
struct Linearization
{
    double value;
    std::vector<double> gradient;
    std::vector<double> jacobian;
    std::vector<double> hessian;
};

Linearization linearize(const Data& data, std::span<const double> budget, std::span<const double> w,
                        bool with_hessian);

// 2 J'J, the Gauss-Newton curvature of f (n x n, row-major).
std::vector<double> gauss_newton(std::span<const double> jacobian, std::size_t n);

// f(w + step) - f(w) and (w + step)'Sigma (w + step) - w'Sigma w, computed from the step itself so that each
// keeps its relative accuracy however small the step is; a difference of two values of f would not.
struct Change
{
    double value;
    double variance;
};

Change value_change(const Data& data, std::span<const double> budget, std::span<const double> w,
                    std::span<const double> step);

} // namespace po::riskparity
