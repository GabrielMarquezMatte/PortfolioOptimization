#pragma once

#include <span>
#include <vector>

#include "core/data.hpp"

namespace po::riskparity
{

struct SpinuResult
{
    std::vector<double> weights;
    int iterations;
    bool converged;
};

// Positive-orthant solution w = y / 1'y of min 1/2 y'Sigma y - b' log y, by damped Newton with Armijo,
// b = budget (positive, summing to 1, one entry per asset). The minimizer satisfies y_i (Sigma y)_i = b_i, so w
// has risk contributions b. Without cash and with Sigma from Data.
SpinuResult spinu(const Data& data, std::span<const double> budget);

} // namespace po::riskparity
