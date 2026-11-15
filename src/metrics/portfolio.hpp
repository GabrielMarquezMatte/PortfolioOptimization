#pragma once

#include <span>

#include "core/data.hpp"

namespace po::metrics
{

// The weights have data.weight_count() entries (cash last); the return needs mu and the volatility the covariance.
double portfolio_return(const Data& data, std::span<const double> weights);

// Over the risky weights only: cash contributes no variance.
double portfolio_volatility(const Data& data, std::span<const double> weights);

// return - delta / 2 w' Sigma w; needs mu and the covariance.
double utility(const Data& data, std::span<const double> weights, double delta);

// sqrt(sum_t (r_t'w + w_0 rf - b_t)^2 / T), the cash term present only with cash; needs returns.
double tracking_error(const Data& data, std::span<const double> weights);

} // namespace po::metrics
