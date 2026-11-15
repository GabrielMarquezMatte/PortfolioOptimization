#pragma once

#include <span>

#include "core/data.hpp"

namespace po::metrics
{

// The weights have data.weight_count() entries (cash last); the return needs mu and the volatility the covariance.
double portfolio_return(const Data& data, std::span<const double> weights);

// Over the risky weights only: cash contributes no variance.
double portfolio_volatility(const Data& data, std::span<const double> weights);

} // namespace po::metrics
