#pragma once

#include <optional>

#include "core/data.hpp"
#include "core/failure.hpp"

namespace po::estimation
{

// Refuses a data object that has no covariance, no mu when needs_mu, or a risk-free rate: the estimation
// utilities work on raw returns, never on excess returns.
std::optional<Failure> check_inputs(const Data& data, bool needs_mu);

} // namespace po::estimation
