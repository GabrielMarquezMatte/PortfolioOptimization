#pragma once

#include <optional>

#include "core/data.hpp"
#include "core/failure.hpp"
#include "estimation/views.hpp"

namespace po::estimation
{

// Refuses a data object that has no covariance, no mu when needs_mu, or a risk-free rate: the estimation
// utilities work on raw returns, never on excess returns.
std::optional<Failure> check_inputs(const Data& data, bool needs_mu);

// Refuses a value that is not finite and positive, with the given message.
std::optional<Failure> check_positive(double value, const char* message);

// Refuses views whose asset count differs from the data's.
std::optional<Failure> check_views(const Data& data, const Views& views);

} // namespace po::estimation
