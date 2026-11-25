#pragma once

#include <cstddef>
#include <expected>
#include <span>
#include <vector>

#include "core/data.hpp"
#include "core/failure.hpp"
#include "estimation/views.hpp"

namespace po::estimation
{

// How far the market weights may be from summing to 1.
inline constexpr double kMarketWeightTolerance = 1e-8;

// Pi = delta Sigma w_mkt, the equilibrium return in excess of the risk-free rate. data.mu() is not used.
std::expected<std::vector<double>, Failure> equilibrium(const Data& data, std::span<const double> w_mkt,
                                                        double delta);

// He and Litterman: Omega = diag(tau p_k Sigma p_k'), k x k row-major. The Omega stored in views is not read.
std::expected<std::vector<double>, Failure> omega_he_litterman(const Data& data, const Views& views, double tau);

// Idzorek (Walters): Omega = diag(((1 - c_k) / c_k) tau p_k Sigma p_k') with 0 < c_k <= 1; c_k = 1 gives a zero
// entry, the full-confidence view.
std::expected<std::vector<double>, Failure> omega_idzorek(const Data& data, const Views& views, double tau,
                                                          std::span<const double> confidence);

} // namespace po::estimation
