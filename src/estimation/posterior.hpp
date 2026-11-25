#pragma once

#include <expected>
#include <vector>

#include "core/data.hpp"
#include "core/failure.hpp"
#include "estimation/views.hpp"

namespace po::estimation
{

struct Posterior
{
    std::vector<double> mu;
    // n x n row-major; empty when not asked for.
    std::vector<double> covariance;
};

// The Black-Litterman posterior with data.mu() as Pi: mu = Pi + B A^-1 (Q - P Pi) and
// Sigma_post = (1 + tau) Sigma - B A^-1 B', with B = tau Sigma P' and A = P B + Omega. The form needs no
// Omega^-1, so Omega = 0 (full confidence) is accepted. Needs the Omega of views when k >= 1.
std::expected<Posterior, Failure> posterior(const Data& data, const Views& views, double tau, bool with_covariance);

} // namespace po::estimation
