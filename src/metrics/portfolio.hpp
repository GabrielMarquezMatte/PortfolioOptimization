#pragma once

#include <span>
#include <vector>

#include "core/data.hpp"

namespace po::metrics
{

// The weights have data.weight_count() entries (cash last); the return needs mu and the volatility the
// covariance.
double portfolio_return(const Data& data, std::span<const double> weights);

// Over the risky weights only: cash contributes no variance.
double portfolio_volatility(const Data& data, std::span<const double> weights);

// return - delta / 2 w' Sigma w; needs mu and the covariance.
double utility(const Data& data, std::span<const double> weights, double delta);

// (return - rf) / volatility with rf = 0 without cash; needs mu and the covariance. NaN when the volatility
// is 0.
double sharpe_ratio(const Data& data, std::span<const double> weights);

// sigma'w / volatility over the risky weights, sigma = sqrt(diag Sigma); needs the covariance. NaN when the
// volatility is 0.
double diversification_ratio(const Data& data, std::span<const double> weights);

// sqrt(sum_t (r_t'w + w_0 rf - b_t)^2 / T), the cash term present only with cash; needs returns.
double tracking_error(const Data& data, std::span<const double> weights);

// c_i = w_i (Sigma w)_i / w'Sigma w over the risky weights; the cash entry (last, when present) is 0. Sums to
// 1; all NaN when the volatility is 0. Needs the covariance.
std::vector<double> risk_contributions(const Data& data, std::span<const double> weights);

// sum_i (c_i - b_i)^2 over the risky weights; NaN when the volatility is 0. budget has asset_count entries.
double risk_budget_residual(const Data& data, std::span<const double> weights,
                            std::span<const double> budget);

} // namespace po::metrics
