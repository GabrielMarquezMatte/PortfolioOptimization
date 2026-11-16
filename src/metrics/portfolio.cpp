#include "metrics/portfolio.hpp"

#include <cmath>
#include <cstddef>
#include <limits>

namespace po::metrics
{
namespace
{

double risky_return(const Data& data, std::span<const double> weights)
{
    double total = 0.0;
    std::span<const double> mu = data.mu();
    for (std::size_t i = 0; i < mu.size(); ++i)
    {
        total += weights[i] * mu[i];
    }
    return total;
}

double row_dot(std::span<const double> row, std::span<const double> weights)
{
    double total = 0.0;
    for (std::size_t j = 0; j < row.size(); ++j)
    {
        total += row[j] * weights[j];
    }
    return total;
}

double variance(const Data& data, std::span<const double> weights)
{
    std::size_t n = data.asset_count();
    double total = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        total += weights[i] * row_dot(data.covariance().subspan(i * n, n), weights);
    }
    return total;
}

double cash_return(const Data& data, std::span<const double> weights)
{
    if (!data.risk_free().has_value())
    {
        return 0.0;
    }
    return weights[data.asset_count()] * *data.risk_free();
}

double weighted_volatilities(const Data& data, std::span<const double> weights)
{
    std::size_t n = data.asset_count();
    double total = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        total += weights[i] * std::sqrt(data.covariance()[i * n + i]);
    }
    return total;
}

double squared_tracking_gaps(const Data& data, std::span<const double> weights)
{
    std::size_t n = data.asset_count();
    double cash = cash_return(data, weights);
    double total = 0.0;
    for (std::size_t t = 0; t < data.observation_count(); ++t)
    {
        double gap = row_dot(data.returns().subspan(t * n, n), weights) + cash - data.benchmark()[t];
        total += gap * gap;
    }
    return total;
}

} // namespace

double portfolio_return(const Data& data, std::span<const double> weights)
{
    return risky_return(data, weights) + cash_return(data, weights);
}

double portfolio_volatility(const Data& data, std::span<const double> weights)
{
    return std::sqrt(variance(data, weights));
}

double utility(const Data& data, std::span<const double> weights, double delta)
{
    return portfolio_return(data, weights) - 0.5 * delta * variance(data, weights);
}

double sharpe_ratio(const Data& data, std::span<const double> weights)
{
    double volatility = portfolio_volatility(data, weights);
    if (volatility == 0.0)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return (portfolio_return(data, weights) - data.risk_free().value_or(0.0)) / volatility;
}

double diversification_ratio(const Data& data, std::span<const double> weights)
{
    double volatility = portfolio_volatility(data, weights);
    if (volatility == 0.0)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return weighted_volatilities(data, weights) / volatility;
}

double tracking_error(const Data& data, std::span<const double> weights)
{
    return std::sqrt(squared_tracking_gaps(data, weights) / static_cast<double>(data.observation_count()));
}

} // namespace po::metrics
