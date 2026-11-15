#include "metrics/portfolio.hpp"

#include <cmath>
#include <cstddef>

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

} // namespace

double portfolio_return(const Data& data, std::span<const double> weights)
{
    double total = risky_return(data, weights);
    if (data.risk_free().has_value())
    {
        total += weights[data.asset_count()] * *data.risk_free();
    }
    return total;
}

double portfolio_volatility(const Data& data, std::span<const double> weights)
{
    return std::sqrt(variance(data, weights));
}

} // namespace po::metrics
