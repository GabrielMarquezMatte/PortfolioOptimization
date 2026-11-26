#include "estimation/moments.hpp"

#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <utility>

namespace po::estimation
{
namespace
{

Failure invalid(std::string message)
{
    return Failure{FailureKind::InvalidArgument, std::move(message)};
}

std::optional<Failure> check_sizes(std::span<const double> returns, std::size_t t, std::size_t n)
{
    if (n == 0)
    {
        return invalid("asset count must be positive");
    }
    if (t < 2)
    {
        return invalid("returns need at least two observations");
    }
    if (t > std::numeric_limits<std::size_t>::max() / n)
    {
        return invalid("returns size is too large");
    }
    if (returns.size() != t * n)
    {
        return invalid("returns has " + std::to_string(returns.size()) + " entries; expected " +
                       std::to_string(t * n));
    }
    return std::nullopt;
}

std::optional<Failure> check_finite(std::span<const double> returns, std::size_t n)
{
    for (std::size_t k = 0; k < returns.size(); ++k)
    {
        if (!std::isfinite(returns[k]))
        {
            return invalid("returns is not finite at (" + std::to_string(k / n) + ", " + std::to_string(k % n) + ")");
        }
    }
    return std::nullopt;
}

std::vector<double> column_means(std::span<const double> returns, std::size_t t, std::size_t n)
{
    std::vector<double> mu(n, 0.0);
    for (std::size_t row = 0; row < t; ++row)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            mu[i] += returns[row * n + i];
        }
    }
    for (double& value : mu)
    {
        value /= static_cast<double>(t);
    }
    return mu;
}

double covariance_entry(std::span<const double> returns, std::span<const double> mu, std::size_t t, std::size_t i,
                        std::size_t j)
{
    std::size_t n = mu.size();
    double sum = 0.0;
    for (std::size_t row = 0; row < t; ++row)
    {
        sum += (returns[row * n + i] - mu[i]) * (returns[row * n + j] - mu[j]);
    }
    return sum / static_cast<double>(t - 1);
}

} // namespace

std::expected<Moments, Failure> sample_moments(std::span<const double> returns, std::size_t t, std::size_t n)
{
    std::optional<Failure> failure =
        check_sizes(returns, t, n).or_else([&] { return check_finite(returns, n); });
    if (failure.has_value())
    {
        return std::unexpected(std::move(*failure));
    }
    Moments moments{column_means(returns, t, n), std::vector<double>(n * n)};
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = i; j < n; ++j)
        {
            double entry = covariance_entry(returns, moments.mu, t, i, j);
            moments.covariance[i * n + j] = entry;
            moments.covariance[j * n + i] = entry;
        }
    }
    return moments;
}

} // namespace po::estimation
