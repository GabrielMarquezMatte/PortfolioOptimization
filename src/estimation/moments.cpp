#include "estimation/moments.hpp"

#include <algorithm>
#include <array>
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

constexpr std::size_t kTile = 4;

using Tile = std::array<double, kTile * kTile>;

// The returns less their means, with zero columns up to a multiple of kTile so that every tile is full.
struct CenteredPanel
{
    std::vector<double> values;
    std::size_t width;
};

CenteredPanel centered(std::span<const double> returns, std::span<const double> mu, std::size_t t)
{
    std::size_t n = mu.size();
    std::size_t width = (n + kTile - 1) / kTile * kTile;
    CenteredPanel panel{std::vector<double>(t * width, 0.0), width};
    for (std::size_t row = 0; row < t; ++row)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            panel.values[row * width + i] = returns[row * n + i] - mu[i];
        }
    }
    return panel;
}

// Every entry of the tile adds its products in row order, from 0.0, as the scalar sum of one entry would.
Tile tile_sums(const CenteredPanel& panel, std::size_t t, std::size_t first_row, std::size_t first_column)
{
    Tile sums{};
    for (std::size_t row = 0; row < t; ++row)
    {
        std::size_t base = row * panel.width;
        for (std::size_t r = 0; r < kTile; ++r)
        {
            for (std::size_t s = 0; s < kTile; ++s)
            {
                sums[r * kTile + s] += panel.values[base + first_row + r] * panel.values[base + first_column + s];
            }
        }
    }
    return sums;
}

void store_tile(std::span<double> covariance, std::size_t n, const Tile& tile, std::size_t first_row,
                std::size_t first_column)
{
    for (std::size_t r = 0; r < kTile && first_row + r < n; ++r)
    {
        std::size_t i = first_row + r;
        for (std::size_t j = std::max(i, first_column); j < first_column + kTile && j < n; ++j)
        {
            double entry = tile[r * kTile + (j - first_column)];
            covariance[i * n + j] = entry;
            covariance[j * n + i] = entry;
        }
    }
}

std::vector<double> sample_covariance(std::span<const double> returns, std::span<const double> mu, std::size_t t)
{
    std::size_t n = mu.size();
    CenteredPanel panel = centered(returns, mu, t);
    double divisor = static_cast<double>(t - 1);
    std::vector<double> covariance(n * n);
    for (std::size_t first_row = 0; first_row < n; first_row += kTile)
    {
        for (std::size_t first_column = first_row; first_column < n; first_column += kTile)
        {
            Tile tile = tile_sums(panel, t, first_row, first_column);
            for (double& entry : tile)
            {
                entry /= divisor;
            }
            store_tile(covariance, n, tile, first_row, first_column);
        }
    }
    return covariance;
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
    Moments moments{column_means(returns, t, n), {}};
    moments.covariance = sample_covariance(returns, moments.mu, t);
    return moments;
}

} // namespace po::estimation
