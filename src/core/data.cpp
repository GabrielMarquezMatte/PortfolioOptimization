#include "core/data.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

#include "linalg/cholesky.hpp"

namespace po
{
namespace
{

Failure invalid(std::string message)
{
    return Failure{FailureKind::InvalidArgument, std::move(message)};
}

std::string position(std::size_t i, std::size_t j)
{
    return "(" + std::to_string(i) + ", " + std::to_string(j) + ")";
}

std::optional<Failure> check_sizes(std::size_t n, std::span<const double> mu, std::span<const double> covariance)
{
    if (n == 0)
    {
        return invalid("asset count must be positive");
    }
    if (!mu.empty() && mu.size() != n)
    {
        return invalid("mu has " + std::to_string(mu.size()) + " entries; expected " + std::to_string(n));
    }
    if (!covariance.empty() && covariance.size() != n * n)
    {
        return invalid("covariance has " + std::to_string(covariance.size()) + " entries; expected " +
                       std::to_string(n * n));
    }
    return std::nullopt;
}

std::optional<Failure> check_mu(std::span<const double> mu)
{
    for (std::size_t i = 0; i < mu.size(); ++i)
    {
        if (!std::isfinite(mu[i]))
        {
            return invalid("mu is not finite at " + std::to_string(i));
        }
    }
    return std::nullopt;
}

std::optional<Failure> check_finite(std::span<const double> covariance, std::size_t n)
{
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            if (!std::isfinite(covariance[i * n + j]))
            {
                return invalid("covariance is not finite at " + position(i, j));
            }
        }
    }
    return std::nullopt;
}

// Runs after check_finite: both magnitudes are finite and carry no sign, so std::max picks the value std::fmax would.
bool symmetric_pair(std::span<const double> covariance, std::size_t n, std::size_t i, std::size_t j)
{
    double scale = std::max(std::fabs(covariance[i * n + i]), std::fabs(covariance[j * n + j]));
    double gap = std::fabs(covariance[i * n + j] - covariance[j * n + i]);
    return gap <= kSymmetryTolerance * scale;
}

std::optional<Failure> check_symmetric(std::span<const double> covariance, std::size_t n)
{
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = i + 1; j < n; ++j)
        {
            if (!symmetric_pair(covariance, n, i, j))
            {
                return invalid("covariance is not symmetric at " + position(i, j));
            }
        }
    }
    return std::nullopt;
}

std::optional<Failure> check_covariance(std::span<const double> covariance, std::size_t n)
{
    if (covariance.empty())
    {
        return std::nullopt;
    }
    return check_finite(covariance, n).or_else([&] { return check_symmetric(covariance, n); });
}

std::optional<Failure> first_failure(std::size_t n, std::span<const double> mu, std::span<const double> covariance)
{
    return check_sizes(n, mu, covariance)
        .or_else([&] { return check_mu(mu); })
        .or_else([&] { return check_covariance(covariance, n); });
}

std::vector<double> symmetrized(std::span<const double> covariance, std::size_t n)
{
    std::vector<double> symmetric(covariance.size());
    if (covariance.empty())
    {
        return symmetric;
    }
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = i; j < n; ++j)
        {
            double entry = 0.5 * (covariance[i * n + j] + covariance[j * n + i]);
            symmetric[i * n + j] = entry;
            symmetric[j * n + i] = entry;
        }
    }
    return symmetric;
}

std::optional<std::size_t> first_non_finite(std::span<const double> values)
{
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (!std::isfinite(values[i]))
        {
            return i;
        }
    }
    return std::nullopt;
}

std::optional<Failure> check_return_sizes(std::size_t n, std::span<const double> returns,
                                          std::span<const double> benchmark, std::size_t t)
{
    if (t == 0)
    {
        return invalid("returns need at least one observation");
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
    if (benchmark.size() != t)
    {
        return invalid("benchmark has " + std::to_string(benchmark.size()) + " entries; expected " +
                       std::to_string(t));
    }
    return std::nullopt;
}

std::optional<Failure> check_return_values(std::size_t n, std::span<const double> returns,
                                           std::span<const double> benchmark)
{
    std::optional<std::size_t> bad_return = first_non_finite(returns);
    if (bad_return.has_value())
    {
        return invalid("returns is not finite at " + position(*bad_return / n, *bad_return % n));
    }
    std::optional<std::size_t> bad_benchmark = first_non_finite(benchmark);
    if (bad_benchmark.has_value())
    {
        return invalid("benchmark is not finite at " + std::to_string(*bad_benchmark));
    }
    return std::nullopt;
}

} // namespace

std::expected<Data, Failure> Data::create(std::size_t n, std::span<const double> mu,
                                          std::span<const double> covariance)
{
    std::optional<Failure> failure = first_failure(n, mu, covariance);
    if (failure.has_value())
    {
        return std::unexpected(std::move(*failure));
    }
    std::vector<double> symmetric = symmetrized(covariance, n);
    std::vector<double> factor;
    if (!symmetric.empty())
    {
        std::expected<std::vector<double>, std::size_t> computed = linalg::cholesky_or_pivot(symmetric, n);
        if (!computed.has_value())
        {
            return std::unexpected(invalid("covariance is not positive definite: pivot " +
                                           std::to_string(computed.error()) + " is not positive"));
        }
        factor = std::move(*computed);
    }
    return Data(n, std::vector<double>(mu.begin(), mu.end()), std::move(symmetric), std::move(factor));
}

Data::Data(std::size_t n, std::vector<double> mu, std::vector<double> covariance, std::vector<double> factor)
    : n_(n), mu_(std::move(mu)), covariance_(std::move(covariance)), covariance_factor_(std::move(factor))
{
}

std::size_t Data::asset_count() const
{
    return n_;
}

std::size_t Data::weight_count() const
{
    if (risk_free_.has_value())
    {
        return n_ + 1;
    }
    return n_;
}

bool Data::has_mu() const
{
    return !mu_.empty();
}

bool Data::has_covariance() const
{
    return !covariance_.empty();
}

std::span<const double> Data::mu() const
{
    return mu_;
}

std::span<const double> Data::covariance() const
{
    return covariance_;
}

std::span<const double> Data::covariance_factor() const
{
    return covariance_factor_;
}

std::optional<double> Data::risk_free() const
{
    return risk_free_;
}

void Data::set_risk_free(double rate)
{
    risk_free_ = rate;
}

std::expected<void, Failure> Data::set_returns(std::span<const double> returns,
                                               std::span<const double> benchmark, std::size_t t)
{
    std::optional<Failure> failure = check_return_sizes(n_, returns, benchmark, t).or_else([&] {
        return check_return_values(n_, returns, benchmark);
    });
    if (failure.has_value())
    {
        return std::unexpected(std::move(*failure));
    }
    std::vector<double> new_returns(returns.begin(), returns.end());
    std::vector<double> new_benchmark(benchmark.begin(), benchmark.end());
    returns_ = std::move(new_returns);
    benchmark_ = std::move(new_benchmark);
    observation_count_ = t;
    return {};
}

bool Data::has_returns() const
{
    return observation_count_ > 0;
}

std::size_t Data::observation_count() const
{
    return observation_count_;
}

std::span<const double> Data::returns() const
{
    return returns_;
}

std::span<const double> Data::benchmark() const
{
    return benchmark_;
}

} // namespace po
