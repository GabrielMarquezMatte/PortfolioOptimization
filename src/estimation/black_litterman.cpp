#include "estimation/black_litterman.hpp"

#include <cmath>
#include <optional>
#include <string>
#include <utility>

#include "estimation/inputs.hpp"

namespace po::estimation
{
namespace
{

Failure invalid(std::string message)
{
    return Failure{FailureKind::InvalidArgument, std::move(message)};
}

std::optional<Failure> check_weights(std::span<const double> w_mkt, std::size_t n)
{
    if (w_mkt.size() != n)
    {
        return invalid("w_mkt has " + std::to_string(w_mkt.size()) + " entries; expected " + std::to_string(n));
    }
    double total = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!std::isfinite(w_mkt[i]))
        {
            return invalid("w_mkt is not finite at " + std::to_string(i));
        }
        total += w_mkt[i];
    }
    if (std::fabs(total - 1.0) > kMarketWeightTolerance)
    {
        return invalid("w_mkt must sum to 1 within 1e-8");
    }
    return std::nullopt;
}

std::optional<Failure> check_confidence(std::span<const double> confidence, std::size_t k)
{
    if (confidence.size() != k)
    {
        return invalid("confidence has " + std::to_string(confidence.size()) + " entries; expected " +
                       std::to_string(k));
    }
    for (std::size_t i = 0; i < k; ++i)
    {
        if (!(confidence[i] > 0.0 && confidence[i] <= 1.0))
        {
            return invalid("confidence must be in (0, 1] at " + std::to_string(i));
        }
    }
    return std::nullopt;
}

// tau p_k Sigma p_k' for each view k.
std::vector<double> view_variances(const Data& data, const Views& views, double tau)
{
    std::size_t n = data.asset_count();
    std::span<const double> covariance = data.covariance();
    std::vector<double> variances(views.view_count());
    for (std::size_t k = 0; k < variances.size(); ++k)
    {
        std::span<const double> row = views.p().subspan(k * n, n);
        double total = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            double sigma_p = 0.0;
            for (std::size_t j = 0; j < n; ++j)
            {
                sigma_p += covariance[i * n + j] * row[j];
            }
            total += row[i] * sigma_p;
        }
        variances[k] = tau * total;
    }
    return variances;
}

std::vector<double> diagonal(const std::vector<double>& entries)
{
    std::size_t k = entries.size();
    std::vector<double> omega(k * k, 0.0);
    for (std::size_t i = 0; i < k; ++i)
    {
        omega[i * k + i] = entries[i];
    }
    return omega;
}

std::optional<Failure> check_omega_inputs(const Data& data, const Views& views, double tau)
{
    return check_inputs(data, false)
        .or_else([&] { return check_views(data, views); })
        .or_else([&] { return check_positive(tau, "tau must be positive and finite"); });
}

} // namespace

std::expected<std::vector<double>, Failure> equilibrium(const Data& data, std::span<const double> w_mkt,
                                                        double delta)
{
    std::size_t n = data.asset_count();
    std::optional<Failure> failure =
        check_inputs(data, false)
            .or_else([&] { return check_weights(w_mkt, n); })
            .or_else([&] { return check_positive(delta, "delta must be positive and finite"); });
    if (failure.has_value())
    {
        return std::unexpected(std::move(*failure));
    }
    std::span<const double> covariance = data.covariance();
    std::vector<double> pi(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        double total = 0.0;
        for (std::size_t j = 0; j < n; ++j)
        {
            total += covariance[i * n + j] * w_mkt[j];
        }
        pi[i] = delta * total;
    }
    return pi;
}

std::expected<std::vector<double>, Failure> omega_he_litterman(const Data& data, const Views& views, double tau)
{
    std::optional<Failure> failure = check_omega_inputs(data, views, tau);
    if (failure.has_value())
    {
        return std::unexpected(std::move(*failure));
    }
    return diagonal(view_variances(data, views, tau));
}

std::expected<std::vector<double>, Failure> omega_idzorek(const Data& data, const Views& views, double tau,
                                                          std::span<const double> confidence)
{
    std::optional<Failure> failure =
        check_omega_inputs(data, views, tau).or_else([&] { return check_confidence(confidence, views.view_count()); });
    if (failure.has_value())
    {
        return std::unexpected(std::move(*failure));
    }
    std::vector<double> entries = view_variances(data, views, tau);
    for (std::size_t k = 0; k < entries.size(); ++k)
    {
        entries[k] *= (1.0 - confidence[k]) / confidence[k];
    }
    return diagonal(entries);
}

} // namespace po::estimation
