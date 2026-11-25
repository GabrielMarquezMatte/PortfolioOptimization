#include "estimation/posterior.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <utility>

#include "estimation/inputs.hpp"
#include "linalg/cholesky.hpp"

namespace po::estimation
{
namespace
{

std::optional<Failure> check_omega(const Views& views)
{
    if (views.has_omega())
    {
        return std::nullopt;
    }
    return Failure{FailureKind::InvalidArgument, "omega is not set"};
}

std::optional<Failure> check_posterior_inputs(const Data& data, const Views& views, double tau)
{
    return check_inputs(data, true)
        .or_else([&] { return check_views(data, views); })
        .or_else([&] { return check_positive(tau, "tau must be positive and finite"); })
        .or_else([&] { return check_omega(views); });
}

// B = tau Sigma P' (n x k).
std::vector<double> view_covariance(const Data& data, const Views& views, double tau)
{
    std::size_t n = data.asset_count();
    std::size_t k = views.view_count();
    std::span<const double> sigma = data.covariance();
    std::span<const double> p = views.p();
    std::vector<double> b(n * k);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t r = 0; r < k; ++r)
        {
            double total = 0.0;
            for (std::size_t j = 0; j < n; ++j)
            {
                total += sigma[i * n + j] * p[r * n + j];
            }
            b[i * k + r] = tau * total;
        }
    }
    return b;
}

// A = P B + Omega (k x k).
std::vector<double> view_system(const Views& views, std::span<const double> b)
{
    std::size_t n = views.asset_count();
    std::size_t k = views.view_count();
    std::span<const double> p = views.p();
    std::span<const double> omega = views.omega();
    std::vector<double> a(k * k);
    for (std::size_t r = 0; r < k; ++r)
    {
        for (std::size_t s = 0; s < k; ++s)
        {
            double total = 0.0;
            for (std::size_t i = 0; i < n; ++i)
            {
                total += p[r * n + i] * b[i * k + s];
            }
            a[r * k + s] = total + omega[r * k + s];
        }
    }
    return a;
}

// mu = Pi + B (L L')^-1 (Q - P Pi).
std::vector<double> posterior_mean(const Data& data, const Views& views, std::span<const double> b,
                                   std::span<const double> factor)
{
    std::size_t n = data.asset_count();
    std::size_t k = views.view_count();
    std::span<const double> pi = data.mu();
    std::span<const double> p = views.p();
    std::vector<double> surprise(k);
    for (std::size_t r = 0; r < k; ++r)
    {
        double implied = 0.0;
        for (std::size_t i = 0; i < n; ++i)
        {
            implied += p[r * n + i] * pi[i];
        }
        surprise[r] = views.q()[r] - implied;
    }
    linalg::solve_lower(factor, k, surprise);
    linalg::solve_lower_transposed(factor, k, surprise);
    std::vector<double> mu(pi.begin(), pi.end());
    for (std::size_t i = 0; i < n; ++i)
    {
        double shift = 0.0;
        for (std::size_t r = 0; r < k; ++r)
        {
            shift += b[i * k + r] * surprise[r];
        }
        mu[i] += shift;
    }
    return mu;
}

// X = L^-1 B' (k x n), held by column: row i of B, solved in place, is column i of X.
std::vector<double> whitened_views(std::span<const double> b, std::span<const double> factor, std::size_t n,
                                   std::size_t k)
{
    std::vector<double> x(b.begin(), b.end());
    for (std::size_t i = 0; i < n; ++i)
    {
        linalg::solve_lower(factor, k, std::span<double>(x).subspan(i * k, k));
    }
    return x;
}

// Sigma_post = (1 + tau) Sigma - X' X, built on the upper triangle and mirrored, so it is exactly symmetric.
std::vector<double> posterior_covariance(const Data& data, const Views& views, std::span<const double> b,
                                         std::span<const double> factor, double tau)
{
    std::size_t n = data.asset_count();
    std::size_t k = views.view_count();
    std::span<const double> sigma = data.covariance();
    std::vector<double> x = whitened_views(b, factor, n, k);
    double scale = 1.0 + tau;
    std::vector<double> covariance(n * n);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = i; j < n; ++j)
        {
            double explained = 0.0;
            for (std::size_t r = 0; r < k; ++r)
            {
                explained += x[i * k + r] * x[j * k + r];
            }
            covariance[i * n + j] = scale * sigma[i * n + j] - explained;
            covariance[j * n + i] = covariance[i * n + j];
        }
    }
    return covariance;
}

} // namespace

std::expected<Posterior, Failure> posterior(const Data& data, const Views& views, double tau, bool with_covariance)
{
    std::optional<Failure> failure = check_posterior_inputs(data, views, tau);
    if (failure.has_value())
    {
        return std::unexpected(std::move(*failure));
    }
    std::vector<double> b = view_covariance(data, views, tau);
    std::expected<std::vector<double>, std::size_t> factor =
        linalg::cholesky_or_pivot(view_system(views, b), views.view_count());
    if (!factor.has_value())
    {
        return std::unexpected(Failure{FailureKind::InvalidArgument,
                                       "P tau Sigma P' + omega is not positive definite: pivot " +
                                           std::to_string(factor.error()) + " is not positive"});
    }
    std::vector<double> mu = posterior_mean(data, views, b, *factor);
    if (!with_covariance)
    {
        return Posterior{std::move(mu), {}};
    }
    return Posterior{std::move(mu), posterior_covariance(data, views, b, *factor, tau)};
}

} // namespace po::estimation
