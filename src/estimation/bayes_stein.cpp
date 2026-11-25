#include "estimation/bayes_stein.hpp"

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

std::optional<Failure> check_observations(std::size_t t, std::size_t n)
{
    if (t > n + 2)
    {
        return std::nullopt;
    }
    return Failure{FailureKind::InvalidArgument, "bayes-stein needs t > n + 2; got t = " + std::to_string(t) +
                                                     " for n = " + std::to_string(n)};
}

// Sigma^-1 x with Sigma = S (t - 1) / (t - n - 2), from the factor of S; the result replaces x.
void solve_sigma(std::span<const double> factor, double scale, std::vector<double>& x)
{
    linalg::solve_lower(factor, x.size(), x);
    linalg::solve_lower_transposed(factor, x.size(), x);
    for (double& value : x)
    {
        value *= scale;
    }
}

double sum(const std::vector<double>& x)
{
    double total = 0.0;
    for (double value : x)
    {
        total += value;
    }
    return total;
}

double dot(const std::vector<double>& a, const std::vector<double>& b)
{
    double total = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        total += a[i] * b[i];
    }
    return total;
}

} // namespace

std::expected<BayesStein, Failure> bayes_stein(const Data& data, std::size_t t)
{
    std::size_t n = data.asset_count();
    std::optional<Failure> failure =
        check_inputs(data, true).or_else([&] { return check_observations(t, n); });
    if (failure.has_value())
    {
        return std::unexpected(std::move(*failure));
    }
    double scale = static_cast<double>(t - n - 2) / static_cast<double>(t - 1);
    std::vector<double> mu(data.mu().begin(), data.mu().end());
    std::vector<double> ones(n, 1.0);
    std::vector<double> solved_mu = mu;
    solve_sigma(data.covariance_factor(), scale, solved_mu);
    solve_sigma(data.covariance_factor(), scale, ones);
    double grand_mean = sum(solved_mu) / sum(ones);
    std::vector<double> deviation(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        deviation[i] = mu[i] - grand_mean;
    }
    std::vector<double> solved = deviation;
    solve_sigma(data.covariance_factor(), scale, solved);
    double shrink = static_cast<double>(n + 2);
    double weight = shrink / (shrink + static_cast<double>(t) * dot(deviation, solved));
    for (std::size_t i = 0; i < n; ++i)
    {
        mu[i] = (1.0 - weight) * mu[i] + weight * grand_mean;
    }
    return BayesStein{std::move(mu), grand_mean, weight};
}

} // namespace po::estimation
