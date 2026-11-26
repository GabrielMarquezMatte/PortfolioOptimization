#include "estimation/views.hpp"

#include <cmath>
#include <limits>
#include <optional>
#include <string>
#include <utility>

#include "core/data.hpp"
#include "linalg/qr.hpp"

namespace po::estimation
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

std::optional<Failure> check_shape(std::size_t k, std::size_t n)
{
    if (n == 0)
    {
        return invalid("asset count must be positive");
    }
    if (k > n)
    {
        return invalid("k must not exceed n");
    }
    return std::nullopt;
}

std::optional<Failure> check_sizes(std::size_t k, std::size_t n, std::span<const double> p,
                                   std::span<const double> q)
{
    if (p.size() != k * n)
    {
        return invalid("p has " + std::to_string(p.size()) + " entries; expected " + std::to_string(k * n));
    }
    if (q.size() != k)
    {
        return invalid("q has " + std::to_string(q.size()) + " entries; expected " + std::to_string(k));
    }
    return std::nullopt;
}

std::optional<Failure> check_finite(std::span<const double> p, std::span<const double> q, std::size_t n)
{
    for (std::size_t i = 0; i < p.size(); ++i)
    {
        if (!std::isfinite(p[i]))
        {
            return invalid("p is not finite at " + position(i / n, i % n));
        }
    }
    for (std::size_t i = 0; i < q.size(); ++i)
    {
        if (!std::isfinite(q[i]))
        {
            return invalid("q is not finite at " + std::to_string(i));
        }
    }
    return std::nullopt;
}

// P has full row rank when P' (n x k) has full column rank.
std::optional<Failure> check_rank(std::size_t k, std::size_t n, std::span<const double> p)
{
    std::vector<double> transposed(n * k);
    for (std::size_t row = 0; row < k; ++row)
    {
        for (std::size_t col = 0; col < n; ++col)
        {
            transposed[col * k + row] = p[row * n + col];
        }
    }
    if (linalg::has_full_column_rank(linalg::householder_qr(transposed, n, k)))
    {
        return std::nullopt;
    }
    return invalid("views are redundant: P does not have full row rank");
}

std::optional<Failure> check_omega_size(std::span<const double> omega, std::size_t k)
{
    if (omega.size() != k * k)
    {
        return invalid("omega has " + std::to_string(omega.size()) + " entries; expected " +
                       std::to_string(k * k));
    }
    return std::nullopt;
}

std::optional<Failure> check_omega_finite(std::span<const double> omega, std::size_t k)
{
    for (std::size_t i = 0; i < omega.size(); ++i)
    {
        if (!std::isfinite(omega[i]))
        {
            return invalid("omega is not finite at " + position(i / k, i % k));
        }
    }
    return std::nullopt;
}

std::optional<Failure> check_omega_symmetric(std::span<const double> omega, std::size_t k)
{
    for (std::size_t i = 0; i < k; ++i)
    {
        for (std::size_t j = i + 1; j < k; ++j)
        {
            double scale = std::fmax(std::fabs(omega[i * k + i]), std::fabs(omega[j * k + j]));
            if (std::fabs(omega[i * k + j] - omega[j * k + i]) > kSymmetryTolerance * scale)
            {
                return invalid("omega is not symmetric at " + position(i, j));
            }
        }
    }
    return std::nullopt;
}

std::optional<Failure> check_omega_diagonal(std::span<const double> omega, std::size_t k)
{
    for (std::size_t i = 0; i < k; ++i)
    {
        if (omega[i * k + i] < 0.0)
        {
            return invalid("omega has a negative diagonal entry at " + std::to_string(i));
        }
    }
    return std::nullopt;
}

} // namespace

Views::Views(std::size_t k, std::size_t n, std::vector<double> p, std::vector<double> q)
    : k_(k), n_(n), p_(std::move(p)), q_(std::move(q))
{
}

std::expected<Views, Failure> Views::create(std::size_t k, std::size_t n, std::span<const double> p,
                                            std::span<const double> q)
{
    std::optional<Failure> failure = check_shape(k, n)
                                         .or_else([&] { return check_sizes(k, n, p, q); })
                                         .or_else([&] { return check_finite(p, q, n); })
                                         .or_else([&] { return check_rank(k, n, p); });
    if (failure.has_value())
    {
        return std::unexpected(std::move(*failure));
    }
    return Views(k, n, std::vector<double>(p.begin(), p.end()), std::vector<double>(q.begin(), q.end()));
}

std::expected<void, Failure> Views::set_omega(std::span<const double> omega)
{
    std::optional<Failure> failure = check_omega_size(omega, k_)
                                         .or_else([&] { return check_omega_finite(omega, k_); })
                                         .or_else([&] { return check_omega_symmetric(omega, k_); })
                                         .or_else([&] { return check_omega_diagonal(omega, k_); });
    if (failure.has_value())
    {
        return std::unexpected(std::move(*failure));
    }
    omega_.assign(omega.begin(), omega.end());
    return {};
}

std::size_t Views::view_count() const
{
    return k_;
}

std::size_t Views::asset_count() const
{
    return n_;
}

std::span<const double> Views::p() const
{
    return p_;
}

std::span<const double> Views::q() const
{
    return q_;
}

bool Views::has_omega() const
{
    return k_ == 0 || !omega_.empty();
}

std::span<const double> Views::omega() const
{
    return omega_;
}

} // namespace po::estimation
