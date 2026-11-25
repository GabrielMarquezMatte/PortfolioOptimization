#pragma once

#include <cstddef>
#include <expected>
#include <span>
#include <vector>

#include "core/failure.hpp"

namespace po::estimation
{

// The k views of a Black-Litterman model: P (k x n, row-major), Q (k entries, in the same unit as the
// equilibrium return, in excess of the risk-free rate) and, once set, Omega (k x k).
class Views
{
public:
    // Needs n >= 1, k <= n, finite p and q, and P of full row rank. k = 0 is the empty set of views, with p
    // and q empty.
    static std::expected<Views, Failure> create(std::size_t k, std::size_t n, std::span<const double> p,
                                                std::span<const double> q);

    // Copies Omega (k x k): finite, symmetric within kSymmetryTolerance of the larger diagonal entry, with a
    // nonnegative diagonal. A failure keeps the previous Omega.
    std::expected<void, Failure> set_omega(std::span<const double> omega);

    std::size_t view_count() const;
    std::size_t asset_count() const;
    std::span<const double> p() const;
    std::span<const double> q() const;
    // True once set_omega has succeeded; vacuously true for k = 0, which has no Omega to set.
    bool has_omega() const;
    // Empty until set_omega succeeds (and for k = 0).
    std::span<const double> omega() const;

private:
    Views(std::size_t k, std::size_t n, std::vector<double> p, std::vector<double> q);

    std::size_t k_;
    std::size_t n_;
    std::vector<double> p_;
    std::vector<double> q_;
    std::vector<double> omega_;
};

} // namespace po::estimation
