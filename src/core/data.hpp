#pragma once

#include <cstddef>
#include <expected>
#include <optional>
#include <span>
#include <vector>

#include "core/failure.hpp"

namespace po
{

class Data
{
public:
    // An empty mu or covariance is absent. The covariance (n x n, row-major) must be finite, symmetric
    // within 1e-10 of the larger diagonal entry and positive definite (linalg::cholesky); it is stored
    // symmetrized, with its Cholesky factor, so that positive definiteness is decided once, here.
    static std::expected<Data, Failure> create(std::size_t n, std::span<const double> mu,
                                               std::span<const double> covariance);

    std::size_t asset_count() const;
    std::size_t weight_count() const;
    bool has_mu() const;
    bool has_covariance() const;
    std::span<const double> mu() const;
    std::span<const double> covariance() const;
    // Lower factor L of the stored covariance, Sigma = L L' (n x n, row-major); empty without a covariance.
    std::span<const double> covariance_factor() const;
    std::optional<double> risk_free() const;
    void set_risk_free(double rate);

    // Copies returns (t x n, row-major) and benchmark (t). Failure when t = 0, a size does not match or a
    // value is not finite; on failure the earlier returns are kept.
    std::expected<void, Failure> set_returns(std::span<const double> returns,
                                             std::span<const double> benchmark, std::size_t t);
    bool has_returns() const;
    std::size_t observation_count() const;
    std::span<const double> returns() const;
    std::span<const double> benchmark() const;

private:
    Data(std::size_t n, std::vector<double> mu, std::vector<double> covariance, std::vector<double> factor);

    std::size_t n_;
    std::vector<double> mu_;
    std::vector<double> covariance_;
    std::vector<double> covariance_factor_;
    std::optional<double> risk_free_;
    std::size_t observation_count_ = 0;
    std::vector<double> returns_;
    std::vector<double> benchmark_;
};

} // namespace po
