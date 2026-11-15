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
    // An empty mu or covariance is absent. The covariance (n x n, row-major) must be finite, symmetric within
    // 1e-10 of the larger diagonal entry and positive definite; it is stored symmetrized.
    static std::expected<Data, Failure> create(std::size_t n, std::span<const double> mu,
                                               std::span<const double> covariance);

    std::size_t asset_count() const;
    std::size_t weight_count() const;
    bool has_mu() const;
    bool has_covariance() const;
    std::span<const double> mu() const;
    std::span<const double> covariance() const;
    std::optional<double> risk_free() const;
    void set_risk_free(double rate);

private:
    Data(std::size_t n, std::vector<double> mu, std::vector<double> covariance);

    std::size_t n_;
    std::vector<double> mu_;
    std::vector<double> covariance_;
    std::optional<double> risk_free_;
};

} // namespace po
