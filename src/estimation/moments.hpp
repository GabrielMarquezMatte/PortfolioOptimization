#pragma once

#include <cstddef>
#include <expected>
#include <span>
#include <vector>

#include "core/failure.hpp"

namespace po::estimation
{

struct Moments
{
    std::vector<double> mu;
    std::vector<double> covariance;
};

// Sample mean and unbiased (t - 1) covariance of a row-major t x n return panel.
std::expected<Moments, Failure> sample_moments(std::span<const double> returns, std::size_t t, std::size_t n);

} // namespace po::estimation
