#pragma once

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace po::linalg
{

std::optional<std::vector<double>> cholesky(std::span<const double> a, std::size_t n);

void solve_lower(std::span<const double> l, std::size_t n, std::span<double> b);

void solve_lower_transposed(std::span<const double> l, std::size_t n, std::span<double> b);

} // namespace po::linalg
