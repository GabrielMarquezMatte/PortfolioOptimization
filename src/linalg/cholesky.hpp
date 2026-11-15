#pragma once

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace po::linalg
{

// Lower factor L of a = L L' (n x n, row-major), or nullopt when a is not positive definite to working
// precision: a pivot at or below kPivotFactor * n * DBL_EPSILON * a_jj, relative to the original diagonal
// entry a_jj, is rounding noise of a singular matrix and is rejected whatever its sign.
inline constexpr double kPivotFactor = 16.0;

std::optional<std::vector<double>> cholesky(std::span<const double> a, std::size_t n);

void solve_lower(std::span<const double> l, std::size_t n, std::span<double> b);

void solve_lower_transposed(std::span<const double> l, std::size_t n, std::span<double> b);

} // namespace po::linalg
