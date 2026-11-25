#include "linalg/cholesky.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

namespace po::linalg
{

namespace
{

// Row k above the diagonal holds column k of L while that column updates the trailing block, so that the
// update reads it with unit stride. The row is written before it is read and cleared after, so the part of l
// above the diagonal, a copy of a on entry, is zero at the end.
void stage_column(std::vector<double>& l, std::size_t n, std::size_t k)
{
    for (std::size_t i = k + 1; i < n; ++i)
    {
        l[k * n + i] = l[i * n + k];
    }
}

void clear_staged_column(std::vector<double>& l, std::size_t n, std::size_t k)
{
    auto row = l.begin() + static_cast<std::ptrdiff_t>(k * n);
    std::fill(row + static_cast<std::ptrdiff_t>(k + 1), row + static_cast<std::ptrdiff_t>(n), 0.0);
}

void update_next_column(std::vector<double>& l, std::size_t n, std::size_t k)
{
    double scale = l[(k + 1) * n + k];
    for (std::size_t i = k + 1; i < n; ++i)
    {
        l[i * n + k + 1] -= l[i * n + k] * scale;
    }
}

void update_trailing_block(std::vector<double>& l, std::size_t n, std::size_t k)
{
    std::span<const double> first = std::span<const double>(l).subspan(k * n, n);
    std::span<const double> second = std::span<const double>(l).subspan((k + 1) * n, n);
    for (std::size_t i = k + 2; i < n; ++i)
    {
        double first_scale = l[i * n + k];
        double second_scale = l[i * n + k + 1];
        std::span<double> row = std::span<double>(l).subspan(i * n, n);
        for (std::size_t j = k + 2; j <= i; ++j)
        {
            double entry = row[j];
            entry -= first_scale * first[j];
            entry -= second_scale * second[j];
            row[j] = entry;
        }
    }
}

std::optional<std::size_t> finish_column(std::span<const double> a, std::size_t n, std::vector<double>& l, std::size_t j,
                   double relative_floor)
{
    double pivot = l[j * n + j];
    if (!(pivot > relative_floor * a[j * n + j]))
    {
        return j;
    }
    l[j * n + j] = std::sqrt(pivot);
    for (std::size_t i = j + 1; i < n; ++i)
    {
        l[i * n + j] /= l[j * n + j];
    }
    return std::nullopt;
}

// Right-looking, two columns at a time: once columns k and k + 1 are final, the entries to their right
// subtract their products. Each entry (i, j) still subtracts L_ik L_jk from a_ij in increasing k, the
// order of the dot product of rows i and j in the left-looking form, and is divided by L_jj (or
// square-rooted) after its last term.
std::optional<std::size_t> factor_in_place(std::span<const double> a, std::size_t n, std::vector<double>& l)
{
    double relative_floor = kPivotFactor * static_cast<double>(n) * std::numeric_limits<double>::epsilon();
    for (std::size_t k = 0; k < n; k += 2)
    {
        if (std::optional<std::size_t> rejected = finish_column(a, n, l, k, relative_floor))
        {
            return *rejected;
        }
        if (k + 1 == n)
        {
            return std::nullopt;
        }
        update_next_column(l, n, k);
        if (std::optional<std::size_t> rejected = finish_column(a, n, l, k + 1, relative_floor))
        {
            return *rejected;
        }
        stage_column(l, n, k);
        stage_column(l, n, k + 1);
        update_trailing_block(l, n, k);
        clear_staged_column(l, n, k);
        clear_staged_column(l, n, k + 1);
    }
    return std::nullopt;
}

} // namespace

std::expected<std::vector<double>, std::size_t> cholesky_or_pivot(std::span<const double> a, std::size_t n)
{
    std::span<const double> square = a.first(n * n);
    std::vector<double> factor(square.begin(), square.end());
    if (std::optional<std::size_t> rejected = factor_in_place(a, n, factor))
    {
        return std::unexpected(*rejected);
    }
    return factor;
}

std::optional<std::vector<double>> cholesky(std::span<const double> a, std::size_t n)
{
    std::expected<std::vector<double>, std::size_t> result = cholesky_or_pivot(a, n);
    if (!result.has_value())
    {
        return std::nullopt;
    }
    return std::move(*result);
}

void solve_lower(std::span<const double> l, std::size_t n, std::span<double> b)
{
    for (std::size_t i = 0; i < n; ++i)
    {
        double sum = b[i];
        for (std::size_t k = 0; k < i; ++k)
        {
            sum -= l[i * n + k] * b[k];
        }
        b[i] = sum / l[i * n + i];
    }
}

void solve_lower_transposed(std::span<const double> l, std::size_t n, std::span<double> b)
{
    for (std::size_t i = n; i-- > 0;)
    {
        double sum = b[i];
        for (std::size_t k = i + 1; k < n; ++k)
        {
            sum -= l[k * n + i] * b[k];
        }
        b[i] = sum / l[i * n + i];
    }
}

} // namespace po::linalg
