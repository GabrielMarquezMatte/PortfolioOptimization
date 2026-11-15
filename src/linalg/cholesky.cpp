#include "linalg/cholesky.hpp"

#include <cmath>

namespace po::linalg
{

namespace
{

bool factor_in_place(std::span<const double> a, std::size_t n, std::vector<double>& l)
{
    for (std::size_t j = 0; j < n; ++j)
    {
        double pivot = a[j * n + j];
        for (std::size_t k = 0; k < j; ++k)
        {
            pivot -= l[j * n + k] * l[j * n + k];
        }
        if (!(pivot > 0.0))
        {
            return false;
        }
        l[j * n + j] = std::sqrt(pivot);
        for (std::size_t i = j + 1; i < n; ++i)
        {
            double sum = a[i * n + j];
            for (std::size_t k = 0; k < j; ++k)
            {
                sum -= l[i * n + k] * l[j * n + k];
            }
            l[i * n + j] = sum / l[j * n + j];
        }
    }
    return true;
}

} // namespace

std::optional<std::vector<double>> cholesky(std::span<const double> a, std::size_t n)
{
    std::optional<std::vector<double>> factor(std::in_place, n * n, 0.0);
    if (!factor_in_place(a, n, *factor))
    {
        factor.reset();
    }
    return factor;
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
