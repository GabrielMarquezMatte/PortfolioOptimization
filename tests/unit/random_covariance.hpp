#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "core/portable_math.hpp"

// The generators of every test instance. They draw from std::mt19937_64, whose output the standard fixes, through
// +, -, *, /, sqrt and po::portable::exp and log only: the <random> distributions are implementation-defined and the
// libm's exp, log, cos and pow differ between platforms, so neither is used, and a seed gives the same bits on
// every platform.
namespace po::testing
{

// 10^x as exp(x ln 10), ln 10 rounded to a double: within a few tens of ulp of 10^x for |x| <= 10, which a
// log-uniform draw does not notice.
inline double exp10(double x)
{
    return po::portable::exp(x * 0x1.26bb1bbb55516p+1);
}

// Uniform in (0, 1): the top 53 bits of the engine plus one half, times 2^-53.
inline double uniform_open(std::mt19937_64& engine)
{
    constexpr double kScale = 1.0 / 9007199254740992.0;
    return (static_cast<double>(engine() >> 11) + 0.5) * kScale;
}

inline double uniform_between(std::mt19937_64& engine, double low, double high)
{
    return low + (high - low) * uniform_open(engine);
}

// Marsaglia's polar method: (u, v) uniform in (-1, 1)^2, drawn again until 0 < s = u^2 + v^2 < 1, then
// u sqrt(-2 log(s) / s). The second normal, v sqrt(-2 log(s) / s), is discarded.
inline double standard_normal(std::mt19937_64& engine)
{
    for (;;)
    {
        double u = 2.0 * uniform_open(engine) - 1.0;
        double v = 2.0 * uniform_open(engine) - 1.0;
        double s = u * u + v * v;
        if (s < 1.0 && s > 0.0)
        {
            return u * std::sqrt(-2.0 * po::portable::log(s) / s);
        }
    }
}

// Uniform in [0, count), count > 0: a draw of the engine below 2^64 mod count is rejected, so that the draws kept
// cover a whole number of copies of [0, count).
inline std::uint64_t uniform_index(std::mt19937_64& engine, std::uint64_t count)
{
    assert(count > 0);
    std::uint64_t rejected = (0 - count) % count;
    for (;;)
    {
        std::uint64_t draw = engine();
        if (draw >= rejected)
        {
            return draw % count;
        }
    }
}

// Uniform in [low, high], low <= high, high - low < 2^64 - 1.
inline std::int64_t uniform_integer(std::mt19937_64& engine, std::int64_t low, std::int64_t high)
{
    assert(low <= high);
    std::uint64_t count = static_cast<std::uint64_t>(high) - static_cast<std::uint64_t>(low) + 1;
    return static_cast<std::int64_t>(static_cast<std::uint64_t>(low) + uniform_index(engine, count));
}

inline std::vector<double> random_vector(std::size_t n, std::uint64_t seed, double low, double high)
{
    std::mt19937_64 engine(seed);
    std::vector<double> values(n);
    for (double& value : values)
    {
        value = uniform_between(engine, low, high);
    }
    return values;
}

inline std::vector<double> normal_matrix(std::size_t n, std::uint64_t seed)
{
    std::mt19937_64 engine(seed);
    std::vector<double> matrix(n * n);
    for (double& value : matrix)
    {
        value = standard_normal(engine);
    }
    return matrix;
}

inline std::vector<double> householder_vector(const std::vector<double>& matrix, std::size_t n, std::size_t k)
{
    std::vector<double> v(n, 0.0);
    double norm_squared = 0.0;
    for (std::size_t i = k; i < n; ++i)
    {
        v[i] = matrix[i * n + k];
        norm_squared += v[i] * v[i];
    }
    double sign = 1.0;
    if (v[k] < 0.0)
    {
        sign = -1.0;
    }
    v[k] += sign * std::sqrt(norm_squared);
    double length = 0.0;
    for (std::size_t i = k; i < n; ++i)
    {
        length += v[i] * v[i];
    }
    length = std::sqrt(length);
    for (std::size_t i = k; i < n; ++i)
    {
        v[i] /= length;
    }
    return v;
}

inline void reflect_rows(std::vector<double>& matrix, const std::vector<double>& v, std::size_t n, std::size_t k)
{
    for (std::size_t j = 0; j < n; ++j)
    {
        double dot = 0.0;
        for (std::size_t i = k; i < n; ++i)
        {
            dot += v[i] * matrix[i * n + j];
        }
        for (std::size_t i = k; i < n; ++i)
        {
            matrix[i * n + j] -= 2.0 * v[i] * dot;
        }
    }
}

inline void reflect_columns(std::vector<double>& matrix, const std::vector<double>& v, std::size_t n, std::size_t k)
{
    for (std::size_t i = 0; i < n; ++i)
    {
        double dot = 0.0;
        for (std::size_t j = k; j < n; ++j)
        {
            dot += matrix[i * n + j] * v[j];
        }
        for (std::size_t j = k; j < n; ++j)
        {
            matrix[i * n + j] -= 2.0 * dot * v[j];
        }
    }
}

inline std::vector<double> orthogonal_factor(std::size_t n, std::uint64_t seed)
{
    std::vector<double> remaining = normal_matrix(n, seed);
    std::vector<double> q(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        q[i * n + i] = 1.0;
    }
    for (std::size_t k = 0; k < n; ++k)
    {
        std::vector<double> v = householder_vector(remaining, n, k);
        reflect_rows(remaining, v, n, k);
        reflect_columns(q, v, n, k);
    }
    return q;
}

inline std::vector<double> random_covariance(std::size_t n, std::uint64_t seed, double condition = 100.0)
{
    std::vector<double> q = orthogonal_factor(n, seed);
    std::vector<double> eigenvalues(n, 0.04);
    for (std::size_t k = 1; k < n; ++k)
    {
        double fraction = static_cast<double>(k) / static_cast<double>(n - 1);
        eigenvalues[k] = 0.04 * po::portable::exp(fraction * po::portable::log(condition));
    }
    std::vector<double> covariance(n * n);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = i; j < n; ++j)
        {
            double sum = 0.0;
            for (std::size_t k = 0; k < n; ++k)
            {
                sum += q[i * n + k] * eigenvalues[k] * q[j * n + k];
            }
            covariance[i * n + j] = sum;
            covariance[j * n + i] = sum;
        }
    }
    return covariance;
}

} // namespace po::testing
