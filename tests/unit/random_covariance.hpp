#pragma once

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numbers>
#include <random>
#include <vector>

namespace po::testing
{

inline double uniform_open(std::mt19937_64& engine)
{
    constexpr double kScale = 1.0 / 9007199254740992.0;
    return (static_cast<double>(engine() >> 11) + 0.5) * kScale;
}

inline double standard_normal(std::mt19937_64& engine)
{
    double radius = std::sqrt(-2.0 * std::log(uniform_open(engine)));
    return radius * std::cos(2.0 * std::numbers::pi * uniform_open(engine));
}

inline std::vector<double> random_vector(std::size_t n, std::uint64_t seed, double low, double high)
{
    std::mt19937_64 engine(seed);
    std::vector<double> values(n);
    for (double& value : values)
    {
        value = low + (high - low) * uniform_open(engine);
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
        eigenvalues[k] = 0.04 * std::exp(fraction * std::log(condition));
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
