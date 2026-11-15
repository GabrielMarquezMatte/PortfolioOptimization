#include "linalg/cholesky.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "random_covariance.hpp"

namespace
{

double max_abs(const std::vector<double>& values)
{
    double largest = 0.0;
    for (double value : values)
    {
        largest = std::max(largest, std::fabs(value));
    }
    return largest;
}

std::vector<double> multiply_vector(const std::vector<double>& a, std::size_t n, const std::vector<double>& x)
{
    std::vector<double> product(n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            product[i] += a[i * n + j] * x[j];
        }
    }
    return product;
}

TEST(Cholesky, ReconstructsMatrix)
{
    constexpr std::size_t n = 5;
    std::vector<double> a = po::testing::random_covariance(n, 1);
    auto l = po::linalg::cholesky(a, n);
    ASSERT_TRUE(l.has_value());
    double scale = max_abs(a);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            double sum = 0.0;
            for (std::size_t k = 0; k < n; ++k)
            {
                sum += (*l)[i * n + k] * (*l)[j * n + k];
            }
            EXPECT_NEAR(sum, a[i * n + j], 1e-14 * scale);
        }
    }
}

TEST(Cholesky, RejectsIndefinite)
{
    std::vector<double> a{1.0, 2.0, 2.0, 1.0};
    EXPECT_FALSE(po::linalg::cholesky(a, 2).has_value());
}

TEST(Cholesky, SolvesSystem)
{
    constexpr std::size_t n = 6;
    std::vector<double> a = po::testing::random_covariance(n, 2);
    std::vector<double> b = po::testing::random_vector(n, 3, -1.0, 1.0);
    auto l = po::linalg::cholesky(a, n);
    ASSERT_TRUE(l.has_value());
    std::vector<double> x = b;
    po::linalg::solve_lower(*l, n, x);
    po::linalg::solve_lower_transposed(*l, n, x);
    std::vector<double> residual = multiply_vector(a, n, x);
    for (std::size_t i = 0; i < n; ++i)
    {
        residual[i] -= b[i];
    }
    EXPECT_LE(max_abs(residual), 1e-12 * max_abs(b));
}

TEST(RandomCovariance, IsDeterministicAndSymmetric)
{
    constexpr std::size_t n = 7;
    std::vector<double> first = po::testing::random_covariance(n, 42);
    std::vector<double> second = po::testing::random_covariance(n, 42);
    EXPECT_EQ(first, second);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            EXPECT_EQ(first[i * n + j], first[j * n + i]);
        }
    }
}

} // namespace
