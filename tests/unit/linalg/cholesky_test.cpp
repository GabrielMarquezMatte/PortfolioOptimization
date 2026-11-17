#include "linalg/cholesky.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <limits>
#include <optional>
#include <span>
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

// The factorization as it was before it became right-looking: each entry is one dot product of two rows,
// summed in order of k.
bool left_looking_factor(std::span<const double> a, std::size_t n, std::vector<double>& l)
{
    double epsilon = std::numeric_limits<double>::epsilon();
    double relative_floor = po::linalg::kPivotFactor * static_cast<double>(n) * epsilon;
    for (std::size_t j = 0; j < n; ++j)
    {
        double pivot = a[j * n + j];
        for (std::size_t k = 0; k < j; ++k)
        {
            pivot -= l[j * n + k] * l[j * n + k];
        }
        if (!(pivot > relative_floor * a[j * n + j]))
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

std::optional<std::vector<double>> left_looking_cholesky(std::span<const double> a, std::size_t n)
{
    std::optional<std::vector<double>> factor(std::in_place, n * n, 0.0);
    if (!left_looking_factor(a, n, *factor))
    {
        factor.reset();
    }
    return factor;
}

bool same_bits(const std::vector<double>& a, const std::vector<double>& b)
{
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
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

TEST(Cholesky, RejectsPivotAtRoundingLevel)
{
    std::vector<double> a{0.04, 0.04, 0.01, 0.04, 0.04, 0.01, 0.01, 0.01, 0.09};
    EXPECT_FALSE(po::linalg::cholesky(a, 3).has_value());
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

TEST(Cholesky, RightLookingIsBitIdenticalToTheLeftLooking)
{
    for (std::size_t n : {1, 2, 3, 7, 16, 17, 64, 100})
    {
        for (double condition : {1e2, 1e8})
        {
            SCOPED_TRACE(testing::Message() << "n = " << n << ", condition = " << condition);
            std::vector<double> a = po::testing::random_covariance(n, 23 + n, condition);
            std::optional<std::vector<double>> factor = po::linalg::cholesky(a, n);
            std::optional<std::vector<double>> reference = left_looking_cholesky(a, n);
            ASSERT_TRUE(factor.has_value());
            ASSERT_TRUE(reference.has_value());
            EXPECT_TRUE(same_bits(*factor, *reference));
        }
    }
}

} // namespace
