#include "linalg/qr.hpp"

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

std::vector<double> column(const std::vector<double>& a, std::size_t rows, std::size_t cols, std::size_t j)
{
    std::vector<double> values(rows);
    for (std::size_t i = 0; i < rows; ++i)
    {
        values[i] = a[i * cols + j];
    }
    return values;
}

// a'a, cols x cols.
std::vector<double> gram(const std::vector<double>& a, std::size_t rows, std::size_t cols)
{
    std::vector<double> product(cols * cols, 0.0);
    for (std::size_t i = 0; i < cols; ++i)
    {
        for (std::size_t j = 0; j < cols; ++j)
        {
            for (std::size_t k = 0; k < rows; ++k)
            {
                product[i * cols + j] += a[k * cols + i] * a[k * cols + j];
            }
        }
    }
    return product;
}

// A rows x cols matrix with singular values from 1 down to 1 / condition, from two random orthogonal factors.
std::vector<double> conditioned_matrix(std::size_t rows, std::size_t cols, double condition)
{
    std::vector<double> left = po::testing::orthogonal_factor(rows, 31);
    std::vector<double> right = po::testing::orthogonal_factor(cols, 32);
    std::vector<double> a(rows * cols, 0.0);
    for (std::size_t k = 0; k < cols; ++k)
    {
        double fraction = static_cast<double>(k) / static_cast<double>(cols - 1);
        double singular = std::exp(-fraction * std::log(condition));
        for (std::size_t i = 0; i < rows; ++i)
        {
            for (std::size_t j = 0; j < cols; ++j)
            {
                a[i * cols + j] += left[i * rows + k] * singular * right[j * cols + k];
            }
        }
    }
    return a;
}

TEST(Qr, ReconstructsTheUpperFactor)
{
    constexpr std::size_t rows = 12;
    constexpr std::size_t cols = 5;
    std::vector<double> a = po::testing::random_vector(rows * cols, 1, -1.0, 1.0);
    po::linalg::Qr qr = po::linalg::householder_qr(a, rows, cols);
    for (std::size_t j = 0; j < cols; ++j)
    {
        SCOPED_TRACE(j);
        std::vector<double> rotated = column(a, rows, cols, j);
        po::linalg::apply_transposed(qr, rotated);
        for (std::size_t i = 0; i < rows; ++i)
        {
            double expected = 0.0;
            if (i < cols)
            {
                expected = qr.upper[i * cols + j];
            }
            EXPECT_NEAR(rotated[i], expected, 1e-14);
        }
    }
    for (std::size_t i = 1; i < cols; ++i)
    {
        for (std::size_t j = 0; j < i; ++j)
        {
            EXPECT_EQ(qr.upper[i * cols + j], 0.0);
        }
    }
    EXPECT_TRUE(po::linalg::has_full_column_rank(qr));
}

TEST(Qr, LowerFactorOfTheGram)
{
    constexpr std::size_t rows = 9;
    constexpr std::size_t cols = 6;
    std::vector<double> a = po::testing::random_vector(rows * cols, 2, -1.0, 1.0);
    std::vector<double> lower = po::linalg::lower_factor(po::linalg::householder_qr(a, rows, cols));
    std::vector<double> expected = gram(a, rows, cols);
    for (std::size_t i = 0; i < cols; ++i)
    {
        EXPECT_GT(lower[i * cols + i], 0.0);
        for (std::size_t j = 0; j < cols; ++j)
        {
            double sum = 0.0;
            for (std::size_t k = 0; k < cols; ++k)
            {
                sum += lower[i * cols + k] * lower[j * cols + k];
            }
            if (j > i)
            {
                EXPECT_EQ(lower[i * cols + j], 0.0);
            }
            EXPECT_NEAR(sum, expected[i * cols + j], 1e-13);
        }
    }
}

TEST(Qr, LeastSquaresHasOrthogonalResidual)
{
    constexpr std::size_t rows = 15;
    constexpr std::size_t cols = 4;
    std::vector<double> a = po::testing::random_vector(rows * cols, 3, -1.0, 1.0);
    std::vector<double> b = po::testing::random_vector(rows, 4, -1.0, 1.0);
    std::vector<double> x = po::linalg::least_squares(po::linalg::householder_qr(a, rows, cols), b);
    ASSERT_EQ(x.size(), cols);
    std::vector<double> residual = b;
    for (std::size_t i = 0; i < rows; ++i)
    {
        for (std::size_t j = 0; j < cols; ++j)
        {
            residual[i] -= a[i * cols + j] * x[j];
        }
    }
    for (std::size_t j = 0; j < cols; ++j)
    {
        double normal = 0.0;
        for (std::size_t i = 0; i < rows; ++i)
        {
            normal += a[i * cols + j] * residual[i];
        }
        EXPECT_NEAR(normal, 0.0, 1e-14);
    }
}

TEST(Qr, DetectsRankDeficiency)
{
    constexpr std::size_t rows = 30;
    constexpr std::size_t cols = 4;
    std::vector<double> duplicated = po::testing::random_vector(rows * cols, 5, -0.05, 0.05);
    for (std::size_t i = 0; i < rows; ++i)
    {
        duplicated[i * cols + 3] = duplicated[i * cols + 1];
    }
    EXPECT_FALSE(po::linalg::has_full_column_rank(po::linalg::householder_qr(duplicated, rows, cols)));

    std::vector<double> combination = po::testing::random_vector(rows * cols, 6, -0.05, 0.05);
    for (std::size_t i = 0; i < rows; ++i)
    {
        combination[i * cols + 2] = 0.3 * combination[i * cols + 0] - 2.0 * combination[i * cols + 3];
    }
    EXPECT_FALSE(po::linalg::has_full_column_rank(po::linalg::householder_qr(combination, rows, cols)));

    std::vector<double> zero = po::testing::random_vector(rows * cols, 7, -0.05, 0.05);
    for (std::size_t i = 0; i < rows; ++i)
    {
        zero[i * cols + 0] = 0.0;
    }
    EXPECT_FALSE(po::linalg::has_full_column_rank(po::linalg::householder_qr(zero, rows, cols)));
}

// The rank decision is on R itself: a condition number of 1e10, whose square no double can resolve, is still
// full rank, and 1e15 is rounding.
TEST(Qr, RankThresholdIsRelativeToTheColumns)
{
    constexpr std::size_t rows = 40;
    constexpr std::size_t cols = 8;
    EXPECT_TRUE(po::linalg::has_full_column_rank(
        po::linalg::householder_qr(conditioned_matrix(rows, cols, 1e10), rows, cols)));
    EXPECT_FALSE(po::linalg::has_full_column_rank(
        po::linalg::householder_qr(conditioned_matrix(rows, cols, 1e15), rows, cols)));
    std::vector<double> scaled = conditioned_matrix(rows, cols, 1e10);
    for (double& entry : scaled)
    {
        entry *= 1e-9;
    }
    EXPECT_TRUE(po::linalg::has_full_column_rank(po::linalg::householder_qr(scaled, rows, cols)));
    EXPECT_LT(max_abs(scaled), 1e-8);
}

} // namespace
