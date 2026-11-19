#include "linalg/qr.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
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
        double fraction = 0.0;
        if (cols > 1)
        {
            fraction = static_cast<double>(k) / static_cast<double>(cols - 1);
        }
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

double reference_column_norm(std::span<const double> a, std::size_t rows, std::size_t cols, std::size_t j)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < rows; ++i)
    {
        sum += a[i * cols + j] * a[i * cols + j];
    }
    return std::sqrt(sum);
}

double reference_make_reflector(const std::vector<double>& work, po::linalg::Qr& qr, std::size_t k)
{
    std::size_t cols = qr.cols;
    std::span<double> v = std::span<double>(qr.reflectors).subspan(k * qr.rows, qr.rows);
    double norm_squared = 0.0;
    for (std::size_t i = k; i < qr.rows; ++i)
    {
        v[i] = work[i * cols + k];
        norm_squared += v[i] * v[i];
    }
    double norm = std::sqrt(norm_squared);
    if (norm == 0.0)
    {
        return 0.0;
    }
    double diagonal = -norm;
    if (v[k] < 0.0)
    {
        diagonal = norm;
    }
    v[k] -= diagonal;
    double length = std::sqrt(2.0 * norm * (norm + std::fabs(work[k * cols + k])));
    for (std::size_t i = k; i < qr.rows; ++i)
    {
        v[i] /= length;
    }
    return diagonal;
}

void reference_reflect(std::span<const double> v, std::size_t k, std::span<double> x, std::size_t stride)
{
    double projection = 0.0;
    for (std::size_t i = k; i < v.size(); ++i)
    {
        projection += v[i] * x[i * stride];
    }
    for (std::size_t i = k; i < v.size(); ++i)
    {
        x[i * stride] -= 2.0 * projection * v[i];
    }
}

// The factorization as it was before the reflections ran by rows: each trailing column reflected on its own,
// read with a stride of cols.
po::linalg::Qr column_order_qr(std::span<const double> a, std::size_t rows, std::size_t cols)
{
    po::linalg::Qr qr{rows, cols, std::vector<double>(cols * rows, 0.0), std::vector<double>(cols * cols, 0.0),
                      0.0};
    std::vector<double> work(a.begin(), a.end());
    for (std::size_t j = 0; j < cols; ++j)
    {
        qr.largest_column_norm = std::max(qr.largest_column_norm, reference_column_norm(a, rows, cols, j));
    }
    for (std::size_t k = 0; k < cols; ++k)
    {
        double diagonal = reference_make_reflector(work, qr, k);
        std::span<const double> v = std::span<const double>(qr.reflectors).subspan(k * rows, rows);
        for (std::size_t j = k + 1; j < cols; ++j)
        {
            reference_reflect(v, k, std::span<double>(work).subspan(j), cols);
        }
        qr.upper[k * cols + k] = diagonal;
        for (std::size_t j = k + 1; j < cols; ++j)
        {
            qr.upper[k * cols + j] = work[k * cols + j];
        }
    }
    return qr;
}

void column_order_apply_transposed(const po::linalg::Qr& qr, std::span<double> b)
{
    for (std::size_t k = 0; k < qr.cols; ++k)
    {
        reference_reflect(std::span<const double>(qr.reflectors).subspan(k * qr.rows, qr.rows), k, b, 1);
    }
}

bool same_bits(const std::vector<double>& a, const std::vector<double>& b)
{
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

void expect_same_factorization(const std::vector<double>& a, std::size_t rows, std::size_t cols)
{
    po::linalg::Qr qr = po::linalg::householder_qr(a, rows, cols);
    po::linalg::Qr reference = column_order_qr(a, rows, cols);
    EXPECT_EQ(qr.rows, reference.rows);
    EXPECT_EQ(qr.cols, reference.cols);
    EXPECT_TRUE(same_bits(qr.reflectors, reference.reflectors));
    EXPECT_TRUE(same_bits(qr.upper, reference.upper));
    EXPECT_EQ(std::memcmp(&qr.largest_column_norm, &reference.largest_column_norm, sizeof(double)), 0);
    std::vector<double> b = po::testing::random_vector(rows, 41 + rows, -1.0, 1.0);
    std::vector<double> reference_b = b;
    po::linalg::apply_transposed(qr, b);
    column_order_apply_transposed(reference, reference_b);
    EXPECT_TRUE(same_bits(b, reference_b));
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

TEST(Qr, ReflectIsBitIdenticalToTheColumnOrder)
{
    for (std::size_t n : {1, 2, 3, 7, 16, 17, 64, 100})
    {
        for (double condition : {1e2, 1e8})
        {
            SCOPED_TRACE(testing::Message() << "n = " << n << ", condition = " << condition);
            expect_same_factorization(conditioned_matrix(n, n, condition), n, n);
            expect_same_factorization(conditioned_matrix(2 * n, n, condition), 2 * n, n);
        }
    }
}

// Q [R; 0] is a, column by column, and Q undoes Q'.
TEST(Qr, ApplyRebuildsTheColumns)
{
    constexpr std::size_t rows = 10;
    constexpr std::size_t cols = 4;
    std::vector<double> a = po::testing::random_vector(rows * cols, 5, -1.0, 1.0);
    po::linalg::Qr qr = po::linalg::householder_qr(a, rows, cols);
    for (std::size_t j = 0; j < cols; ++j)
    {
        SCOPED_TRACE(j);
        std::vector<double> rebuilt(rows, 0.0);
        for (std::size_t i = 0; i <= j; ++i)
        {
            rebuilt[i] = qr.upper[i * cols + j];
        }
        po::linalg::apply(qr, rebuilt);
        std::vector<double> expected = column(a, rows, cols, j);
        for (std::size_t i = 0; i < rows; ++i)
        {
            EXPECT_NEAR(rebuilt[i], expected[i], 1e-14);
        }
    }
    std::vector<double> b = po::testing::random_vector(rows, 6, -1.0, 1.0);
    std::vector<double> round_trip = b;
    po::linalg::apply_transposed(qr, round_trip);
    po::linalg::apply(qr, round_trip);
    for (std::size_t i = 0; i < rows; ++i)
    {
        EXPECT_NEAR(round_trip[i], b[i], 1e-15);
    }
}

// R x = b and R'y = b, checked by multiplying back.
TEST(Qr, SolvesWithTheUpperFactorAndItsTranspose)
{
    constexpr std::size_t rows = 8;
    constexpr std::size_t cols = 5;
    po::linalg::Qr qr =
        po::linalg::householder_qr(po::testing::random_vector(rows * cols, 7, -1.0, 1.0), rows, cols);
    std::vector<double> b = po::testing::random_vector(cols, 8, -1.0, 1.0);
    std::vector<double> x = b;
    po::linalg::solve_upper(qr, x);
    std::vector<double> y = b;
    po::linalg::solve_upper_transposed(qr, y);
    for (std::size_t i = 0; i < cols; ++i)
    {
        SCOPED_TRACE(i);
        double upper_times_x = 0.0;
        double transposed_times_y = 0.0;
        for (std::size_t k = 0; k < cols; ++k)
        {
            upper_times_x += qr.upper[i * cols + k] * x[k];
            transposed_times_y += qr.upper[k * cols + i] * y[k];
        }
        EXPECT_NEAR(upper_times_x, b[i], 1e-13);
        EXPECT_NEAR(transposed_times_y, b[i], 1e-13);
    }
}

} // namespace
