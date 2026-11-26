#include "qp/active_set.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <optional>
#include <span>
#include <vector>

#include "linalg/cholesky.hpp"
#include "random_covariance.hpp"

namespace
{

void reference_solve_lower_transposed(std::span<const double> l, std::size_t n, std::span<double> b)
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

// The base as it was built before it was computed by rows: one solve with L' per unit column, every column
// summed over all n rows, zeros included.
po::qp::ActiveSet reference_empty_active_set(std::span<const double> factor, std::size_t n)
{
    po::qp::ActiveSet set;
    set.n = n;
    set.basis.assign(n * n, 0.0);
    for (std::size_t column = 0; column < n; ++column)
    {
        std::span<double> unit = po::qp::basis_column(set, column);
        unit[column] = 1.0;
        reference_solve_lower_transposed(factor, n, unit);
    }
    set.triangle.assign(n * n, 0.0);
    set.projection.assign(n, 0.0);
    set.primal_direction.assign(n, 0.0);
    set.dual_direction.assign(n, 0.0);
    return set;
}

bool same_bits(const std::vector<double>& a, const std::vector<double>& b)
{
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

TEST(ActiveSet, EmptyBaseIsBitIdenticalToTheOldOrder)
{
    for (std::size_t n : {1, 2, 3, 7, 16, 17, 64, 100})
    {
        for (double condition : {1e2, 1e8})
        {
            SCOPED_TRACE(testing::Message() << "n = " << n << ", condition = " << condition);
            std::vector<double> covariance = po::testing::random_covariance(n, 11 + n, condition);
            std::optional<std::vector<double>> factor = po::linalg::cholesky(covariance, n);
            ASSERT_TRUE(factor.has_value());
            po::qp::ActiveSet set = po::qp::empty_active_set(*factor, n);
            po::qp::ActiveSet reference = reference_empty_active_set(*factor, n);
            EXPECT_EQ(set.n, reference.n);
            EXPECT_TRUE(same_bits(set.basis, reference.basis));
            EXPECT_TRUE(same_bits(set.triangle, reference.triangle));
            EXPECT_TRUE(same_bits(set.projection, reference.projection));
            EXPECT_TRUE(same_bits(set.primal_direction, reference.primal_direction));
            EXPECT_TRUE(same_bits(set.dual_direction, reference.dual_direction));
            EXPECT_TRUE(set.active.empty());
        }
    }
}

// A basis rotated by q random active rows, the state compute_directions sees after q additions.
po::qp::ActiveSet rotated_set(std::span<const double> factor, std::size_t n, std::size_t q,
                              std::uint64_t seed)
{
    po::qp::ActiveSet set = po::qp::empty_active_set(factor, n);
    std::vector<double> rows = po::testing::normal_matrix(n, seed);
    for (std::size_t k = 0; k < q; ++k)
    {
        po::qp::compute_directions(set, std::span<const double>(rows).subspan(k * n, n), 1.0);
        po::qp::append_active(set, po::qp::ActiveConstraint{k, 1.0, 0.0});
    }
    return set;
}

// The unit row e_column * entry, with -0 in some of its zeros: dot() adds them like any other zero.
std::vector<double> unit_row(std::size_t n, std::size_t column, double entry)
{
    std::vector<double> row(n, 0.0);
    for (std::size_t i = 0; i < n; i += 3)
    {
        row[i] = -0.0;
    }
    row[column] = entry;
    return row;
}

void expect_unit_directions_match(const po::qp::ActiveSet& start, std::size_t column, double entry,
                                  double sign)
{
    SCOPED_TRACE(testing::Message() << "column = " << column << ", entry = " << entry << ", sign = " << sign);
    po::qp::ActiveSet dense = start;
    po::qp::ActiveSet unit = start;
    po::qp::compute_directions(dense, unit_row(start.n, column, entry), sign);
    po::qp::compute_unit_directions(unit, column, entry, sign);
    EXPECT_TRUE(same_bits(unit.projection, dense.projection));
    EXPECT_TRUE(same_bits(unit.primal_direction, dense.primal_direction));
    EXPECT_TRUE(same_bits(unit.dual_direction, dense.dual_direction));
}

// The basis of q = 0 is J0, upper triangular: its +0 entries meet the entry -1 in products of -0.
TEST(ActiveSet, UnitRowDirectionsAreBitIdenticalToTheDenseRow)
{
    for (std::size_t n : {1, 2, 3, 7, 16, 17, 64, 100})
    {
        for (double condition : {1e2, 1e8})
        {
            std::vector<double> covariance = po::testing::random_covariance(n, 23 + n, condition);
            std::optional<std::vector<double>> factor = po::linalg::cholesky(covariance, n);
            ASSERT_TRUE(factor.has_value());
            for (std::size_t q : {std::size_t{0}, n / 2, n - 1})
            {
                SCOPED_TRACE(testing::Message()
                             << "n = " << n << ", condition = " << condition << ", q = " << q);
                po::qp::ActiveSet set = rotated_set(*factor, n, q, 31 + n + q);
                std::vector<double> entries = po::testing::random_vector(n, 41 + n + q, -3.0, 3.0);
                for (std::size_t column = 0; column < n; ++column)
                {
                    expect_unit_directions_match(set, column, 1.0, 1.0);
                    expect_unit_directions_match(set, column, -1.0, -1.0);
                    expect_unit_directions_match(set, column, entries[column], 1.0);
                }
            }
        }
    }
}

// The directions as they were computed before the columns were blocked: one dot() per projection, one axpy
// per free column.
void reference_compute_directions(po::qp::ActiveSet& set, std::span<const double> row, double sign)
{
    for (std::size_t column = 0; column < set.n; ++column)
    {
        set.projection[column] = sign * po::qp::dot(po::qp::basis_column(set, column), row);
    }
    std::ranges::fill(set.primal_direction, 0.0);
    for (std::size_t column = set.active.size(); column < set.n; ++column)
    {
        std::span<const double> direction = po::qp::basis_column(set, column);
        for (std::size_t i = 0; i < set.n; ++i)
        {
            set.primal_direction[i] += set.projection[column] * direction[i];
        }
    }
    po::qp::back_substitute(set);
}

// A dense row with +0 and -0 among its entries.
std::vector<double> dense_row(std::size_t n, std::uint64_t seed)
{
    std::vector<double> row = po::testing::random_vector(n, seed, -2.0, 2.0);
    for (std::size_t i = 1; i < n; i += 5)
    {
        row[i] = 0.0;
    }
    for (std::size_t i = 3; i < n; i += 5)
    {
        row[i] = -0.0;
    }
    return row;
}

TEST(ActiveSet, BlockedDirectionsAreBitIdenticalToOneColumnAtATime)
{
    for (std::size_t n : {1, 2, 3, 4, 5, 7, 8, 9, 12, 16, 17, 64, 100})
    {
        std::vector<double> covariance = po::testing::random_covariance(n, 61 + n, 1e6);
        std::optional<std::vector<double>> factor = po::linalg::cholesky(covariance, n);
        ASSERT_TRUE(factor.has_value());
        for (std::size_t q : {std::size_t{0}, std::size_t{1}, n / 2, n - 1})
        {
            SCOPED_TRACE(testing::Message() << "n = " << n << ", q = " << q);
            po::qp::ActiveSet start = rotated_set(*factor, n, q, 67 + n + q);
            std::vector<double> row = dense_row(n, 71 + n + q);
            for (double sign : {1.0, -1.0})
            {
                po::qp::ActiveSet blocked = start;
                po::qp::ActiveSet reference = start;
                po::qp::compute_directions(blocked, row, sign);
                reference_compute_directions(reference, row, sign);
                EXPECT_TRUE(same_bits(blocked.projection, reference.projection));
                EXPECT_TRUE(same_bits(blocked.primal_direction, reference.primal_direction));
                EXPECT_TRUE(same_bits(blocked.dual_direction, reference.dual_direction));
            }
        }
    }
}

TEST(ActiveSet, UnitRowDotIsBitIdenticalToTheDenseDot)
{
    for (std::size_t n : {1, 2, 3, 7, 16, 17, 64, 100})
    {
        SCOPED_TRACE(testing::Message() << "n = " << n);
        std::vector<double> x = po::testing::random_vector(n, 53 + n, -1.0, 1.0);
        for (std::size_t i = 0; i < n; i += 4)
        {
            x[i] = 0.0;
        }
        for (std::size_t i = 2; i < n; i += 4)
        {
            x[i] = -0.0;
        }
        for (std::size_t column = 0; column < n; ++column)
        {
            for (double entry : {1.0, -1.0, 0.5, -2.75})
            {
                double dense = po::qp::dot(unit_row(n, column, entry), x);
                double unit = po::qp::unit_row_dot(x, column, entry);
                EXPECT_EQ(std::memcmp(&dense, &unit, sizeof(double)), 0) << column << " " << entry;
            }
        }
    }
}

} // namespace
