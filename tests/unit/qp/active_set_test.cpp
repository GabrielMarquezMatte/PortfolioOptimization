#include "qp/active_set.hpp"

#include <gtest/gtest.h>

#include <cstddef>
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

} // namespace
