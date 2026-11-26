#include "qp/flat_row.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <vector>

#include "qp/active_set.hpp"
#include "qp/goldfarb_idnani.hpp"
#include "random_covariance.hpp"

namespace
{

constexpr double kInfinity = std::numeric_limits<double>::infinity();

// The homogenized bound row scale 1'y + sign y_column of problems/homogeneous.cpp, built the same way.
std::vector<double> bound_row(std::size_t n, std::size_t column, double scale, double sign)
{
    std::vector<double> row(n, scale);
    row[column] += sign;
    return row;
}

double max_abs(std::span<const double> values)
{
    double largest = 0.0;
    for (double value : values)
    {
        largest = std::max(largest, std::fabs(value));
    }
    return largest;
}

// The scan's own tolerance and slack for the row.
struct ScannedRow
{
    double limit = 0.0;
    double slack = 0.0;
};

ScannedRow scanned(std::span<const double> row, std::span<const double> x, double rhs)
{
    double norm = std::sqrt(po::qp::dot(row, row));
    return ScannedRow{po::qp::feasibility_tolerance(rhs, norm, max_abs(x)), po::qp::dot(row, x) - rhs};
}

bool screened(std::span<const double> row, std::span<const double> x, double rhs, double limit)
{
    std::optional<po::qp::FlatRow> flat = po::qp::flat_row_of(row);
    return flat.has_value() && po::qp::holds_by_bound(*flat, po::qp::point_sums(x), x, rhs, limit);
}

TEST(FlatRow, ClassifiesTheHomogenizedBoundRows)
{
    std::optional<po::qp::FlatRow> upper = po::qp::flat_row_of(bound_row(5, 3, 0.3, -1.0));
    ASSERT_TRUE(upper.has_value());
    EXPECT_EQ(upper->column, 3U);
    EXPECT_EQ(upper->flat, 0.3);
    EXPECT_EQ(upper->offset, (0.3 - 1.0) - 0.3);

    std::optional<po::qp::FlatRow> first = po::qp::flat_row_of(bound_row(4, 0, -0.1, 1.0));
    ASSERT_TRUE(first.has_value());
    EXPECT_EQ(first->column, 0U);
    EXPECT_EQ(first->flat, -0.1);
    EXPECT_EQ(first->offset, (-0.1 + 1.0) + 0.1);

    std::optional<po::qp::FlatRow> constant = po::qp::flat_row_of(std::vector<double>{2.0, 2.0, 2.0});
    ASSERT_TRUE(constant.has_value());
    EXPECT_EQ(constant->flat, 2.0);
    EXPECT_EQ(constant->offset, 0.0);

    EXPECT_TRUE(po::qp::flat_row_of(std::vector<double>{1.0, -3.0}).has_value());
    EXPECT_FALSE(po::qp::flat_row_of(std::vector<double>{1.0, 2.0, 3.0}).has_value());
    EXPECT_FALSE(po::qp::flat_row_of(std::vector<double>{0.3, 0.3, -0.7, 0.3, 0.5}).has_value());
    EXPECT_FALSE(po::qp::flat_row_of(std::vector<double>{}).has_value());
}

// A bound row that holds with room is screened, and one that is violated is left to the dot.
TEST(FlatRow, ScreensRowsThatHoldAndLeavesViolatedRowsToTheDot)
{
    std::size_t n = 1000;
    std::vector<double> x = po::testing::random_vector(n, 71, 0.0, 1.0);
    x[7] = 400.0;
    for (std::size_t column : {std::size_t{0}, std::size_t{7}, std::size_t{500}, n - 1})
    {
        SCOPED_TRACE(testing::Message() << "column = " << column);
        std::vector<double> upper = bound_row(n, column, 0.3, -1.0);
        ScannedRow row = scanned(upper, x, 0.0);
        bool holds = row.slack >= -row.limit;
        EXPECT_EQ(holds, column != 7);
        EXPECT_EQ(screened(upper, x, 0.0, row.limit), holds);
    }
}

// The screen never passes a row whose computed slack is below -limit: the right-hand side is swept by units in the
// last place across the point where the slack meets -limit, with the scan's tolerance and with tolerances far
// below it, where only the error bound separates the two.
void expect_screen_is_sound(std::span<const double> row, std::span<const double> x, std::size_t& screened_count)
{
    double edge = po::qp::dot(row, x) + scanned(row, x, 0.0).limit;
    for (double limit_scale : {1.0, 1e-6, 0.0})
    {
        double rhs = edge;
        for (int step = 0; step < 64; ++step)
        {
            rhs = std::nextafter(rhs, -kInfinity);
        }
        for (int step = 0; step < 128; ++step)
        {
            ScannedRow scan = scanned(row, x, rhs);
            double limit = limit_scale * scan.limit;
            if (screened(row, x, rhs, limit))
            {
                ++screened_count;
                ASSERT_GE(scan.slack, -limit) << "rhs = " << rhs << ", limit = " << limit;
            }
            rhs = std::nextafter(rhs, kInfinity);
        }
    }
}

TEST(FlatRow, NeverScreensARowTheScanWouldTake)
{
    std::size_t screened_count = 0;
    for (std::size_t n : {1, 2, 3, 17, 100, 1000})
    {
        for (double spread : {1.0, 1e-8, 1e8})
        {
            SCOPED_TRACE(testing::Message() << "n = " << n << ", spread = " << spread);
            std::vector<double> x = po::testing::random_vector(n, 83 + n, -spread, spread);
            for (std::size_t column : {std::size_t{0}, n / 2, n - 1})
            {
                expect_screen_is_sound(bound_row(n, column, 0.3, -1.0), x, screened_count);
                expect_screen_is_sound(bound_row(n, column, -0.05, 1.0), x, screened_count);
                expect_screen_is_sound(bound_row(n, column, 1.0 / 3.0, 0.0), x, screened_count);
            }
        }
    }
    EXPECT_GT(screened_count, 0U);
}

TEST(FlatRow, DoesNotScreenAPointOrRowThatIsNotFinite)
{
    std::vector<double> row = bound_row(4, 1, 0.3, -1.0);
    std::vector<double> x{0.1, 0.0, 0.2, 0.3};
    EXPECT_TRUE(screened(row, x, -1.0, 1e-12));
    for (double bad : {kInfinity, -kInfinity, std::numeric_limits<double>::quiet_NaN()})
    {
        std::vector<double> point = x;
        point[2] = bad;
        EXPECT_FALSE(screened(row, point, -1.0, 1e-12));
        std::vector<double> bad_row = row;
        bad_row[1] = bad;
        EXPECT_FALSE(screened(bad_row, x, -1.0, 1e-12));
        EXPECT_FALSE(screened(row, x, bad, 1e-12));
    }
    std::vector<double> huge{0x1p1000, 0x1p1000, 0.0, 0.0};
    EXPECT_FALSE(screened(bound_row(4, 1, 0x1p20, 1.0), huge, -1.0, 1e-12));
}

} // namespace
