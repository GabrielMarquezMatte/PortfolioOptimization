#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>

namespace po::qp
{

// A row a = c 1 + (e - c) e_j: every entry is c except at most the one at column j, which is e. The homogenized
// bounds of problems/homogeneous.cpp have this shape (hi_j 1'y - y_j >= 0, y_j - lo_j 1'y >= 0). offset is
// fl(e - c).
struct FlatRow
{
    std::size_t column = 0;
    double flat = 0.0;
    double offset = 0.0;
};

// The sums 1'x and 1'|x| of the point, each in index order, computed once for all the flat rows of a scan.
struct PointSums
{
    double sum = 0.0;
    double absolute_sum = 0.0;
};

inline std::size_t count_unequal(std::span<const double> normal, double value, std::size_t first)
{
    std::size_t count = 0;
    for (std::size_t i = first; i < normal.size(); ++i)
    {
        if (normal[i] != value)
        {
            ++count;
        }
    }
    return count;
}

inline std::size_t first_unequal(std::span<const double> normal, double value)
{
    std::size_t column = 0;
    while (column < normal.size() && normal[column] == value)
    {
        ++column;
    }
    return column;
}

// The row as c 1 + (e - c) e_j, or nothing when two of its entries differ from the rest.
inline std::optional<FlatRow> flat_row_of(std::span<const double> normal)
{
    if (normal.empty())
    {
        return std::nullopt;
    }
    double first = normal[0];
    std::size_t unequal = count_unequal(normal, first, 1);
    if (unequal == 0)
    {
        return FlatRow{0, first, 0.0};
    }
    if (unequal == 1)
    {
        std::size_t column = first_unequal(normal, first);
        return FlatRow{column, first, normal[column] - first};
    }
    if (count_unequal(normal, normal[1], 2) == 0)
    {
        return FlatRow{0, normal[1], first - normal[1]};
    }
    return std::nullopt;
}

inline PointSums point_sums(std::span<const double> x)
{
    PointSums sums;
    for (double value : x)
    {
        sums.sum += value;
        sums.absolute_sum += std::fabs(value);
    }
    return sums;
}

// True only when the slack v = fl(dot(a, x) - b) of the flat row a is proved above -limit without computing it,
// so the scan can skip the row it would skip anyway. With T the exact a'x - b, u = eps / 2 and
// gamma_k = k u / (1 - k u) (Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed., 3.1 and (3.5)):
// - v sums n products and -b in order: |v - T| <= gamma_{n+1} (|a|'|x| + |b|), with
//   |a|'|x| <= |c| 1'|x| + |e - c| |x_j|.
// - w = fl(fl(c fl(1'x)) + fl(d x_j)) - b, with d = fl(e - c) = (e - c)(1 + delta), takes every term through at
//   most n + 3 roundings: |w - T| <= gamma_{n+3} (|c| 1'|x| + |e - c| |x_j| + |b|).
// - |e - c| <= |d| / (1 - u), so |v - w| <= 2 gamma_{n+3} M / (1 - u) with M = |c| 1'|x| + |d| |x_j| + |b|.
// The computed M^ uses fl(1'|x|) >= (1 - gamma_{n-1}) 1'|x| and loses at most gamma_{n+2} more, and n <= 2^32 (the
// n x n Hessian fits in memory), so (2n + 8) eps M^ is nearly twice that bound. A product that underflows adds at
// most 2^-1075 to it, and the n + 4 products of v, w and M^ stay below DBL_MIN / 2, the other term of error.
// Then fl(w - error) > -limit gives w - error > -limit (rounding is monotone and -limit is a double), and
// v >= w - error > -limit. M^ <= 2^1020 keeps every partial sum of both computations finite; a NaN or infinite
// entry of x or of the row fails it.
inline bool holds_by_bound(const FlatRow& row, const PointSums& sums, std::span<const double> x, double rhs,
                           double limit)
{
    constexpr double kLargestMagnitude = 0x1p1020;
    constexpr double kEpsilon = std::numeric_limits<double>::epsilon();
    double entry = x[row.column];
    double magnitude = std::fabs(row.flat) * sums.absolute_sum + std::fabs(row.offset) * std::fabs(entry) +
                       std::fabs(rhs);
    if (!(magnitude <= kLargestMagnitude))
    {
        return false;
    }
    double roundings = static_cast<double>(2 * x.size() + 8);
    double error = roundings * kEpsilon * magnitude + std::numeric_limits<double>::min();
    double estimate = row.flat * sums.sum + row.offset * entry - rhs;
    return estimate - error > -limit;
}

} // namespace po::qp
