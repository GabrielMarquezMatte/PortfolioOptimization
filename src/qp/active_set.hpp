#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <vector>

#include "linalg/cholesky.hpp"

namespace po::qp
{

// index numbers the explicit rows, equalities first, and after them the generated rows of an implicit family.
struct ActiveConstraint
{
    std::size_t index = 0;
    double sign = 1.0;
    double multiplier = 0.0;
};

// Goldfarb and Idnani (1983): for the q active signed normals N, J' N = [R; 0] with R upper triangular
// and J J' = G^{-1}. basis holds J by columns; triangle holds R by rows in its leading q x q block. For the
// normal n+ being added, projection = J' n+, primal_direction = J2 J2' n+, dual_direction = R^{-1} J1' n+.
struct ActiveSet
{
    std::size_t n = 0;
    std::vector<double> basis;
    std::vector<double> triangle;
    std::vector<ActiveConstraint> active;
    std::vector<double> projection;
    std::vector<double> primal_direction;
    std::vector<double> dual_direction;
};

struct Givens
{
    double cosine = 1.0;
    double sine = 0.0;
};

inline double dot(std::span<const double> a, std::span<const double> b)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        sum += a[i] * b[i];
    }
    return sum;
}

inline std::span<double> basis_column(ActiveSet& set, std::size_t column)
{
    return std::span<double>(set.basis).subspan(column * set.n, set.n);
}

inline ActiveSet empty_active_set(std::span<const double> factor, std::size_t n)
{
    ActiveSet set;
    set.n = n;
    set.basis.assign(n * n, 0.0);
    for (std::size_t column = 0; column < n; ++column)
    {
        std::span<double> unit = basis_column(set, column);
        unit[column] = 1.0;
        linalg::solve_lower_transposed(factor, n, unit);
    }
    set.triangle.assign(n * n, 0.0);
    set.projection.assign(n, 0.0);
    set.primal_direction.assign(n, 0.0);
    set.dual_direction.assign(n, 0.0);
    return set;
}

inline Givens make_givens(double first, double second)
{
    double length = std::hypot(first, second);
    if (length == 0.0)
    {
        return Givens{};
    }
    return Givens{first / length, second / length};
}

inline void rotate(Givens rotation, double& first, double& second)
{
    double rotated_first = rotation.cosine * first + rotation.sine * second;
    second = rotation.cosine * second - rotation.sine * first;
    first = rotated_first;
}

inline void rotate_basis(ActiveSet& set, std::size_t column, Givens rotation)
{
    std::span<double> first = basis_column(set, column);
    std::span<double> second = basis_column(set, column + 1);
    for (std::size_t i = 0; i < set.n; ++i)
    {
        rotate(rotation, first[i], second[i]);
    }
}

inline void solve_triangle(const ActiveSet& set, std::span<double> values)
{
    std::size_t n = set.n;
    for (std::size_t i = values.size(); i-- > 0;)
    {
        double sum = values[i];
        for (std::size_t j = i + 1; j < values.size(); ++j)
        {
            sum -= set.triangle[i * n + j] * values[j];
        }
        values[i] = sum / set.triangle[i * n + i];
    }
}

inline void solve_triangle_transposed(const ActiveSet& set, std::span<double> values)
{
    std::size_t n = set.n;
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        double sum = values[i];
        for (std::size_t j = 0; j < i; ++j)
        {
            sum -= set.triangle[j * n + i] * values[j];
        }
        values[i] = sum / set.triangle[i * n + i];
    }
}

inline void back_substitute(ActiveSet& set)
{
    std::size_t q = set.active.size();
    std::copy_n(set.projection.begin(), q, set.dual_direction.begin());
    solve_triangle(set, std::span<double>(set.dual_direction).first(q));
}

// One step of iterative refinement for the active residuals r = b+ - N' x. Since N' J1 = R' and
// G J1 = N R^{-1}, x += J1 R^{-T} r restores N' x = b+ and u += R^{-1} R^{-T} r keeps G x + g = N u.
inline void refine_toward_active(ActiveSet& set, std::span<double> residual, std::span<double> x)
{
    solve_triangle_transposed(set, residual);
    for (std::size_t k = 0; k < residual.size(); ++k)
    {
        std::span<const double> column = basis_column(set, k);
        for (std::size_t i = 0; i < x.size(); ++i)
        {
            x[i] += residual[k] * column[i];
        }
    }
    solve_triangle(set, residual);
    for (std::size_t k = 0; k < residual.size(); ++k)
    {
        set.active[k].multiplier += residual[k];
    }
}

inline void compute_directions(ActiveSet& set, std::span<const double> row, double sign)
{
    for (std::size_t column = 0; column < set.n; ++column)
    {
        set.projection[column] = sign * dot(basis_column(set, column), row);
    }
    std::ranges::fill(set.primal_direction, 0.0);
    for (std::size_t column = set.active.size(); column < set.n; ++column)
    {
        std::span<const double> direction = basis_column(set, column);
        for (std::size_t i = 0; i < set.n; ++i)
        {
            set.primal_direction[i] += set.projection[column] * direction[i];
        }
    }
    back_substitute(set);
}

inline void append_active(ActiveSet& set, const ActiveConstraint& constraint)
{
    std::size_t n = set.n;
    std::size_t q = set.active.size();
    for (std::size_t column = n - 1; column > q; --column)
    {
        Givens rotation = make_givens(set.projection[column - 1], set.projection[column]);
        rotate(rotation, set.projection[column - 1], set.projection[column]);
        rotate_basis(set, column - 1, rotation);
    }
    for (std::size_t i = 0; i <= q; ++i)
    {
        set.triangle[i * n + q] = set.projection[i];
    }
    set.active.push_back(constraint);
}

inline void retriangulate_column(ActiveSet& set, std::size_t column, std::size_t width)
{
    std::size_t n = set.n;
    double& below = set.triangle[(column + 1) * n + column];
    Givens rotation = make_givens(set.triangle[column * n + column], below);
    for (std::size_t j = column; j < width; ++j)
    {
        rotate(rotation, set.triangle[column * n + j], set.triangle[(column + 1) * n + j]);
    }
    below = 0.0;
    rotate_basis(set, column, rotation);
}

inline void drop_active(ActiveSet& set, std::size_t position)
{
    std::size_t n = set.n;
    std::size_t q = set.active.size();
    set.active.erase(set.active.begin() + static_cast<std::ptrdiff_t>(position));
    for (std::size_t column = position; column + 1 < q; ++column)
    {
        for (std::size_t i = 0; i <= column + 1; ++i)
        {
            set.triangle[i * n + column] = set.triangle[i * n + column + 1];
        }
    }
    for (std::size_t i = 0; i < q; ++i)
    {
        set.triangle[i * n + q - 1] = 0.0;
    }
    for (std::size_t column = position; column + 1 < q; ++column)
    {
        retriangulate_column(set, column, q - 1);
    }
}

} // namespace po::qp
