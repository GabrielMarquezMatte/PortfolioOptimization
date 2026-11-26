#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include "core/portable_math.hpp"

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

// dot(a, r) for a row r whose only nonzero entry is r_column = entry, with the same bits for a finite a: the
// terms a_i * 0 that dot() adds are zeros, which leave its sum unchanged, and + 0.0 turns a product of -0
// into the +0 that dot() returns for it.
inline double unit_row_dot(std::span<const double> a, std::size_t column, double entry)
{
    return a[column] * entry + 0.0;
}

inline std::span<double> basis_column(ActiveSet& set, std::size_t column)
{
    return std::span<double>(set.basis).subspan(column * set.n, set.n);
}

// rows (n x n, zero on entry) becomes J0 = L^{-T} by rows. Entry (i, c) is the back substitution of L' J0 = I
// in its own order, (delta_ic - sum over k = i+1..c of L_ki J0_kc) / L_ii with k increasing: the order of
// linalg::solve_lower_transposed on the unit column c, less the terms L_ki * 0 of k > c, which leave the sum
// unchanged for a finite factor.
inline void inverse_transposed_by_rows(std::span<const double> factor, std::size_t n, std::span<double> rows)
{
    for (std::size_t i = n; i-- > 0;)
    {
        std::span<double> row = rows.subspan(i * n, n);
        row[i] = 1.0;
        for (std::size_t k = i + 1; k < n; ++k)
        {
            double scale = factor[k * n + i];
            std::span<const double> later = rows.subspan(k * n, n);
            for (std::size_t c = k; c < n; ++c)
            {
                row[c] -= scale * later[c];
            }
        }
        for (std::size_t c = i; c < n; ++c)
        {
            row[c] /= factor[i * n + i];
        }
    }
}

inline void transpose_in_place(std::span<double> square, std::size_t n)
{
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = i + 1; j < n; ++j)
        {
            std::swap(square[i * n + j], square[j * n + i]);
        }
    }
}

// J0 = L^{-T} by columns. The entries below its diagonal stay +0, which is 0 / L_ii for the positive diagonal
// of a Cholesky factor.
inline std::vector<double> initial_basis(std::span<const double> factor, std::size_t n)
{
    std::vector<double> basis(n * n, 0.0);
    inverse_transposed_by_rows(factor, n, basis);
    transpose_in_place(basis, n);
    return basis;
}

// No active constraint, with the basis J0 of initial_basis.
inline ActiveSet active_set_from(std::vector<double> basis, std::size_t n)
{
    ActiveSet set;
    set.n = n;
    set.basis = std::move(basis);
    set.triangle.assign(n * n, 0.0);
    set.active.reserve(n);
    set.projection.assign(n, 0.0);
    set.primal_direction.assign(n, 0.0);
    set.dual_direction.assign(n, 0.0);
    return set;
}

inline ActiveSet empty_active_set(std::span<const double> factor, std::size_t n)
{
    return active_set_from(initial_basis(factor, n), n);
}

inline Givens make_givens(double first, double second)
{
    double length = portable::hypot(first, second);
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

inline constexpr std::size_t axpy_columns = 4;
inline constexpr std::size_t dot_columns = 8;

// primal_direction += projection_c J_c for the axpy_columns columns from first, each entry adding them in
// increasing c: the order and the bits of one axpy per column.
inline void add_column_block(ActiveSet& set, std::size_t first)
{
    std::size_t n = set.n;
    std::span<const double> block = std::span<const double>(set.basis).subspan(first * n, axpy_columns * n);
    std::array<double, axpy_columns> scales{};
    std::copy_n(set.projection.begin() + static_cast<std::ptrdiff_t>(first), axpy_columns, scales.begin());
    for (std::size_t i = 0; i < n; ++i)
    {
        double entry = set.primal_direction[i];
        for (std::size_t k = 0; k < axpy_columns; ++k)
        {
            entry += scales[k] * block[k * n + i];
        }
        set.primal_direction[i] = entry;
    }
}

inline void add_column(ActiveSet& set, std::size_t column)
{
    std::span<const double> direction = basis_column(set, column);
    for (std::size_t i = 0; i < set.n; ++i)
    {
        set.primal_direction[i] += set.projection[column] * direction[i];
    }
}

inline void directions_from_projection(ActiveSet& set)
{
    std::ranges::fill(set.primal_direction, 0.0);
    std::size_t column = set.active.size();
    for (; column + axpy_columns <= set.n; column += axpy_columns)
    {
        add_column_block(set, column);
    }
    for (; column < set.n; ++column)
    {
        add_column(set, column);
    }
    back_substitute(set);
}

// dot(J_c, row) for the dot_columns columns from first, each its own sum from +0 with i increasing, the order
// and the bits of dot().
inline std::array<double, dot_columns> dot_column_block(const ActiveSet& set, std::size_t first,
                                                        std::span<const double> row)
{
    std::size_t n = set.n;
    std::span<const double> block = std::span<const double>(set.basis).subspan(first * n, dot_columns * n);
    std::array<double, dot_columns> sums{};
    for (std::size_t i = 0; i < n; ++i)
    {
        double entry = row[i];
        for (std::size_t k = 0; k < dot_columns; ++k)
        {
            sums[k] += block[k * n + i] * entry;
        }
    }
    return sums;
}

inline void compute_directions(ActiveSet& set, std::span<const double> row, double sign)
{
    std::size_t column = 0;
    for (; column + dot_columns <= set.n; column += dot_columns)
    {
        std::array<double, dot_columns> sums = dot_column_block(set, column, row);
        for (std::size_t k = 0; k < dot_columns; ++k)
        {
            set.projection[column + k] = sign * sums[k];
        }
    }
    for (; column < set.n; ++column)
    {
        set.projection[column] = sign * dot(basis_column(set, column), row);
    }
    directions_from_projection(set);
}

// compute_directions for a row whose only nonzero entry is r_column = entry, with the same bits while the
// basis is finite (see unit_row_dot): it reads one entry of each basis column instead of the column.
inline void compute_unit_directions(ActiveSet& set, std::size_t column, double entry, double sign)
{
    std::size_t n = set.n;
    std::span<const double> basis(set.basis);
    for (std::size_t k = 0; k < n; ++k)
    {
        set.projection[k] = sign * unit_row_dot(basis.subspan(k * n, n), column, entry);
    }
    directions_from_projection(set);
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
