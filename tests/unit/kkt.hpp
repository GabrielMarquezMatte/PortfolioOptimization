#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <vector>

#include "qp/problem.hpp"

namespace po::testing
{

struct KktResiduals
{
    double stationarity;
    double primal;
    double dual;
    double complementarity;
};

// NaN must win, or a NaN solution would pass every residual check.
inline double worst_of(double current, double candidate)
{
    if (std::isnan(candidate) || candidate > current)
    {
        return candidate;
    }
    return current;
}

inline double largest_magnitude(std::span<const double> values)
{
    double largest = 0.0;
    for (double value : values)
    {
        largest = worst_of(largest, std::fabs(value));
    }
    return largest;
}

inline double dot(std::span<const double> a, std::span<const double> b)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        sum += a[i] * b[i];
    }
    return sum;
}

inline double row_dot(std::span<const double> rows, std::size_t row, std::span<const double> x)
{
    double sum = 0.0;
    for (std::size_t j = 0; j < x.size(); ++j)
    {
        sum += rows[row * x.size() + j] * x[j];
    }
    return sum;
}

// The magnitude of the computed C_i x: sum_j |C_ij x_j|.
inline double row_magnitude(std::span<const double> rows, std::size_t row, std::span<const double> x)
{
    double sum = 0.0;
    for (std::size_t j = 0; j < x.size(); ++j)
    {
        sum += std::fabs(rows[row * x.size() + j] * x[j]);
    }
    return sum;
}

inline std::vector<double> transposed_product(std::span<const double> rows, std::span<const double> weights,
                                              std::size_t n)
{
    std::vector<double> product(n, 0.0);
    for (std::size_t k = 0; k < weights.size(); ++k)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            product[j] += rows[k * n + j] * weights[k];
        }
    }
    return product;
}

// A scale of 0 only comes from all-zero terms, whose residual is 0 too.
inline double relative_to(double value, double scale)
{
    if (scale == 0.0)
    {
        return value;
    }
    return value / scale;
}

// The term sum_j zeta_j (a - s_j) of the generated rows of the gross exposure family; 0 without the family.
inline std::vector<double> family_term(const po::qp::Problem& problem, const po::qp::Solution& solution)
{
    std::size_t n = problem.n;
    std::vector<double> term(n, 0.0);
    if (!problem.gross_exposure.has_value())
    {
        return term;
    }
    std::span<const double> offset = problem.gross_exposure->offset;
    for (std::size_t j = 0; j < solution.gross_exposure_multipliers.size(); ++j)
    {
        for (std::size_t i = 0; i < n; ++i)
        {
            double normal = offset[i] - solution.gross_exposure_signs[j * n + i];
            term[i] += solution.gross_exposure_multipliers[j] * normal;
        }
    }
    return term;
}

// Stationarity G x + g - E'y - C'z - sum_j zeta_j (a - s_j), relative to the largest of its five terms, with no
// floor: a problem at a small scale is held to the same relative accuracy as one at scale 1.
inline double stationarity_residual(const po::qp::Problem& problem, const po::qp::Solution& solution)
{
    std::size_t n = problem.n;
    std::vector<double> curvature(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        curvature[i] = row_dot(problem.hessian, i, solution.x);
    }
    std::vector<double> equality =
        transposed_product(problem.equality_rows, solution.equality_multipliers, n);
    std::vector<double> inequality =
        transposed_product(problem.inequality_rows, solution.inequality_multipliers, n);
    std::vector<double> family = family_term(problem, solution);
    std::vector<double> residual(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        residual[i] = curvature[i] + problem.linear[i] - equality[i] - inequality[i] - family[i];
    }
    double scale = std::fmax(largest_magnitude(curvature), largest_magnitude(problem.linear));
    scale = std::fmax(scale, std::fmax(largest_magnitude(equality), largest_magnitude(inequality)));
    scale = std::fmax(scale, largest_magnitude(family));
    return relative_to(largest_magnitude(residual), scale);
}

// Complementarity of one row, |z s| / (max(|z|, tiny) row_scale): the slack of a row whose multiplier is not
// zero, relative to the magnitude of the row (row_scale = max(1, |c_i|, sum_j |C_ij x_j|)), and 0 for a row
// whose multiplier is exactly zero. NaN stays NaN.
inline double complementarity_of_row(double multiplier, double slack, double row_scale)
{
    double weight = std::fmax(std::fabs(multiplier), std::numeric_limits<double>::min());
    return std::fabs(multiplier * slack) / (weight * row_scale);
}

// The unscaled residuals of the explicit inequality rows: primal (the largest violation), dual (the most negative
// multiplier, negated) and complementarity (the per-row measure).
inline KktResiduals inequality_residuals(const po::qp::Problem& problem, const po::qp::Solution& solution)
{
    KktResiduals residuals{0.0, 0.0, 0.0, 0.0};
    for (std::size_t k = 0; k < problem.inequality_rhs.size(); ++k)
    {
        double slack = row_dot(problem.inequality_rows, k, solution.x) - problem.inequality_rhs[k];
        double multiplier = solution.inequality_multipliers[k];
        double row_scale = std::fmax(1.0, std::fmax(std::fabs(problem.inequality_rhs[k]),
                                                    row_magnitude(problem.inequality_rows, k, solution.x)));
        residuals.primal = worst_of(residuals.primal, -slack);
        residuals.dual = worst_of(residuals.dual, -multiplier);
        residuals.complementarity =
            worst_of(residuals.complementarity, complementarity_of_row(multiplier, slack, row_scale));
    }
    return residuals;
}

// The same for the gross exposure family a'x - ||x||_1 >= b: primal is its violation b - a'x + ||x||_1, and each
// generated row (a - s_j)'x >= b enters dual and complementarity like an explicit row. Since s_j'x <= ||x||_1, a
// generated row that holds with equality where the family holds has s_j = sign(x) wherever x is not 0.
inline KktResiduals family_residuals(const po::qp::Problem& problem, const po::qp::Solution& solution)
{
    KktResiduals residuals{0.0, 0.0, 0.0, 0.0};
    if (!problem.gross_exposure.has_value())
    {
        return residuals;
    }
    std::size_t n = problem.n;
    const po::qp::GrossExposure& family = *problem.gross_exposure;
    double gross = 0.0;
    for (double value : solution.x)
    {
        gross += std::fabs(value);
    }
    residuals.primal = family.rhs - (dot(family.offset, solution.x) - gross);
    for (std::size_t j = 0; j < solution.gross_exposure_multipliers.size(); ++j)
    {
        std::vector<double> row(family.offset);
        for (std::size_t i = 0; i < n; ++i)
        {
            row[i] -= solution.gross_exposure_signs[j * n + i];
        }
        double slack = dot(row, solution.x) - family.rhs;
        double multiplier = solution.gross_exposure_multipliers[j];
        double magnitude = row_magnitude(row, 0, solution.x);
        double row_scale = std::fmax(1.0, std::fmax(std::fabs(family.rhs), magnitude));
        residuals.dual = worst_of(residuals.dual, -multiplier);
        residuals.complementarity =
            worst_of(residuals.complementarity, complementarity_of_row(multiplier, slack, row_scale));
    }
    return residuals;
}

inline double family_rhs_magnitude(const po::qp::Problem& problem)
{
    if (!problem.gross_exposure.has_value())
    {
        return 0.0;
    }
    return std::fabs(problem.gross_exposure->rhs);
}

inline KktResiduals kkt_residuals(const po::qp::Problem& problem, const po::qp::Solution& solution)
{
    double primal = 0.0;
    for (std::size_t k = 0; k < problem.equality_rhs.size(); ++k)
    {
        double gap = row_dot(problem.equality_rows, k, solution.x) - problem.equality_rhs[k];
        primal = worst_of(primal, std::fabs(gap));
    }
    KktResiduals inequality = inequality_residuals(problem, solution);
    KktResiduals family = family_residuals(problem, solution);
    primal = worst_of(primal, worst_of(inequality.primal, family.primal));
    double rhs_scale = std::fmax(1.0, std::fmax(largest_magnitude(problem.equality_rhs),
                                                largest_magnitude(problem.inequality_rhs)));
    rhs_scale = std::fmax(rhs_scale, family_rhs_magnitude(problem));
    double multiplier_scale = std::fmax(largest_magnitude(solution.inequality_multipliers),
                                        largest_magnitude(solution.gross_exposure_multipliers));
    return KktResiduals{stationarity_residual(problem, solution), primal / rhs_scale,
                        relative_to(worst_of(inequality.dual, family.dual), multiplier_scale),
                        worst_of(inequality.complementarity, family.complementarity)};
}

inline void expect_kkt(const po::qp::Problem& problem, const po::qp::Solution& solution,
                       double tolerance = 1e-10)
{
    ASSERT_EQ(solution.status, po::qp::Status::Optimal);
    ASSERT_EQ(solution.x.size(), problem.n);
    ASSERT_EQ(solution.equality_multipliers.size(), problem.equality_rhs.size());
    ASSERT_EQ(solution.inequality_multipliers.size(), problem.inequality_rhs.size());
    ASSERT_EQ(solution.gross_exposure_signs.size(), solution.gross_exposure_multipliers.size() * problem.n);
    if (!problem.gross_exposure.has_value())
    {
        ASSERT_TRUE(solution.gross_exposure_multipliers.empty());
    }
    KktResiduals residuals = kkt_residuals(problem, solution);
    EXPECT_LE(residuals.stationarity, tolerance);
    EXPECT_LE(residuals.primal, tolerance);
    EXPECT_LE(residuals.dual, tolerance);
    EXPECT_LE(residuals.complementarity, tolerance);
}

} // namespace po::testing
