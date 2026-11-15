#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
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

inline double row_dot(std::span<const double> rows, std::size_t row, std::span<const double> x)
{
    double sum = 0.0;
    for (std::size_t j = 0; j < x.size(); ++j)
    {
        sum += rows[row * x.size() + j] * x[j];
    }
    return sum;
}

inline void subtract_transposed(std::span<const double> rows, std::span<const double> weights,
                                std::vector<double>& residual)
{
    std::size_t n = residual.size();
    for (std::size_t k = 0; k < weights.size(); ++k)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            residual[j] -= rows[k * n + j] * weights[k];
        }
    }
}

inline double stationarity_residual(const po::qp::Problem& problem, const po::qp::Solution& solution)
{
    std::vector<double> residual(problem.n);
    for (std::size_t i = 0; i < problem.n; ++i)
    {
        residual[i] = row_dot(problem.hessian, i, solution.x) + problem.linear[i];
    }
    subtract_transposed(problem.equality_rows, solution.equality_multipliers, residual);
    subtract_transposed(problem.inequality_rows, solution.inequality_multipliers, residual);
    double hessian_scale = largest_magnitude(problem.hessian) * largest_magnitude(solution.x);
    double scale = std::fmax(1.0, std::fmax(hessian_scale, largest_magnitude(problem.linear)));
    return largest_magnitude(residual) / scale;
}

inline KktResiduals kkt_residuals(const po::qp::Problem& problem, const po::qp::Solution& solution)
{
    double primal = 0.0;
    for (std::size_t k = 0; k < problem.equality_rhs.size(); ++k)
    {
        double gap = row_dot(problem.equality_rows, k, solution.x) - problem.equality_rhs[k];
        primal = worst_of(primal, std::fabs(gap));
    }
    double dual = 0.0;
    double complementarity = 0.0;
    for (std::size_t k = 0; k < problem.inequality_rhs.size(); ++k)
    {
        double slack = row_dot(problem.inequality_rows, k, solution.x) - problem.inequality_rhs[k];
        double multiplier = solution.inequality_multipliers[k];
        primal = worst_of(primal, -slack);
        dual = worst_of(dual, -multiplier);
        complementarity = worst_of(complementarity, std::fabs(multiplier * slack));
    }
    double rhs_scale = std::fmax(1.0, std::fmax(largest_magnitude(problem.equality_rhs),
                                                largest_magnitude(problem.inequality_rhs)));
    return KktResiduals{stationarity_residual(problem, solution), primal / rhs_scale, dual / rhs_scale,
                        complementarity / rhs_scale};
}

inline void expect_kkt(const po::qp::Problem& problem, const po::qp::Solution& solution,
                       double tolerance = 1e-10)
{
    ASSERT_EQ(solution.status, po::qp::Status::Optimal);
    ASSERT_EQ(solution.x.size(), problem.n);
    ASSERT_EQ(solution.equality_multipliers.size(), problem.equality_rhs.size());
    ASSERT_EQ(solution.inequality_multipliers.size(), problem.inequality_rhs.size());
    KktResiduals residuals = kkt_residuals(problem, solution);
    EXPECT_LE(residuals.stationarity, tolerance);
    EXPECT_LE(residuals.primal, tolerance);
    EXPECT_LE(residuals.dual, tolerance);
    EXPECT_LE(residuals.complementarity, tolerance);
}

} // namespace po::testing
