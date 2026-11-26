#pragma once

#include <cmath>
#include <cstddef>
#include <optional>
#include <span>
#include <vector>

#include "kkt.hpp"
#include "qp/goldfarb_idnani.hpp"
#include "qp/problem.hpp"

namespace po::testing
{

// 0.5 x'Gx + g'x.
inline double qp_objective(const po::qp::Problem& problem, std::span<const double> x)
{
    std::vector<double> hessian = hessian_of(problem);
    double value = 0.0;
    for (std::size_t i = 0; i < problem.n; ++i)
    {
        value += (0.5 * row_dot(hessian, i, x) + problem.linear[i]) * x[i];
    }
    return value;
}

// |0.5 x'Gx| + |g'x|, the magnitude of the two terms of the objective: a value that cancels to near 0 is compared
// relative to its terms, not to its own rounding.
inline double qp_objective_scale(const po::qp::Problem& problem, std::span<const double> x)
{
    std::vector<double> hessian = hessian_of(problem);
    double curvature = 0.0;
    double linear = 0.0;
    for (std::size_t i = 0; i < problem.n; ++i)
    {
        curvature += 0.5 * row_dot(hessian, i, x) * x[i];
        linear += problem.linear[i] * x[i];
    }
    return std::fabs(curvature) + std::fabs(linear);
}

inline void add_sign_pattern(po::qp::Problem& problem, std::size_t pattern, double leverage)
{
    std::size_t n = problem.n;
    std::vector<double> signs(n, 1.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        if ((pattern >> i) & 1U)
        {
            signs[i] = -1.0;
        }
        std::vector<double> row(n, 0.0);
        row[i] = signs[i];
        problem.inequality_rows.insert(problem.inequality_rows.end(), row.begin(), row.end());
        problem.inequality_rhs.push_back(0.0);
    }
    for (double sign : signs)
    {
        problem.inequality_rows.push_back(-sign);
    }
    problem.inequality_rhs.push_back(-leverage);
}

// The optimum of problem_without_family under sum_i |x_i| <= leverage, by brute force for n <= 8: on the orthant of
// a sign pattern s, s_i x_i >= 0, the constraint is the linear row s'x <= leverage, so the optimum is the best
// Optimal value over the 2^n patterns. nullopt when no pattern is Optimal. An independent route to the
// leverage-constrained optimum: no generated row, only explicit ones.
inline std::optional<double> sign_pattern_optimum(const po::qp::Problem& problem_without_family,
                                                  double leverage)
{
    std::optional<double> best;
    for (std::size_t pattern = 0; pattern < (std::size_t{1} << problem_without_family.n); ++pattern)
    {
        po::qp::Problem orthant = problem_without_family;
        add_sign_pattern(orthant, pattern, leverage);
        po::qp::Solution solution = po::qp::solve(orthant);
        if (solution.status != po::qp::Status::Optimal)
        {
            continue;
        }
        double value = qp_objective(orthant, solution.x);
        if (!best.has_value() || value < *best)
        {
            best = value;
        }
    }
    return best;
}

} // namespace po::testing
