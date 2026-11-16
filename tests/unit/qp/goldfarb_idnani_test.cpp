#include "qp/goldfarb_idnani.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

#include "kkt.hpp"
#include "linalg/cholesky.hpp"
#include "random_covariance.hpp"

namespace
{

using po::qp::Problem;
using po::qp::Solution;
using po::qp::Status;

Problem unconstrained_problem(std::size_t n, std::uint64_t seed, double linear_scale)
{
    Problem problem;
    problem.n = n;
    problem.hessian = po::testing::random_covariance(n, seed);
    problem.linear = po::testing::random_vector(n, seed + 1000, -linear_scale, linear_scale);
    return problem;
}

void scale_hessian(Problem& problem, double factor)
{
    for (double& entry : problem.hessian)
    {
        entry *= factor;
    }
}

void scale_objective(Problem& problem, double factor)
{
    scale_hessian(problem, factor);
    for (double& entry : problem.linear)
    {
        entry *= factor;
    }
}

std::vector<double> solve_hessian(const Problem& problem, std::vector<double> rhs)
{
    std::vector<double> factor = po::linalg::cholesky(problem.hessian, problem.n).value();
    po::linalg::solve_lower(factor, problem.n, rhs);
    po::linalg::solve_lower_transposed(factor, problem.n, rhs);
    return rhs;
}

std::vector<double> unconstrained_minimizer(const Problem& problem)
{
    std::vector<double> negated_linear(problem.n);
    for (std::size_t i = 0; i < problem.n; ++i)
    {
        negated_linear[i] = -problem.linear[i];
    }
    return solve_hessian(problem, negated_linear);
}

std::vector<double> unit_row(std::size_t n, std::size_t index, double value)
{
    std::vector<double> row(n, 0.0);
    row[index] = value;
    return row;
}

void add_equality(Problem& problem, const std::vector<double>& row, double rhs)
{
    problem.equality_rows.insert(problem.equality_rows.end(), row.begin(), row.end());
    problem.equality_rhs.push_back(rhs);
}

void add_inequality(Problem& problem, const std::vector<double>& row, double rhs)
{
    problem.inequality_rows.insert(problem.inequality_rows.end(), row.begin(), row.end());
    problem.inequality_rhs.push_back(rhs);
}

void add_bounds(Problem& problem, const std::vector<double>& lower, const std::vector<double>& upper)
{
    for (std::size_t i = 0; i < problem.n; ++i)
    {
        add_inequality(problem, unit_row(problem.n, i, 1.0), lower[i]);
        add_inequality(problem, unit_row(problem.n, i, -1.0), -upper[i]);
    }
}

void expect_all_nan(const std::vector<double>& values, std::size_t size)
{
    ASSERT_EQ(values.size(), size);
    for (double value : values)
    {
        EXPECT_TRUE(std::isnan(value));
    }
}

TEST(GoldfarbIdnani, Unconstrained)
{
    constexpr std::size_t n = 6;
    Problem problem = unconstrained_problem(n, 3, 1.0);
    std::vector<double> expected = unconstrained_minimizer(problem);
    Solution solution = po::qp::solve(problem);
    ASSERT_EQ(solution.status, Status::Optimal);
    ASSERT_EQ(solution.x.size(), n);
    double scale = po::testing::largest_magnitude(expected);
    for (std::size_t i = 0; i < n; ++i)
    {
        EXPECT_NEAR(solution.x[i], expected[i], 1e-12 * scale);
    }
    po::testing::expect_kkt(problem, solution);
}

TEST(GoldfarbIdnani, EqualityOnly)
{
    constexpr std::size_t n = 4;
    Problem problem = unconstrained_problem(n, 5, 1.0);
    problem.linear.assign(n, 0.0);
    add_equality(problem, std::vector<double>(n, 1.0), 1.0);
    std::vector<double> direction = solve_hessian(problem, std::vector<double>(n, 1.0));
    double total = 0.0;
    for (double value : direction)
    {
        total += value;
    }
    Solution solution = po::qp::solve(problem);
    ASSERT_EQ(solution.status, Status::Optimal);
    ASSERT_EQ(solution.x.size(), n);
    double scale = po::testing::largest_magnitude(direction) / total;
    for (std::size_t i = 0; i < n; ++i)
    {
        EXPECT_NEAR(solution.x[i], direction[i] / total, 1e-12 * scale);
    }
    ASSERT_EQ(solution.equality_multipliers.size(), 1u);
    EXPECT_GT(solution.equality_multipliers[0], 0.0);
    EXPECT_NEAR(solution.equality_multipliers[0], 1.0 / total, 1e-12 / total);
    po::testing::expect_kkt(problem, solution);
}

TEST(GoldfarbIdnani, ActiveBounds)
{
    constexpr std::size_t n = 5;
    Problem problem = unconstrained_problem(n, 3, 1.0);
    const std::vector<double> optimum = {0.3, 0.0, 0.4, 0.0, 0.3};
    const std::vector<double> bound_multipliers = {0.0, 0.05, 0.0, 0.05, 0.0};
    for (std::size_t i = 0; i < n; ++i)
    {
        problem.linear[i] = bound_multipliers[i] - po::testing::row_dot(problem.hessian, i, optimum);
        add_inequality(problem, unit_row(n, i, 1.0), 0.0);
    }
    std::vector<double> unconstrained = unconstrained_minimizer(problem);
    for (std::size_t i = 0; i < n; ++i)
    {
        ASSERT_EQ(unconstrained[i] < 0.0, bound_multipliers[i] > 0.0) << "fixture precondition at " << i;
    }
    Solution solution = po::qp::solve(problem);
    po::testing::expect_kkt(problem, solution);
    ASSERT_EQ(solution.x.size(), n);
    ASSERT_EQ(solution.inequality_multipliers.size(), n);
    for (std::size_t i = 0; i < n; ++i)
    {
        EXPECT_NEAR(solution.x[i], optimum[i], 1e-12);
        EXPECT_NEAR(solution.inequality_multipliers[i], bound_multipliers[i], 1e-12);
    }
    EXPECT_GT(solution.inequality_multipliers[1], 0.0);
    EXPECT_GT(solution.inequality_multipliers[3], 0.0);
}

TEST(GoldfarbIdnani, ReportsActiveInequalities)
{
    constexpr std::size_t n = 5;
    Problem problem = unconstrained_problem(n, 3, 1.0);
    const std::vector<double> optimum = {0.3, 0.0, 0.4, 0.0, 0.3};
    const std::vector<double> bound_multipliers = {0.0, 0.05, 0.0, 0.05, 0.0};
    for (std::size_t i = 0; i < n; ++i)
    {
        problem.linear[i] = bound_multipliers[i] - po::testing::row_dot(problem.hessian, i, optimum);
        add_inequality(problem, unit_row(n, i, -1.0), -1.0);
        add_inequality(problem, unit_row(n, i, 1.0), 0.0);
    }
    Solution solution = po::qp::solve(problem);
    ASSERT_EQ(solution.status, Status::Optimal);
    EXPECT_EQ(solution.active_inequalities, (std::vector<std::size_t>{3, 7}));

    add_equality(problem, std::vector<double>(n, 1.0), 6.0);
    Solution infeasible = po::qp::solve(problem);
    ASSERT_EQ(infeasible.status, Status::Infeasible);
    EXPECT_TRUE(infeasible.active_inequalities.empty());
}

TEST(GoldfarbIdnani, InfeasibleBounds)
{
    constexpr std::size_t n = 3;
    Problem problem = unconstrained_problem(n, 7, 1.0);
    add_equality(problem, std::vector<double>(n, 1.0), 1.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        add_inequality(problem, unit_row(n, i, -1.0), -0.3);
    }
    Solution solution = po::qp::solve(problem);
    EXPECT_EQ(solution.status, Status::Infeasible);
    expect_all_nan(solution.x, n);
}

TEST(GoldfarbIdnani, RedundantConsistentEquality)
{
    constexpr std::size_t n = 4;
    Problem problem = unconstrained_problem(n, 13, 1.0);
    add_equality(problem, std::vector<double>(n, 1.0), 1.0);
    add_equality(problem, std::vector<double>(n, 2.0), 2.0);
    Solution solution = po::qp::solve(problem);
    EXPECT_EQ(solution.status, Status::Optimal);
    po::testing::expect_kkt(problem, solution);
}

TEST(GoldfarbIdnani, InconsistentEqualities)
{
    constexpr std::size_t n = 4;
    Problem problem = unconstrained_problem(n, 17, 1.0);
    add_equality(problem, std::vector<double>(n, 1.0), 1.0);
    add_equality(problem, std::vector<double>(n, 1.0), 2.0);
    Solution solution = po::qp::solve(problem);
    EXPECT_EQ(solution.status, Status::Infeasible);
    expect_all_nan(solution.x, n);
}

TEST(GoldfarbIdnani, EqualBoundsFixVariable)
{
    constexpr std::size_t n = 3;
    Problem problem = unconstrained_problem(n, 19, 1.0);
    add_inequality(problem, unit_row(n, 0, 1.0), 0.2);
    add_inequality(problem, unit_row(n, 0, -1.0), -0.2);
    Solution solution = po::qp::solve(problem);
    po::testing::expect_kkt(problem, solution);
    ASSERT_EQ(solution.x.size(), n);
    EXPECT_NEAR(solution.x[0], 0.2, 1e-14);
}

TEST(GoldfarbIdnani, RejectsIndefiniteHessian)
{
    Problem problem;
    problem.n = 2;
    problem.hessian = {1.0, 2.0, 2.0, 1.0};
    problem.linear = {0.0, 0.0};
    EXPECT_THROW(po::qp::solve(problem), std::invalid_argument);
}

TEST(GoldfarbIdnani, RejectsInconsistentSizes)
{
    Problem problem = unconstrained_problem(3, 23, 1.0);
    problem.equality_rows = {1.0, 1.0};
    problem.equality_rhs = {1.0};
    EXPECT_THROW(po::qp::solve(problem), std::invalid_argument);
    problem = unconstrained_problem(3, 23, 1.0);
    problem.linear.pop_back();
    EXPECT_THROW(po::qp::solve(problem), std::invalid_argument);
}

Problem random_budget_problem(std::uint64_t seed)
{
    std::size_t n = 2 + static_cast<std::size_t>(seed % 29);
    Problem problem = unconstrained_problem(n, seed, 1.0);
    add_equality(problem, std::vector<double>(n, 1.0), 1.0);
    double equal_weight = 1.0 / static_cast<double>(n);
    std::vector<double> lower = po::testing::random_vector(n, seed + 2000, -0.5, equal_weight);
    std::vector<double> upper = po::testing::random_vector(n, seed + 3000, equal_weight, 1.0);
    add_bounds(problem, lower, upper);
    return problem;
}

TEST(GoldfarbIdnani, RandomProblemsSatisfyKkt)
{
    std::size_t active_bounds = 0;
    for (std::uint64_t seed = 1; seed <= 50; ++seed)
    {
        SCOPED_TRACE(seed);
        Problem problem = random_budget_problem(seed);
        Solution solution = po::qp::solve(problem);
        po::testing::expect_kkt(problem, solution);
        for (double multiplier : solution.inequality_multipliers)
        {
            active_bounds += multiplier > 0.0 ? 1 : 0;
        }
    }
    EXPECT_GT(active_bounds, 50u);
}

TEST(GoldfarbIdnani, RedundantConsistentEqualityAtSmallScale)
{
    for (std::uint64_t seed = 1; seed <= 50; ++seed)
    {
        SCOPED_TRACE(seed);
        Problem problem = random_budget_problem(seed);
        add_equality(problem, std::vector<double>(problem.n, 2.0), 2.0);
        scale_objective(problem, 1e-4);
        Solution solution = po::qp::solve(problem);
        EXPECT_EQ(solution.status, Status::Optimal);
        po::testing::expect_kkt(problem, solution);
    }
}

TEST(GoldfarbIdnani, UsesProvidedFactor)
{
    for (std::uint64_t seed = 1; seed <= 10; ++seed)
    {
        SCOPED_TRACE(seed);
        Problem problem = random_budget_problem(seed);
        Solution internal = po::qp::solve(problem);
        problem.hessian_factor = po::linalg::cholesky(problem.hessian, problem.n).value();
        Solution provided = po::qp::solve(problem);
        po::testing::expect_kkt(problem, provided);
        EXPECT_EQ(provided.x, internal.x);
        EXPECT_EQ(provided.equality_multipliers, internal.equality_multipliers);
        EXPECT_EQ(provided.inequality_multipliers, internal.inequality_multipliers);
    }
}

TEST(GoldfarbIdnani, ProvidedFactorIsNotRefactored)
{
    Problem problem;
    problem.n = 2;
    problem.hessian = {1.0, 1.0, 1.0, 1.0};
    problem.linear = {0.0, 0.0};
    EXPECT_THROW(po::qp::solve(problem), std::invalid_argument);
    problem.hessian_factor = {1.0, 0.0, 1.0, 1e-3};
    EXPECT_EQ(po::qp::solve(problem).status, Status::Optimal);
    problem.hessian_factor.pop_back();
    EXPECT_THROW(po::qp::solve(problem), std::invalid_argument);
}

// G = Sigma over n - 1 assets and a nearly riskless one of variance 1e-10, and g = -t mu~ with nearly equal means
// at a daily scale, under the budget and long-only bounds: Q(t) of MAX_RETURN far along its path.
Problem riskless_asset_problem(std::size_t n, std::uint64_t seed, double t)
{
    std::size_t risky = n - 1;
    std::vector<double> covariance = po::testing::random_covariance(risky, seed);
    std::vector<double> noise = po::testing::random_vector(risky, seed + 4000, -1.0, 1.0);
    Problem problem;
    problem.n = n;
    problem.hessian.assign(n * n, 0.0);
    problem.linear.assign(n, -t * 3e-5);
    for (std::size_t i = 0; i < risky; ++i)
    {
        std::copy_n(covariance.begin() + static_cast<std::ptrdiff_t>(i * risky), risky,
                    problem.hessian.begin() + static_cast<std::ptrdiff_t>(i * n));
        problem.linear[i] = -t * (8e-5 + 1e-10 * noise[i]);
    }
    problem.hessian[n * n - 1] = 1e-10;
    add_equality(problem, std::vector<double>(n, 1.0), 1.0);
    add_bounds(problem, std::vector<double>(n, 0.0), std::vector<double>(n, 1.0));
    return problem;
}

// The solver starts from t Sigma^-1 mu~, near 3e17 in the riskless weight, and the rounding of that start leaves
// the active rows off by O(1) until the refinement at the optimum moves x onto them. Every row must hold at that
// final x to the solver's own tolerance, 1e-12 max(1, |c_k|, ||C_k|| |x|max).
TEST(GoldfarbIdnani, FinalFeasibilityAtLargeScale)
{
    constexpr std::size_t n = 30;
    Problem problem = riskless_asset_problem(n, 4, 1e12);
    Solution solution = po::qp::solve(problem);
    ASSERT_EQ(solution.status, Status::Optimal);
    double x_scale = po::testing::largest_magnitude(solution.x);
    for (std::size_t k = 0; k < problem.inequality_rhs.size(); ++k)
    {
        SCOPED_TRACE(k);
        std::span<const double> row = std::span<const double>(problem.inequality_rows).subspan(k * n, n);
        double row_norm = std::sqrt(po::testing::row_dot(problem.inequality_rows, k, row));
        double rhs = problem.inequality_rhs[k];
        double slack = po::testing::row_dot(problem.inequality_rows, k, solution.x) - rhs;
        EXPECT_GE(slack, -1e-12 * std::fmax(1.0, std::fmax(std::fabs(rhs), row_norm * x_scale)));
    }
    po::testing::expect_kkt(problem, solution);
}

TEST(GoldfarbIdnani, RandomProblemsAtExtremeScales)
{
    for (double factor : {1e-4, 1e8})
    {
        for (std::uint64_t seed = 1; seed <= 50; ++seed)
        {
            SCOPED_TRACE(testing::Message() << "factor " << factor << " seed " << seed);
            Problem problem = random_budget_problem(seed);
            scale_objective(problem, factor);
            po::testing::expect_kkt(problem, po::qp::solve(problem));
        }
    }
}

} // namespace
