#include "qp/goldfarb_idnani.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <new>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

#include "fuzz/known_shape_fixture.hpp"
#include "fuzz/random_problems.hpp"
#include "fuzz/sign_pattern_reference.hpp"
#include "kkt.hpp"
#include "linalg/cholesky.hpp"
#include "random_covariance.hpp"
#include "user_kkt.hpp"

// GCC cannot see that the replaced operator new below is the one that hands out the pointers freed here.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif

// A failure of the global operator new on demand, for the tests of exception safety: with the countdown at -1
// (its value outside those tests) every allocation goes to malloc, and with it at k >= 0 the allocation after
// the next k throws std::bad_alloc, once.
std::atomic<long> allocations_until_failure{-1};

void* operator new(std::size_t size)
{
    long remaining = allocations_until_failure.load(std::memory_order_relaxed);
    if (remaining >= 0 && allocations_until_failure.fetch_sub(1, std::memory_order_relaxed) == 0)
    {
        allocations_until_failure.store(-1, std::memory_order_relaxed);
        throw std::bad_alloc();
    }
    void* pointer = std::malloc(size == 0 ? 1 : size);
    if (pointer == nullptr)
    {
        throw std::bad_alloc();
    }
    return pointer;
}

void* operator new[](std::size_t size)
{
    return operator new(size);
}

void operator delete(void* pointer) noexcept
{
    std::free(pointer);
}

void operator delete(void* pointer, std::size_t) noexcept
{
    std::free(pointer);
}

void operator delete[](void* pointer) noexcept
{
    std::free(pointer);
}

void operator delete[](void* pointer, std::size_t) noexcept
{
    std::free(pointer);
}

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

// Seed 2: the first whose unconstrained minimizer is negative exactly where the bounds bind.
TEST(GoldfarbIdnani, ActiveBounds)
{
    constexpr std::size_t n = 5;
    Problem problem = unconstrained_problem(n, 2, 1.0);
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

double gross_exposure(std::span<const double> x)
{
    double total = 0.0;
    for (double value : x)
    {
        total += std::fabs(value);
    }
    return total;
}

// 1'x = 1 over 2 to 8 variables with G and g random.
Problem random_simplex_problem(std::uint64_t seed)
{
    std::size_t n = 2 + static_cast<std::size_t>(seed % 7);
    Problem problem = unconstrained_problem(n, seed, 1.0);
    add_equality(problem, std::vector<double>(n, 1.0), 1.0);
    return problem;
}

TEST(GoldfarbIdnani, GrossExposureMatchesSignPatterns)
{
    std::size_t active = 0;
    for (std::uint64_t seed = 1; seed <= 100; ++seed)
    {
        SCOPED_TRACE(seed);
        Problem problem = random_simplex_problem(seed);
        double leverage = po::testing::random_vector(1, seed + 5000, 1.0, 3.0)[0];
        std::optional<double> reference = po::testing::sign_pattern_optimum(problem, leverage);
        ASSERT_TRUE(reference.has_value());
        problem.gross_exposure = po::qp::GrossExposure{std::vector<double>(problem.n, 0.0), -leverage};
        Solution solution = po::qp::solve(problem);
        po::testing::expect_kkt(problem, solution);
        EXPECT_NEAR(po::testing::qp_objective(problem, solution.x), *reference,
                    1e-10 * po::testing::qp_objective_scale(problem, solution.x));
        EXPECT_LE(gross_exposure(solution.x), leverage * (1.0 + 1e-12));
        active += solution.gross_exposure_multipliers.empty() ? 0 : 1;
    }
    EXPECT_GT(active, 50u);
}

TEST(GoldfarbIdnani, GrossExposureInactive)
{
    for (std::uint64_t seed = 1; seed <= 20; ++seed)
    {
        SCOPED_TRACE(seed);
        Problem problem = random_budget_problem(seed);
        Solution without = po::qp::solve(problem);
        problem.gross_exposure = po::qp::GrossExposure{std::vector<double>(problem.n, 0.0), -1e6};
        Solution with = po::qp::solve(problem);
        po::testing::expect_kkt(problem, with);
        EXPECT_EQ(with.x, without.x);
        EXPECT_EQ(with.equality_multipliers, without.equality_multipliers);
        EXPECT_EQ(with.inequality_multipliers, without.inequality_multipliers);
        EXPECT_EQ(with.active_inequalities, without.active_inequalities);
        EXPECT_EQ(with.iterations, without.iterations);
        EXPECT_TRUE(with.gross_exposure_signs.empty());
        EXPECT_TRUE(with.gross_exposure_multipliers.empty());
    }
}

// Under 1'x = 1 the homogeneous form L 1'x - ||x||_1 >= 0 is ||x||_1 <= L: the family with an offset and a zero
// right-hand side gives the optimum of the family with offset 0 and right-hand side -L.
TEST(GoldfarbIdnani, GrossExposureWithOffset)
{
    for (std::uint64_t seed = 1; seed <= 20; ++seed)
    {
        SCOPED_TRACE(seed);
        Problem problem = random_simplex_problem(seed);
        double leverage = po::testing::random_vector(1, seed + 5000, 1.0, 3.0)[0];
        problem.gross_exposure = po::qp::GrossExposure{std::vector<double>(problem.n, 0.0), -leverage};
        Solution plain = po::qp::solve(problem);
        problem.gross_exposure = po::qp::GrossExposure{std::vector<double>(problem.n, leverage), 0.0};
        Solution homogeneous = po::qp::solve(problem);
        po::testing::expect_kkt(problem, homogeneous);
        ASSERT_EQ(plain.status, Status::Optimal);
        double scale = po::testing::largest_magnitude(plain.x);
        for (std::size_t i = 0; i < problem.n; ++i)
        {
            EXPECT_NEAR(homogeneous.x[i], plain.x[i], 1e-12 * scale);
        }
    }
}

// The MAX_RETURN path holds the face 1'x <= L of the family as the equality -1'x = -L (Q(t) at t = 1e10 of a shape
// with cash and nearly equal means, shifted along that face). The steps from the start far away leave a residual on
// the equality well above the tolerance of the final x, and the family, evaluated at sign(x) = 1, regenerated the
// equality's own row: dependent and violated, the problem was declared infeasible. A generated row equal to an
// active explicit row is never generated.
TEST(GoldfarbIdnani, FamilyRowHeldAsEquality)
{
    po::testing::ProblemShape shape{5,    1e6,  1.0, po::testing::ProblemShape::Means::NearlyEqual,
                                    po::testing::ProblemShape::Box::None, true, 5609355878685293435ULL,
                                    1.1024736465921829};
    po::Data data = po::testing::generate(shape).data;
    std::size_t n = shape.n;
    constexpr double t = 1e10;
    Problem problem;
    problem.n = n;
    problem.hessian.assign(data.covariance().begin(), data.covariance().end());
    problem.hessian_factor.assign(data.covariance_factor().begin(), data.covariance_factor().end());
    std::vector<double> excess = po::testing::excess_mu(data);
    auto [low, high] = std::ranges::minmax(excess);
    for (double value : excess)
    {
        problem.linear.push_back(-t * (value - 0.5 * (low + high)));
    }
    problem.equality_rows.assign(n, -1.0);
    problem.equality_rhs.push_back(-*shape.leverage);
    problem.gross_exposure = po::qp::GrossExposure{std::vector<double>(n, 0.0), -*shape.leverage};
    po::testing::expect_kkt(problem, po::qp::solve(problem), 1e-9);
}

// A point of the MAX_RETURN path of the known shape under its Linux covariance (n = 24, cond(Sigma) = 1e8, nearly
// equal means, cash, the leverage only): Q(t) at the t of the climb where 22 generated rows end active, after
// generated rows have entered and left many times. It takes more iterations than a family counted as n rows allows.
TEST(GoldfarbIdnani, GrossExposureFamilyWithManyExchanges)
{
    po::testing::GeneratedProblem fixture = po::testing::fixture_problem(po::testing::known_shape_linux());
    const po::Data& data = fixture.data;
    std::size_t n = data.asset_count();
    constexpr double t = 0x1.0bf0eddf334c2p+43;
    Problem problem;
    problem.n = n;
    problem.hessian.assign(data.covariance().begin(), data.covariance().end());
    problem.hessian_factor.assign(data.covariance_factor().begin(), data.covariance_factor().end());
    for (double value : po::testing::excess_mu(data))
    {
        problem.linear.push_back(-t * value);
    }
    problem.gross_exposure = po::qp::GrossExposure{std::vector<double>(n, 0.0), -*fixture.spec.leverage};
    Solution solution = po::qp::solve(problem);
    ASSERT_EQ(solution.status, Status::Optimal);
    EXPECT_EQ(solution.gross_exposure_multipliers.size(), 22U);
    po::testing::expect_kkt(problem, solution, 1e-9);
}

bool same_bits(const std::vector<double>& a, const std::vector<double>& b)
{
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

void expect_same_solution(const Solution& cached, const Solution& fresh)
{
    EXPECT_EQ(cached.status, fresh.status);
    EXPECT_EQ(cached.iterations, fresh.iterations);
    EXPECT_TRUE(same_bits(cached.x, fresh.x));
    EXPECT_TRUE(same_bits(cached.equality_multipliers, fresh.equality_multipliers));
    EXPECT_TRUE(same_bits(cached.inequality_multipliers, fresh.inequality_multipliers));
    EXPECT_EQ(cached.active_inequalities, fresh.active_inequalities);
    EXPECT_TRUE(same_bits(cached.gross_exposure_signs, fresh.gross_exposure_signs));
    EXPECT_TRUE(same_bits(cached.gross_exposure_multipliers, fresh.gross_exposure_multipliers));
}

// x >= 0 and 1'x = 1 over G = Sigma with its factor, the shape of a MAX_RETURN path.
Problem path_problem(std::size_t n, double condition)
{
    Problem problem;
    problem.n = n;
    problem.hessian = po::testing::random_covariance(n, n + 3, condition);
    problem.hessian_factor = po::linalg::cholesky(problem.hessian, n).value();
    problem.linear.assign(n, 0.0);
    problem.equality_rows.assign(n, 1.0);
    problem.equality_rhs = {1.0};
    problem.inequality_rows.assign(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        problem.inequality_rows[i * n + i] = 1.0;
    }
    problem.inequality_rhs.assign(n, 0.0);
    return problem;
}

// The points of a path share one factor: the linear term moves with t and a row is held as an equality and
// released again. A second problem with another factor (2 L) between them must not reuse the first one's
// basis.
TEST(GoldfarbIdnani, BasisCacheIsBitIdenticalToAFreshBasis)
{
    for (std::size_t n : {1, 2, 3, 7, 16, 17, 64, 100})
    {
        for (double condition : {1e2, 1e8})
        {
            SCOPED_TRACE(testing::Message() << "n = " << n << ", condition = " << condition);
            Problem problem = path_problem(n, condition);
            Problem other = problem;
            for (double& entry : other.hessian_factor)
            {
                entry *= 2.0;
            }
            std::vector<double> excess = po::testing::random_vector(n, n + 5, 0.02, 0.12);
            po::qp::BasisCache cache;
            for (int point = 0; point < 6; ++point)
            {
                double t = 0.05 * point;
                for (std::size_t i = 0; i < n; ++i)
                {
                    problem.linear[i] = -t * excess[i];
                }
                problem.equality_rows.resize(n * (1 + point % 2), 1.0);
                problem.equality_rhs.resize(1 + point % 2, 1.0);
                expect_same_solution(po::qp::solve(problem, cache), po::qp::solve(problem));
                EXPECT_TRUE(same_bits(cache.factor, problem.hessian_factor));
                expect_same_solution(po::qp::solve(problem, cache), po::qp::solve(problem));
                expect_same_solution(po::qp::solve(other, cache), po::qp::solve(other));
                EXPECT_TRUE(same_bits(cache.factor, other.hessian_factor));
            }
        }
    }
}

// A solve that throws while it builds the basis of a new factor must leave the cache as it was, not with the
// new key over the old basis: the next solve with that factor would start from the wrong J0.
TEST(GoldfarbIdnani, BasisCacheKeepsItsPairWhenAnAllocationFails)
{
    std::size_t n = 12;
    Problem problem = path_problem(n, 1e2);
    Problem other = problem;
    for (double& entry : other.hessian_factor)
    {
        entry *= 2.0;
    }
    po::qp::BasisCache before;
    po::qp::solve(problem, before);
    po::qp::BasisCache after_other;
    po::qp::solve(other, after_other);
    int failures = 0;
    bool succeeded = false;
    for (long countdown = 0; countdown < 200; ++countdown)
    {
        SCOPED_TRACE(testing::Message() << "countdown = " << countdown);
        po::qp::BasisCache cache = before;
        allocations_until_failure.store(countdown);
        bool threw = false;
        try
        {
            po::qp::solve(other, cache);
        }
        catch (const std::bad_alloc&)
        {
            threw = true;
        }
        allocations_until_failure.store(-1);
        if (!threw)
        {
            EXPECT_TRUE(same_bits(cache.factor, after_other.factor));
            EXPECT_TRUE(same_bits(cache.basis, after_other.basis));
            succeeded = true;
            break;
        }
        ++failures;
        bool old_pair = same_bits(cache.factor, before.factor) && same_bits(cache.basis, before.basis);
        bool new_pair = same_bits(cache.factor, after_other.factor) && same_bits(cache.basis, after_other.basis);
        EXPECT_TRUE(old_pair || new_pair);
        expect_same_solution(po::qp::solve(other, cache), po::qp::solve(other));
    }
    EXPECT_GT(failures, 0);
    EXPECT_TRUE(succeeded);
}

} // namespace
