#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <vector>

#include "fuzz/random_problems.hpp"
#include "fuzz/ratio_reference.hpp"
#include "kkt.hpp"
#include "linalg/cholesky.hpp"
#include "problems/solve.hpp"
#include "random_covariance.hpp"
#include "user_kkt.hpp"

namespace
{

using po::Data;
using po::Objective;
using po::Result;
using po::SolveStatus;
using po::Spec;
using po::testing::largest_magnitude;
using po::testing::ratio_of;
using po::testing::ratio_vector;

constexpr double kRiskFree = 0.02;
constexpr double kRiskFreeAboveMeans = 0.2;
// The closed forms, relative to the largest weight.
constexpr double kClosedForm = 1e-12;

Data random_data(std::size_t n, std::uint64_t seed)
{
    std::vector<double> mu = po::testing::random_vector(n, seed + 500, 0.03, 0.15);
    return Data::create(n, mu, po::testing::random_covariance(n, seed)).value();
}

std::vector<double> covariance_solve(const Data& data, std::vector<double> rhs)
{
    std::size_t n = data.asset_count();
    std::vector<double> factor = po::linalg::cholesky(data.covariance(), n).value();
    po::linalg::solve_lower(factor, n, rhs);
    po::linalg::solve_lower_transposed(factor, n, rhs);
    return rhs;
}

double sum(std::span<const double> values)
{
    double total = 0.0;
    for (double value : values)
    {
        total += value;
    }
    return total;
}

Spec ratio_spec(Objective objective, double lower, double upper)
{
    Spec spec;
    spec.objective = objective;
    spec.lower = lower;
    spec.upper = upper;
    return spec;
}

// No bounds, the default.
Spec unbounded_spec(Objective objective)
{
    Spec spec;
    spec.objective = objective;
    return spec;
}

// OPTIMAL with KKT in the user's variables for f = -v'w / vol(w), the gradient written here.
Result solve_checked(const Spec& spec, const Data& data)
{
    Result result = po::solve(spec, data).value();
    std::vector<double> v = ratio_vector(data, spec.objective);
    po::testing::expect_user_kkt(po::testing::with_cash_fixed(spec, data), data, result,
                                 po::testing::ratio_gradient(data, v, result));
    return result;
}

// The closed form without bounds: w = Sigma^-1 v / (1' Sigma^-1 v) over the risky weights.
void expect_closed_form(const Data& data, const Result& result, std::span<const double> v)
{
    std::vector<double> direction = covariance_solve(data, std::vector<double>(v.begin(), v.end()));
    double total = sum(direction);
    ASSERT_GT(total, 0.0);
    for (double& entry : direction)
    {
        entry /= total;
    }
    double scale = largest_magnitude(direction);
    for (std::size_t i = 0; i < direction.size(); ++i)
    {
        EXPECT_NEAR(result.weights[i], direction[i], kClosedForm * scale);
    }
}

void expect_no_point(const Result& result, SolveStatus status)
{
    EXPECT_EQ(result.status, status);
    for (double weight : result.weights)
    {
        EXPECT_TRUE(std::isnan(weight));
    }
    for (double multiplier : result.bound_multipliers)
    {
        EXPECT_TRUE(std::isnan(multiplier));
    }
    EXPECT_TRUE(std::isnan(result.budget_multiplier));
    EXPECT_TRUE(std::isnan(result.leverage_multiplier));
    EXPECT_TRUE(std::isnan(result.residual));
}

void expect_cash_at_zero(const Data& data, const Result& result)
{
    ASSERT_EQ(result.weights.size(), data.asset_count() + 1);
    EXPECT_EQ(result.weights.back(), 0.0);
    EXPECT_TRUE(std::isnan(result.budget_multiplier));
    EXPECT_TRUE(std::isfinite(result.bound_multipliers.back()));
}

TEST(MaxSharpe, UnconstrainedIsTangency)
{
    for (std::uint64_t seed = 1; seed <= 5; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = random_data(6, seed);
        Result result = solve_checked(unbounded_spec(Objective::MaxSharpe), data);
        expect_closed_form(data, result, data.mu());
    }
}

// With rf the tangency of mu - rf 1, the cash weight exactly 0 and the budget term in beta_0: 0 at the
// tangency, where w'grad f = 0 (f is homogeneous of degree 0) leaves no budget term, and not 0 under an
// active cap.
TEST(MaxSharpe, CashWeightZero)
{
    Data data = random_data(6, 7);
    data.set_risk_free(kRiskFree);
    Result tangency = solve_checked(unbounded_spec(Objective::MaxSharpe), data);
    expect_cash_at_zero(data, tangency);
    expect_closed_form(data, tangency, po::testing::excess_mu(data));
    double multiplier_scale = largest_magnitude(tangency.bound_multipliers);
    EXPECT_NEAR(tangency.bound_multipliers.back(), 0.0, 1e-12 * multiplier_scale);
    Result capped = solve_checked(ratio_spec(Objective::MaxSharpe, 0.0, 0.3), data);
    expect_cash_at_zero(data, capped);
    EXPECT_NE(capped.bound_multipliers.back(), 0.0);
    Result excluded = po::solve(ratio_spec(Objective::MaxSharpe, 0.05, 0.3), data).value();
    expect_no_point(excluded, SolveStatus::Infeasible);
}

TEST(MaxSharpe, RiskFreeAboveAllMeansLongOnly)
{
    Data data = random_data(6, 8);
    Spec long_only = ratio_spec(Objective::MaxSharpe, 0.0, 1.0);
    std::vector<double> mu(data.mu().begin(), data.mu().end());
    Data negative = Data::create(6, po::testing::negated(mu), data.covariance()).value();
    expect_no_point(po::solve(long_only, negative).value(), SolveStatus::Infeasible);
    data.set_risk_free(kRiskFreeAboveMeans);
    expect_no_point(po::solve(long_only, data).value(), SolveStatus::Infeasible);
    long_only.leverage = 2.0;
    expect_no_point(po::solve(long_only, data).value(), SolveStatus::Infeasible);
}

// Without bounds and 1' Sigma^-1 mu~ < 0, the tangency direction has a negative sum: along the efficient
// branch the ratio increases toward its supremum, the slope of the frontier's asymptote, reached only as
// 1'y -> 0.
TEST(MaxSharpe, UnboundedWithShorts)
{
    Data data = random_data(6, 9);
    data.set_risk_free(kRiskFreeAboveMeans);
    ASSERT_LT(sum(covariance_solve(data, po::testing::excess_mu(data))), 0.0);
    Result result = po::solve(unbounded_spec(Objective::MaxSharpe), data).value();
    ASSERT_EQ(result.weights.size(), 7u);
    expect_no_point(result, SolveStatus::Unbounded);
    Data negative = Data::create(6, po::testing::excess_mu(data), data.covariance()).value();
    Result without_cash = po::solve(unbounded_spec(Objective::MaxSharpe), negative).value();
    expect_no_point(without_cash, SolveStatus::Unbounded);
}

// Long-only, a leverage with shorts, and both with cash: every case OPTIMAL with KKT, the bounds and the
// leverage binding in some.
TEST(MaxSharpe, LongOnlyAndLeverageKkt)
{
    std::size_t active_bounds = 0;
    std::size_t binding_leverage = 0;
    for (std::uint64_t seed = 1; seed <= 10; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = random_data(8, seed);
        if (seed % 2 == 0)
        {
            data.set_risk_free(kRiskFree);
        }
        Result long_only = solve_checked(ratio_spec(Objective::MaxSharpe, 0.0, 1.0), data);
        Spec leveraged = ratio_spec(Objective::MaxSharpe, -0.5, 1.0);
        leveraged.leverage = 1.4;
        Result with_leverage = solve_checked(leveraged, data);
        for (std::size_t i = 0; i < data.asset_count(); ++i)
        {
            active_bounds += long_only.bound_multipliers[i] != 0.0 ? 1 : 0;
        }
        EXPECT_LE(with_leverage.leverage_multiplier, 0.0);
        binding_leverage += with_leverage.leverage_multiplier < 0.0 ? 1 : 0;
    }
    EXPECT_GE(active_bounds, 10u);
    EXPECT_GE(binding_leverage, 5u);
}

// A random point of the long-only box with upper 0.4: uniform weights normalized, drawn again until they fit.
std::vector<double> random_feasible(std::size_t n, std::mt19937_64& engine)
{
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::vector<double> weights(n);
    for (int draw = 0; draw < 1000; ++draw)
    {
        for (double& weight : weights)
        {
            weight = uniform(engine);
        }
        double total = sum(weights);
        bool fits = true;
        for (double& weight : weights)
        {
            weight /= total;
            fits = fits && weight <= 0.4;
        }
        if (fits)
        {
            return weights;
        }
    }
    ADD_FAILURE() << "no feasible draw";
    return weights;
}

// 10 000 feasible portfolios, the convex combinations (1 - t) w* + t u of the optimum with random feasible u
// at t log-uniform in [1e-6, 1], so that many lie close to the optimum: none has a larger Sharpe ratio.
TEST(MaxSharpe, NoSharperFeasiblePortfolio)
{
    constexpr std::size_t kAssets = 6;
    Data data = random_data(kAssets, 21);
    Result result = solve_checked(ratio_spec(Objective::MaxSharpe, 0.0, 0.4), data);
    double best = ratio_of(data, data.mu(), result.weights);
    std::mt19937_64 engine(2026);
    std::uniform_real_distribution<double> exponent(-6.0, 0.0);
    for (int sample = 0; sample < 10000; ++sample)
    {
        std::vector<double> point = random_feasible(kAssets, engine);
        double t = std::pow(10.0, exponent(engine));
        for (std::size_t i = 0; i < kAssets; ++i)
        {
            point[i] = (1.0 - t) * result.weights[i] + t * point[i];
        }
        double ratio = ratio_of(data, data.mu(), point);
        ASSERT_LE(ratio, best * (1.0 + 1e-12)) << "sample " << sample << " t " << t;
    }
}

TEST(MaxDiversification, UnconstrainedClosedForm)
{
    for (std::uint64_t seed = 1; seed <= 5; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = random_data(6, seed);
        Result result = solve_checked(unbounded_spec(Objective::MaxDiversification), data);
        expect_closed_form(data, result, ratio_vector(data, Objective::MaxDiversification));
    }
}

TEST(MaxDiversification, LongOnlyKkt)
{
    std::size_t active_bounds = 0;
    for (std::uint64_t seed = 1; seed <= 10; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = Data::create(8, {}, po::testing::random_covariance(8, seed, 1e4)).value();
        Result result = solve_checked(ratio_spec(Objective::MaxDiversification, 0.0, 0.3), data);
        for (double beta : result.bound_multipliers)
        {
            active_bounds += beta != 0.0 ? 1 : 0;
        }
    }
    EXPECT_GE(active_bounds, 10u);
}

// The ratio does not depend on cash: with cash fixed at 0 the risky weights are those without cash, bit for
// bit, the budget multiplier becomes beta_0 = -lambda_budget and mu is not needed.
TEST(MaxDiversification, CashWeightZero)
{
    Data data = Data::create(6, {}, po::testing::random_covariance(6, 4, 1e4)).value();
    Spec spec = ratio_spec(Objective::MaxDiversification, 0.0, 0.3);
    Result without_cash = solve_checked(spec, data);
    data.set_risk_free(kRiskFree);
    Result with_cash = solve_checked(spec, data);
    expect_cash_at_zero(data, with_cash);
    for (std::size_t i = 0; i < data.asset_count(); ++i)
    {
        EXPECT_EQ(with_cash.weights[i], without_cash.weights[i]);
        EXPECT_EQ(with_cash.bound_multipliers[i], without_cash.bound_multipliers[i]);
    }
    EXPECT_EQ(with_cash.bound_multipliers.back(), -without_cash.budget_multiplier);
    EXPECT_NE(with_cash.bound_multipliers.back(), 0.0);
}

} // namespace
