#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "kkt.hpp"
#include "problems/mean_variance.hpp"
#include "qp/goldfarb_idnani.hpp"
#include "random_covariance.hpp"
#include "solve/solve.hpp"
#include "user_kkt.hpp"

namespace
{

using po::Data;
using po::Objective;
using po::Result;
using po::SolveStatus;
using po::Spec;
using po::testing::covariance_times;
using po::testing::dot;
using po::testing::excess_mu;

constexpr double kRiskFree = 0.02;

// f = (quadratic / 2) w'Sigma w - linear (mu - rf 1)'w, written here independently of the builder.
struct Scales
{
    double quadratic;
    double linear;
};

constexpr Scales kMinVariance{2.0, 0.0};
constexpr double kDelta = 5.0;
constexpr Scales kUtility{kDelta, 1.0};

Data random_data(std::size_t n, std::uint64_t seed)
{
    std::vector<double> mu = po::testing::random_vector(n, seed + 500, 0.03, 0.15);
    return Data::create(n, mu, po::testing::random_covariance(n, seed)).value();
}

Spec with_leverage(Objective objective, double leverage)
{
    Spec spec;
    spec.objective = objective;
    spec.leverage = leverage;
    if (objective == Objective::Utility)
    {
        spec.risk_aversion = kDelta;
    }
    return spec;
}

double gross_exposure(const Data& data, const Result& result)
{
    double total = 0.0;
    for (double weight : std::span<const double>(result.weights).first(data.asset_count()))
    {
        total += std::fabs(weight);
    }
    return total;
}

po::testing::ObjectiveGradient objective_gradient(const Data& data, const Result& result, Scales scales)
{
    std::vector<double> curvature = covariance_times(data, result.weights);
    std::vector<double> excess = excess_mu(data);
    for (std::size_t i = 0; i < curvature.size(); ++i)
    {
        curvature[i] *= scales.quadratic;
        excess[i] *= -scales.linear;
    }
    return po::testing::sum_of({curvature, excess});
}

double objective_value(const Data& data, const Result& result, Scales scales)
{
    std::span<const double> risky = std::span<const double>(result.weights).first(data.asset_count());
    return 0.5 * scales.quadratic * dot(risky, covariance_times(data, risky)) -
           scales.linear * dot(risky, excess_mu(data));
}

// The QP of the builder satisfies KKT with the family, and the result satisfies KKT in the user's variables.
Result solve_checked(const Spec& spec, const Data& data, Scales scales)
{
    po::problems::Translation translation =
        po::problems::mean_variance_problem(spec, data, po::problems::mean_variance_of(spec, data));
    po::testing::expect_kkt(translation.problem, po::qp::solve(translation.problem));
    Result result = po::solve(spec, data).value();
    po::testing::expect_user_kkt(spec, data, result, objective_gradient(data, result, scales));
    return result;
}

// At L = 1.3 the leverage binds in 12 of the 20 cases; both the active and the inactive one are checked.
TEST(Leverage, MinVarianceAndUtilityKkt)
{
    std::size_t active = 0;
    for (std::uint64_t seed = 1; seed <= 10; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = random_data(8, seed);
        for (Objective objective : {Objective::MinVariance, Objective::Utility})
        {
            Scales scales = objective == Objective::Utility ? kUtility : kMinVariance;
            Result result = solve_checked(with_leverage(objective, 1.3), data, scales);
            EXPECT_LE(result.leverage_multiplier, 0.0);
            active += result.leverage_multiplier < 0.0 ? 1 : 0;
        }
    }
    EXPECT_GE(active, 5u);
    EXPECT_LE(active, 15u);
}

// w* = (0.9, 0.4, -0.3, 0) is the optimum of UTILITY with L = ||w*||_1 = 1.6 when mu is chosen so that
// delta Sigma w* - mu = budget 1 + lambda_lev g with g = (1, 1, -1, 0.5): g_4 = 0.5 is a subgradient strictly
// inside [-1, 1] at w_4 = 0, so the optimum keeps w_4 = 0 and no sign is right for it.
TEST(Leverage, ZeroWeightsSubgradient)
{
    constexpr std::size_t n = 4;
    constexpr double kBudget = 0.01;
    constexpr double kLambda = -0.02;
    const std::vector<double> optimum = {0.9, 0.4, -0.3, 0.0};
    const std::vector<double> subgradient = {1.0, 1.0, -1.0, 0.5};
    std::vector<double> covariance = po::testing::random_covariance(n, 31);
    Data unit = Data::create(n, {}, covariance).value();
    std::vector<double> mu = covariance_times(unit, optimum);
    for (std::size_t i = 0; i < n; ++i)
    {
        mu[i] = kDelta * mu[i] - kBudget - kLambda * subgradient[i];
    }
    Data data = Data::create(n, mu, covariance).value();
    Result result = solve_checked(with_leverage(Objective::Utility, 1.6), data, kUtility);
    for (std::size_t i = 0; i < n; ++i)
    {
        EXPECT_NEAR(result.weights[i], optimum[i], 1e-12);
    }
    EXPECT_NEAR(result.leverage_multiplier, kLambda, 1e-12);
    EXPECT_NEAR(result.budget_multiplier, kBudget, 1e-12);
}

TEST(Leverage, BelowOneWithoutCashIsInfeasible)
{
    Data data = random_data(6, 4);
    for (Objective objective : {Objective::MinVariance, Objective::Utility})
    {
        Result result = po::solve(with_leverage(objective, 0.9), data).value();
        EXPECT_EQ(result.status, SolveStatus::Infeasible);
        for (double weight : result.weights)
        {
            EXPECT_TRUE(std::isnan(weight));
        }
        EXPECT_TRUE(std::isnan(result.leverage_multiplier));
    }
}

// With cash, L bounds the risky weights only: at L = 0.2 the risky weights reach 0.2 and the cash weight, at least
// 0.8, is not counted.
TEST(Leverage, CashExcluded)
{
    constexpr double kLeverage = 0.2;
    for (std::uint64_t seed = 1; seed <= 5; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = random_data(8, seed);
        data.set_risk_free(kRiskFree);
        Result result = solve_checked(with_leverage(Objective::Utility, kLeverage), data, kUtility);
        EXPECT_LT(result.leverage_multiplier, 0.0);
        EXPECT_NEAR(gross_exposure(data, result), kLeverage, 1e-12);
        EXPECT_GE(result.weights.back(), 1.0 - kLeverage - 1e-12);
    }
}

// Seeds where L = 1.3 binds for both objectives.
TEST(Leverage, MultiplierIsDerivativeOfOptimum)
{
    constexpr double kLeverage = 1.3;
    constexpr double kStep = 1e-4 * kLeverage;
    for (std::uint64_t seed : {4, 6, 7, 9})
    {
        SCOPED_TRACE(seed);
        Data data = random_data(8, seed);
        for (Objective objective : {Objective::MinVariance, Objective::Utility})
        {
            Scales scales = objective == Objective::Utility ? kUtility : kMinVariance;
            Result result = po::solve(with_leverage(objective, kLeverage), data).value();
            ASSERT_EQ(result.status, SolveStatus::Optimal);
            EXPECT_LT(result.leverage_multiplier, 0.0);
            Result above = po::solve(with_leverage(objective, kLeverage + kStep), data).value();
            Result below = po::solve(with_leverage(objective, kLeverage - kStep), data).value();
            double difference =
                (objective_value(data, above, scales) - objective_value(data, below, scales)) / (2.0 * kStep);
            EXPECT_NEAR(difference, result.leverage_multiplier, 1e-6 * std::fabs(result.leverage_multiplier));
        }
    }
}

TEST(Leverage, MultiplierNaNWithoutLeverage)
{
    Data data = random_data(8, 2);
    Spec spec;
    Result without = po::solve(spec, data).value();
    EXPECT_TRUE(std::isnan(without.leverage_multiplier));
    spec.leverage = 1e6;
    Result inactive = po::solve(spec, data).value();
    EXPECT_EQ(inactive.leverage_multiplier, 0.0);
    EXPECT_FALSE(std::signbit(inactive.leverage_multiplier));
    EXPECT_EQ(inactive.weights, without.weights);
    EXPECT_EQ(inactive.bound_multipliers, without.bound_multipliers);
}

} // namespace
