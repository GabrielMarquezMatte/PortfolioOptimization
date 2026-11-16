#include "riskparity/risk_parity.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <vector>

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "fuzz/risk_parity/risk_parity_reference.hpp"
#include "random_covariance.hpp"
#include "solve/solve.hpp"

namespace
{

using po::Data;
using po::Failure;
using po::FailureKind;
using po::Objective;
using po::Result;
using po::SolveStatus;
using po::Spec;
using po::testing::random_covariance;
using po::testing::random_vector;

Data sample_data(std::size_t n, std::uint64_t seed)
{
    return Data::create(n, random_vector(n, seed + 500, 0.03, 0.15), random_covariance(n, seed, 1e4)).value();
}

Spec risk_parity_spec()
{
    Spec spec;
    spec.objective = Objective::RiskParity;
    return spec;
}

Result solved(const Spec& spec, const Data& data)
{
    std::expected<Result, Failure> result = po::solve(spec, data);
    EXPECT_TRUE(result.has_value());
    return result.value_or(Result{});
}

TEST(RiskParity, BudgetOnlyIsSpinu)
{
    Data data = sample_data(6, 3);
    Spec spec = risk_parity_spec();
    spec.risk_budget = {0.3, 0.2, 0.1, 0.1, 0.2, 0.1};
    Result result = solved(spec, data);
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    EXPECT_EQ(result.budget_multiplier, 0.0);
    EXPECT_EQ(result.bound_multipliers, std::vector<double>(6, 0.0));
    EXPECT_TRUE(std::isnan(result.target_return_multiplier));
    EXPECT_TRUE(std::isnan(result.volatility_multiplier));
    EXPECT_TRUE(std::isnan(result.leverage_multiplier));
    EXPECT_LE(result.residual, 1e-18);
    po::testing::expect_risk_parity_kkt(spec, data, result, 1e-10);
}

TEST(RiskParity, DefaultBudgetIsOneOverN)
{
    Data data = sample_data(4, 5);
    Spec spec = risk_parity_spec();
    EXPECT_EQ(po::riskparity::risk_budget_of(spec, data), std::vector<double>(4, 0.25));
    spec.risk_budget = {0.4, 0.3, 0.2, 0.1};
    EXPECT_EQ(po::riskparity::risk_budget_of(spec, data), spec.risk_budget);
}

TEST(RiskParity, BudgetOnlyInsideTheBoxIsSpinu)
{
    Data data = sample_data(5, 7);
    Spec spec = risk_parity_spec();
    spec.lower = 0.0;
    spec.upper = 1.0;
    Result result = solved(spec, data);
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    EXPECT_EQ(result.weights, solved(risk_parity_spec(), data).weights);
    po::testing::expect_risk_parity_kkt(spec, data, result, 1e-10);
}

TEST(RiskParity, NegativeLowerBoundKeepsPositiveOrthant)
{
    Data data = sample_data(5, 9);
    Spec bounded = risk_parity_spec();
    bounded.lower = -0.5;
    bounded.upper = 1.0;
    Result open = solved(risk_parity_spec(), data);
    Result result = solved(bounded, data);
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    EXPECT_EQ(result.weights, open.weights);
    for (double weight : result.weights)
    {
        EXPECT_GT(weight, 0.0);
    }
}

TEST(RiskParity, SingleAssetIsTheWholeBudget)
{
    Data data = Data::create(1, {}, std::vector<double>{0.04}).value();
    Result result = solved(risk_parity_spec(), data);
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    ASSERT_EQ(result.weights.size(), 1u);
    EXPECT_NEAR(result.weights[0], 1.0, 1e-15);
}

TEST(RiskParity, ActiveConstraintsUnsupportedForNow)
{
    Data data = sample_data(4, 11);
    Spec target = risk_parity_spec();
    target.target_return = 0.08;
    std::expected<Result, Failure> result = po::solve(target, data);
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error().kind, FailureKind::Unsupported);
    EXPECT_EQ(result.error().message, "RISK_PARITY with active constraints not implemented yet");
    Spec outside = risk_parity_spec();
    outside.upper = 0.01;
    result = po::solve(outside, data);
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error().kind, FailureKind::Unsupported);
}

} // namespace
