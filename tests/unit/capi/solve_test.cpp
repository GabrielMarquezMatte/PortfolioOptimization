#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <portopt/portopt.h>
#include <string>
#include <thread>
#include <vector>

#include "capi_scope.hpp"
#include "metrics/portfolio.hpp"
#include "problems/constraints.hpp"
#include "random_covariance.hpp"
#include "qp/problem.hpp"
#include "solve/solve.hpp"

namespace
{

using po::testing::DataPtr;
using po::testing::ResultPtr;
using po::testing::SpecPtr;

constexpr std::size_t kAssets = 4;

std::vector<double> sample_mu()
{
    return po::testing::random_vector(kAssets, 11, 0.03, 0.15);
}

std::vector<double> sample_covariance()
{
    return po::testing::random_covariance(kAssets, 5);
}

DataPtr sample_data()
{
    return po::testing::make_data(kAssets, sample_mu(), sample_covariance());
}

std::vector<double> bound_multipliers_of(const po_result* result, std::size_t count)
{
    std::vector<double> multipliers(count);
    EXPECT_EQ(po_result_bound_multipliers(result, multipliers.data(), count), PO_OK);
    return multipliers;
}

double multiplier_of(const po_result* result, int32_t constraint)
{
    double value = 0.0;
    EXPECT_EQ(po_result_multiplier(result, constraint, &value), PO_OK);
    return value;
}

TEST(CapiSolve, MinVarianceEndToEnd)
{
    DataPtr data = sample_data();
    SpecPtr spec = po::testing::make_spec();
    ASSERT_EQ(po_spec_set_bounds(spec.get(), 0.0, 1.0), PO_OK);
    ResultPtr result = po::testing::solve_ok(spec.get(), data.get());

    po::Spec core_spec;
    core_spec.lower = 0.0;
    core_spec.upper = 1.0;
    po::Data core_data = po::Data::create(kAssets, sample_mu(), sample_covariance()).value();
    po::Result expected = po::solve(core_spec, core_data).value();

    EXPECT_EQ(po::testing::status_of(result.get()), PO_SOLVE_OPTIMAL);
    EXPECT_EQ(po::testing::weights_of(result.get()), expected.weights);
    EXPECT_EQ(bound_multipliers_of(result.get(), kAssets), expected.bound_multipliers);
    EXPECT_EQ(multiplier_of(result.get(), PO_CONSTRAINT_BUDGET), expected.budget_multiplier);
    EXPECT_TRUE(std::isnan(multiplier_of(result.get(), PO_CONSTRAINT_TARGET_RETURN)));
    EXPECT_TRUE(std::isnan(multiplier_of(result.get(), PO_CONSTRAINT_VOLATILITY)));
    EXPECT_TRUE(std::isnan(multiplier_of(result.get(), PO_CONSTRAINT_LEVERAGE)));
    int32_t iterations = -1;
    ASSERT_EQ(po_result_iterations(result.get(), &iterations), PO_OK);
    EXPECT_EQ(iterations, expected.iterations);
    double residual = -1.0;
    ASSERT_EQ(po_result_residual(result.get(), &residual), PO_OK);
    EXPECT_EQ(residual, expected.residual);
}

constexpr std::size_t kPeriods = 10;
constexpr double kRiskFree = 0.01;

std::vector<double> sample_returns()
{
    return po::testing::random_vector(kPeriods * kAssets, 41, -0.05, 0.05);
}

std::vector<double> sample_benchmark()
{
    return po::testing::random_vector(kPeriods, 42, -0.03, 0.03);
}

struct DataPair
{
    DataPtr handle;
    po::Data core;
};

DataPair full_data(bool with_cash)
{
    DataPair pair{sample_data(), po::Data::create(kAssets, sample_mu(), sample_covariance()).value()};
    po::testing::set_returns_ok(pair.handle.get(), sample_returns(), sample_benchmark());
    EXPECT_TRUE(pair.core.set_returns(sample_returns(), sample_benchmark(), kPeriods).has_value());
    if (with_cash)
    {
        EXPECT_EQ(po_data_set_risk_free(pair.handle.get(), kRiskFree), PO_OK);
        pair.core.set_risk_free(kRiskFree);
    }
    return pair;
}

bool same_value(double left, double right)
{
    return (std::isnan(left) && std::isnan(right)) || left == right;
}

void expect_same_as_core(const po_spec* spec, const po::Spec& core_spec, const DataPair& pair)
{
    ResultPtr result = po::testing::solve_ok(spec, pair.handle.get());
    po::Result expected = po::solve(core_spec, pair.core).value();
    EXPECT_EQ(po::testing::status_of(result.get()), static_cast<int32_t>(expected.status));
    EXPECT_EQ(po::testing::weights_of(result.get()), expected.weights);
    const po_result* handle = result.get();
    EXPECT_EQ(bound_multipliers_of(handle, expected.bound_multipliers.size()), expected.bound_multipliers);
    EXPECT_TRUE(same_value(multiplier_of(handle, PO_CONSTRAINT_BUDGET), expected.budget_multiplier));
    double target = multiplier_of(handle, PO_CONSTRAINT_TARGET_RETURN);
    EXPECT_TRUE(same_value(target, expected.target_return_multiplier));
    EXPECT_TRUE(same_value(multiplier_of(handle, PO_CONSTRAINT_VOLATILITY), expected.volatility_multiplier));
    EXPECT_TRUE(same_value(multiplier_of(handle, PO_CONSTRAINT_LEVERAGE), expected.leverage_multiplier));
}

po::Spec core_spec_of(po::Objective objective)
{
    po::Spec spec;
    spec.objective = objective;
    spec.lower = 0.0;
    spec.upper = 1.0;
    return spec;
}

SpecPtr long_only_spec(int32_t objective)
{
    SpecPtr spec = po::testing::make_spec(objective);
    EXPECT_EQ(po_spec_set_bounds(spec.get(), 0.0, 1.0), PO_OK);
    return spec;
}

TEST(CapiSolve, EachQpObjectiveEndToEnd)
{
    for (bool with_cash : {false, true})
    {
        SCOPED_TRACE(with_cash);
        DataPair pair = full_data(with_cash);

        SpecPtr utility = long_only_spec(PO_OBJECTIVE_UTILITY);
        ASSERT_EQ(po_spec_set_risk_aversion(utility.get(), 3.0), PO_OK);
        po::Spec core_utility = core_spec_of(po::Objective::Utility);
        core_utility.risk_aversion = 3.0;
        expect_same_as_core(utility.get(), core_utility, pair);

        SpecPtr importance = long_only_spec(PO_OBJECTIVE_IMPORTANCE);
        ASSERT_EQ(po_spec_set_importance(importance.get(), 0.4), PO_OK);
        po::Spec core_importance = core_spec_of(po::Objective::Importance);
        core_importance.importance = 0.4;
        expect_same_as_core(importance.get(), core_importance, pair);

        SpecPtr kelly = long_only_spec(PO_OBJECTIVE_KELLY);
        expect_same_as_core(kelly.get(), core_spec_of(po::Objective::Kelly), pair);
        SpecPtr tracking = long_only_spec(PO_OBJECTIVE_TRACKING_ERROR);
        expect_same_as_core(tracking.get(), core_spec_of(po::Objective::TrackingError), pair);

        SpecPtr max_return = long_only_spec(PO_OBJECTIVE_MAX_RETURN);
        ASSERT_EQ(po_spec_set_max_volatility(max_return.get(), 0.3), PO_OK);
        po::Spec core_max_return = core_spec_of(po::Objective::MaxReturn);
        core_max_return.max_volatility = 0.3;
        expect_same_as_core(max_return.get(), core_max_return, pair);
    }
}

constexpr double kLeverage = 1.4;

TEST(CapiSolve, SharpeAndDiversificationEndToEnd)
{
    for (bool with_cash : {false, true})
    {
        DataPair pair = full_data(with_cash);
        for (po::Objective objective : {po::Objective::MaxSharpe, po::Objective::MaxDiversification})
        {
            for (bool with_leverage : {false, true})
            {
                SCOPED_TRACE(with_cash);
                SCOPED_TRACE(static_cast<int32_t>(objective));
                SCOPED_TRACE(with_leverage);
                SpecPtr spec = po::testing::make_spec(static_cast<int32_t>(objective));
                ASSERT_EQ(po_spec_set_bounds(spec.get(), -1.0, 1.0), PO_OK);
                po::Spec core_spec = core_spec_of(objective);
                core_spec.lower = -1.0;
                if (with_leverage)
                {
                    ASSERT_EQ(po_spec_set_leverage(spec.get(), kLeverage), PO_OK);
                    core_spec.leverage = kLeverage;
                }
                expect_same_as_core(spec.get(), core_spec, pair);
                ResultPtr result = po::testing::solve_ok(spec.get(), pair.handle.get());
                EXPECT_EQ(po::testing::status_of(result.get()), PO_SOLVE_OPTIMAL);
            }
        }
    }
}

DataPtr cash_above_every_mean()
{
    DataPtr data = po::testing::make_data(6, po::testing::random_vector(6, 509, 0.03, 0.15),
                                          po::testing::random_covariance(6, 9));
    EXPECT_EQ(po_data_set_risk_free(data.get(), 0.2), PO_OK);
    return data;
}

TEST(CapiSolve, UnboundedIsResultNotError)
{
    DataPtr data = cash_above_every_mean();
    SpecPtr spec = po::testing::make_spec(PO_OBJECTIVE_MAX_SHARPE);
    po_result* raw = nullptr;
    ASSERT_EQ(po_solve(spec.get(), data.get(), &raw), PO_OK);
    ResultPtr result(raw);
    EXPECT_EQ(po::testing::status_of(result.get()), PO_SOLVE_UNBOUNDED);
    std::vector<double> weights = po::testing::weights_of(result.get());
    ASSERT_EQ(weights.size(), 7u);
    for (double weight : weights)
    {
        EXPECT_TRUE(std::isnan(weight));
    }
    for (double multiplier : bound_multipliers_of(result.get(), 7))
    {
        EXPECT_TRUE(std::isnan(multiplier));
    }
    for (int32_t constraint : {PO_CONSTRAINT_BUDGET, PO_CONSTRAINT_TARGET_RETURN, PO_CONSTRAINT_VOLATILITY,
                               PO_CONSTRAINT_LEVERAGE})
    {
        EXPECT_TRUE(std::isnan(multiplier_of(result.get(), constraint)));
    }
}

double gross_of(const po_result* result)
{
    double gross = 0.0;
    for (double weight : po::testing::weights_of(result))
    {
        gross += std::abs(weight);
    }
    return gross;
}

// Correlations that make the unrestricted optimum of every objective short some assets.
DataPtr hedged_data()
{
    const double volatility[kAssets] = {0.1, 0.2, 0.3, 0.25};
    const double correlation[kAssets][kAssets] = {
        {1.0, 0.7, 0.1, 0.5}, {0.7, 1.0, -0.1, 0.2}, {0.1, -0.1, 1.0, -0.1}, {0.5, 0.2, -0.1, 1.0}};
    std::vector<double> covariance;
    for (std::size_t i = 0; i < kAssets; ++i)
    {
        for (std::size_t j = 0; j < kAssets; ++j)
        {
            covariance.push_back(volatility[i] * volatility[j] * correlation[i][j]);
        }
    }
    return po::testing::make_data(kAssets, sample_mu(), covariance);
}

SpecPtr boxed_spec(int32_t objective)
{
    SpecPtr spec = po::testing::make_spec(objective);
    EXPECT_EQ(po_spec_set_bounds(spec.get(), -3.0, 3.0), PO_OK);
    if (objective == PO_OBJECTIVE_MAX_RETURN)
    {
        EXPECT_EQ(po_spec_set_max_volatility(spec.get(), 0.3), PO_OK);
    }
    return spec;
}

TEST(CapiSolve, LeverageMultiplierOnlyWithLeverage)
{
    DataPtr data = hedged_data();
    for (int32_t objective : {PO_OBJECTIVE_MIN_VARIANCE, PO_OBJECTIVE_UTILITY, PO_OBJECTIVE_MAX_RETURN,
                              PO_OBJECTIVE_MAX_SHARPE, PO_OBJECTIVE_MAX_DIVERSIFICATION})
    {
        SCOPED_TRACE(objective);
        SpecPtr spec = boxed_spec(objective);
        ResultPtr unrestricted = po::testing::solve_ok(spec.get(), data.get());
        ASSERT_EQ(po::testing::status_of(unrestricted.get()), PO_SOLVE_OPTIMAL);
        EXPECT_TRUE(std::isnan(multiplier_of(unrestricted.get(), PO_CONSTRAINT_LEVERAGE)));
        double gross = gross_of(unrestricted.get());
        ASSERT_GT(gross, 1.0);

        ASSERT_EQ(po_spec_set_leverage(spec.get(), 2.0 * gross), PO_OK);
        ResultPtr inactive = po::testing::solve_ok(spec.get(), data.get());
        EXPECT_EQ(multiplier_of(inactive.get(), PO_CONSTRAINT_LEVERAGE), 0.0);

        ASSERT_EQ(po_spec_set_leverage(spec.get(), 0.5 * (1.0 + gross)), PO_OK);
        ResultPtr binding = po::testing::solve_ok(spec.get(), data.get());
        ASSERT_EQ(po::testing::status_of(binding.get()), PO_SOLVE_OPTIMAL);
        EXPECT_LT(multiplier_of(binding.get(), PO_CONSTRAINT_LEVERAGE), 0.0);
    }
}

TEST(CapiSolve, VolatilityMultiplierIsOnlyForMaxReturn)
{
    DataPair pair = full_data(false);
    for (int32_t objective :
         {PO_OBJECTIVE_UTILITY, PO_OBJECTIVE_IMPORTANCE, PO_OBJECTIVE_KELLY, PO_OBJECTIVE_TRACKING_ERROR})
    {
        SCOPED_TRACE(objective);
        ResultPtr result = po::testing::solve_ok(long_only_spec(objective).get(), pair.handle.get());
        EXPECT_TRUE(std::isnan(multiplier_of(result.get(), PO_CONSTRAINT_VOLATILITY)));
    }
    SpecPtr active = long_only_spec(PO_OBJECTIVE_MAX_RETURN);
    ASSERT_EQ(po_spec_set_max_volatility(active.get(), 0.3), PO_OK);
    ResultPtr capped = po::testing::solve_ok(active.get(), pair.handle.get());
    EXPECT_LT(multiplier_of(capped.get(), PO_CONSTRAINT_VOLATILITY), 0.0);
    SpecPtr loose = long_only_spec(PO_OBJECTIVE_MAX_RETURN);
    ASSERT_EQ(po_spec_set_max_volatility(loose.get(), 100.0), PO_OK);
    ResultPtr uncapped = po::testing::solve_ok(loose.get(), pair.handle.get());
    EXPECT_EQ(multiplier_of(uncapped.get(), PO_CONSTRAINT_VOLATILITY), 0.0);
}

TEST(CapiSolve, SetterOutsideMatrixIsInvalid)
{
    DataPair pair = full_data(false);
    po_result* raw = nullptr;
    SpecPtr utility = po::testing::make_spec(PO_OBJECTIVE_UTILITY);
    ASSERT_EQ(po_spec_set_target_return(utility.get(), 0.05), PO_OK);
    EXPECT_EQ(po_solve(utility.get(), pair.handle.get(), &raw), PO_INVALID_ARGUMENT);
    EXPECT_EQ(raw, nullptr);
    EXPECT_STREQ(po_last_error(), "target_return does not apply to UTILITY");

    SpecPtr importance = po::testing::make_spec(PO_OBJECTIVE_IMPORTANCE);
    ASSERT_EQ(po_spec_set_risk_aversion(importance.get(), 2.0), PO_OK);
    EXPECT_EQ(po_solve(importance.get(), pair.handle.get(), &raw), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "risk_aversion does not apply to IMPORTANCE");

    SpecPtr min_variance = po::testing::make_spec();
    ASSERT_EQ(po_spec_set_max_volatility(min_variance.get(), 0.2), PO_OK);
    EXPECT_EQ(po_solve(min_variance.get(), pair.handle.get(), &raw), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "max_volatility does not apply to MIN_VARIANCE");
}

TEST(CapiSolve, MissingRequirementsAreInvalid)
{
    po_result* raw = nullptr;
    DataPtr without_returns = sample_data();
    SpecPtr tracking = po::testing::make_spec(PO_OBJECTIVE_TRACKING_ERROR);
    EXPECT_EQ(po_solve(tracking.get(), without_returns.get(), &raw), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "TRACKING_ERROR requires returns");

    SpecPtr max_return = po::testing::make_spec(PO_OBJECTIVE_MAX_RETURN);
    EXPECT_EQ(po_solve(max_return.get(), without_returns.get(), &raw), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "MAX_RETURN requires max_volatility");

    DataPtr without_mu = po::testing::make_data(kAssets, {}, sample_covariance());
    SpecPtr utility = po::testing::make_spec(PO_OBJECTIVE_UTILITY);
    EXPECT_EQ(po_solve(utility.get(), without_mu.get(), &raw), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "UTILITY requires mu");

    SpecPtr kelly = po::testing::make_spec(PO_OBJECTIVE_KELLY);
    DataPtr ruinous = sample_data();
    ASSERT_EQ(po_data_set_risk_free(ruinous.get(), -1.5), PO_OK);
    EXPECT_EQ(po_solve(kelly.get(), ruinous.get(), &raw), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "KELLY requires risk_free > -1");
    EXPECT_EQ(raw, nullptr);
}

TEST(CapiSolve, TargetReturnHasMultiplier)
{
    DataPtr data = sample_data();
    SpecPtr spec = po::testing::make_spec();
    ASSERT_EQ(po_spec_set_target_return(spec.get(), 0.09), PO_OK);
    ResultPtr result = po::testing::solve_ok(spec.get(), data.get());
    EXPECT_EQ(po::testing::status_of(result.get()), PO_SOLVE_OPTIMAL);
    EXPECT_TRUE(std::isfinite(multiplier_of(result.get(), PO_CONSTRAINT_TARGET_RETURN)));
    for (double multiplier : bound_multipliers_of(result.get(), kAssets))
    {
        EXPECT_EQ(multiplier, 0.0);
    }
}

TEST(CapiSolve, InfeasibleIsResultNotError)
{
    DataPtr data = sample_data();
    SpecPtr spec = po::testing::make_spec();
    ASSERT_EQ(po_spec_set_bounds(spec.get(), 0.0, 0.1), PO_OK);
    po_result* raw = nullptr;
    ASSERT_EQ(po_solve(spec.get(), data.get(), &raw), PO_OK);
    ResultPtr result(raw);
    EXPECT_EQ(po::testing::status_of(result.get()), PO_SOLVE_INFEASIBLE);
    for (double weight : po::testing::weights_of(result.get()))
    {
        EXPECT_TRUE(std::isnan(weight));
    }
    double residual = 0.0;
    ASSERT_EQ(po_result_residual(result.get(), &residual), PO_OK);
    EXPECT_TRUE(std::isnan(residual));
}

TEST(CapiSolve, RiskParityIsSolved)
{
    DataPtr data = sample_data();
    SpecPtr spec = po::testing::make_spec(PO_OBJECTIVE_RISK_PARITY);
    ResultPtr result = po::testing::solve_ok(spec.get(), data.get());
    EXPECT_EQ(po::testing::status_of(result.get()), PO_SOLVE_OPTIMAL);
    double total = 0.0;
    for (double weight : po::testing::weights_of(result.get()))
    {
        total += weight;
    }
    EXPECT_NEAR(total, 1.0, 1e-12);
}

DataPair risk_parity_pair()
{
    return DataPair{sample_data(), po::Data::create(kAssets, sample_mu(), sample_covariance()).value()};
}

double residual_of(const po_result* result)
{
    double residual = -1.0;
    EXPECT_EQ(po_result_residual(result, &residual), PO_OK);
    return residual;
}

void expect_identical_to_core(const po_spec* spec, const po::Spec& core_spec, const DataPair& pair)
{
    expect_same_as_core(spec, core_spec, pair);
    ResultPtr result = po::testing::solve_ok(spec, pair.handle.get());
    po::Result expected = po::solve(core_spec, pair.core).value();
    EXPECT_EQ(po::testing::status_of(result.get()), PO_SOLVE_OPTIMAL);
    EXPECT_TRUE(same_value(residual_of(result.get()), expected.residual));
    int32_t iterations = -1;
    ASSERT_EQ(po_result_iterations(result.get(), &iterations), PO_OK);
    EXPECT_EQ(iterations, expected.iterations);
}

po::Spec core_risk_parity()
{
    po::Spec spec;
    spec.objective = po::Objective::RiskParity;
    return spec;
}

TEST(CapiSolve, RiskParityEndToEnd)
{
    DataPair pair = risk_parity_pair();
    constexpr double kTargetReturn = 0.09;

    SpecPtr budget_only = po::testing::make_spec(PO_OBJECTIVE_RISK_PARITY);
    expect_identical_to_core(budget_only.get(), core_risk_parity(), pair);

    SpecPtr with_return = po::testing::make_spec(PO_OBJECTIVE_RISK_PARITY);
    ASSERT_EQ(po_spec_set_bounds(with_return.get(), 0.0, 1.0), PO_OK);
    ASSERT_EQ(po_spec_set_target_return(with_return.get(), kTargetReturn), PO_OK);
    po::Spec core_return = core_risk_parity();
    core_return.lower = 0.0;
    core_return.upper = 1.0;
    core_return.target_return = kTargetReturn;
    expect_identical_to_core(with_return.get(), core_return, pair);

    ResultPtr spinu_result = po::testing::solve_ok(budget_only.get(), pair.handle.get());
    std::vector<double> spinu = po::testing::weights_of(spinu_result.get());
    double volatility = po::metrics::portfolio_volatility(pair.core, spinu);
    SpecPtr with_volatility = po::testing::make_spec(PO_OBJECTIVE_RISK_PARITY);
    ASSERT_EQ(po_spec_set_target_volatility(with_volatility.get(), 1.1 * volatility), PO_OK);
    po::Spec core_volatility = core_risk_parity();
    core_volatility.target_volatility = 1.1 * volatility;
    expect_identical_to_core(with_volatility.get(), core_volatility, pair);

    SpecPtr capped = po::testing::make_spec(PO_OBJECTIVE_RISK_PARITY);
    ASSERT_EQ(po_spec_set_bounds(capped.get(), 0.0, 0.3), PO_OK);
    const double budget[kAssets] = {0.4, 0.3, 0.2, 0.1};
    ASSERT_EQ(po_spec_set_risk_budget(capped.get(), budget, kAssets), PO_OK);
    po::Spec core_capped = core_risk_parity();
    core_capped.lower = 0.0;
    core_capped.upper = 0.3;
    core_capped.risk_budget.assign(budget, budget + kAssets);
    expect_identical_to_core(capped.get(), core_capped, pair);
}

TEST(CapiSolve, RiskParityRejectsCashAndLeverage)
{
    DataPair pair = risk_parity_pair();
    po_result* raw = nullptr;

    SpecPtr leverage = po::testing::make_spec(PO_OBJECTIVE_RISK_PARITY);
    ASSERT_EQ(po_spec_set_leverage(leverage.get(), 2.0), PO_OK);
    EXPECT_EQ(po_solve(leverage.get(), pair.handle.get(), &raw), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "leverage does not apply to RISK_PARITY");
    EXPECT_EQ(raw, nullptr);

    SpecPtr plain = po::testing::make_spec(PO_OBJECTIVE_RISK_PARITY);
    ASSERT_EQ(po_data_set_risk_free(pair.handle.get(), kRiskFree), PO_OK);
    EXPECT_EQ(po_solve(plain.get(), pair.handle.get(), &raw), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "risk_free does not apply to RISK_PARITY");
    EXPECT_EQ(raw, nullptr);
}

TEST(CapiSolve, RiskParitySetterMatrix)
{
    DataPair pair = risk_parity_pair();
    po_result* raw = nullptr;
    SpecPtr volatility = po::testing::make_spec(PO_OBJECTIVE_UTILITY);
    ASSERT_EQ(po_spec_set_target_volatility(volatility.get(), 0.2), PO_OK);
    EXPECT_EQ(po_solve(volatility.get(), pair.handle.get(), &raw), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "target_volatility does not apply to UTILITY");

    SpecPtr budget = po::testing::make_spec(PO_OBJECTIVE_RISK_PARITY);
    const double three[3] = {0.5, 0.3, 0.2};
    ASSERT_EQ(po_spec_set_risk_budget(budget.get(), three, 3), PO_OK);
    EXPECT_EQ(po_solve(budget.get(), pair.handle.get(), &raw), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "risk_budget has 3 entries; expected 4");

    SpecPtr other = po::testing::make_spec(PO_OBJECTIVE_MIN_VARIANCE);
    ASSERT_EQ(po_spec_set_risk_budget(other.get(), three, 3), PO_OK);
    EXPECT_EQ(po_solve(other.get(), pair.handle.get(), &raw), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "risk_budget does not apply to MIN_VARIANCE");

    SpecPtr target = po::testing::make_spec(PO_OBJECTIVE_RISK_PARITY);
    ASSERT_EQ(po_spec_set_target_return(target.get(), 0.09), PO_OK);
    DataPtr without_mu = po::testing::make_data(kAssets, {}, sample_covariance());
    EXPECT_EQ(po_solve(target.get(), without_mu.get(), &raw), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "RISK_PARITY requires mu");
}

TEST(CapiSolve, RiskParityResidualAndMultipliers)
{
    DataPair pair = risk_parity_pair();
    SpecPtr spec = po::testing::make_spec(PO_OBJECTIVE_RISK_PARITY);
    const double budget[kAssets] = {0.4, 0.3, 0.2, 0.1};
    ASSERT_EQ(po_spec_set_risk_budget(spec.get(), budget, kAssets), PO_OK);
    ResultPtr result = po::testing::solve_ok(spec.get(), pair.handle.get());
    ASSERT_EQ(po::testing::status_of(result.get()), PO_SOLVE_OPTIMAL);
    std::vector<double> weights = po::testing::weights_of(result.get());
    EXPECT_EQ(residual_of(result.get()), po::metrics::risk_budget_residual(pair.core, weights, budget));
    EXPECT_LE(residual_of(result.get()), 1e-18);
    EXPECT_TRUE(std::isfinite(multiplier_of(result.get(), PO_CONSTRAINT_BUDGET)));
    EXPECT_TRUE(std::isnan(multiplier_of(result.get(), PO_CONSTRAINT_TARGET_RETURN)));
    EXPECT_TRUE(std::isnan(multiplier_of(result.get(), PO_CONSTRAINT_VOLATILITY)));
    EXPECT_TRUE(std::isnan(multiplier_of(result.get(), PO_CONSTRAINT_LEVERAGE)));
    for (double multiplier : bound_multipliers_of(result.get(), kAssets))
    {
        EXPECT_EQ(multiplier, 0.0);
    }

    double volatility = po::metrics::portfolio_volatility(pair.core, weights);
    ASSERT_EQ(po_spec_set_target_volatility(spec.get(), 1.1 * volatility), PO_OK);
    ResultPtr targeted = po::testing::solve_ok(spec.get(), pair.handle.get());
    ASSERT_EQ(po::testing::status_of(targeted.get()), PO_SOLVE_OPTIMAL);
    EXPECT_TRUE(std::isfinite(multiplier_of(targeted.get(), PO_CONSTRAINT_VOLATILITY)));
    EXPECT_TRUE(std::isnan(multiplier_of(targeted.get(), PO_CONSTRAINT_TARGET_RETURN)));
    EXPECT_GT(residual_of(targeted.get()), 0.0);
    std::vector<double> moved = po::testing::weights_of(targeted.get());
    EXPECT_EQ(residual_of(targeted.get()), po::metrics::risk_budget_residual(pair.core, moved, budget));
    EXPECT_NEAR(po::metrics::portfolio_volatility(pair.core, moved), 1.1 * volatility, 1e-9 * volatility);
}

void expect_no_point(const po_result* result, std::size_t count)
{
    EXPECT_EQ(po::testing::status_of(result), PO_SOLVE_NOT_CONVERGED);
    for (double weight : po::testing::weights_of(result))
    {
        EXPECT_TRUE(std::isnan(weight));
    }
    for (double multiplier : bound_multipliers_of(result, count))
    {
        EXPECT_TRUE(std::isnan(multiplier));
    }
    for (int32_t constraint : {PO_CONSTRAINT_BUDGET, PO_CONSTRAINT_TARGET_RETURN, PO_CONSTRAINT_VOLATILITY,
                               PO_CONSTRAINT_LEVERAGE})
    {
        EXPECT_TRUE(std::isnan(multiplier_of(result, constraint)));
    }
    EXPECT_TRUE(std::isnan(residual_of(result)));
}

// Sigma = diag(0.04, 0.09), long only: no portfolio has a volatility above 0.3, so s = 0.5 is out of reach.
TEST(CapiSolve, NotConvergedHasNoPoint)
{
    const std::vector<double> covariance = {0.04, 0.0, 0.0, 0.09};
    DataPtr data = po::testing::make_data(2, {}, covariance);
    SpecPtr spec = po::testing::make_spec(PO_OBJECTIVE_RISK_PARITY);
    ASSERT_EQ(po_spec_set_bounds(spec.get(), 0.0, 1.0), PO_OK);
    ASSERT_EQ(po_spec_set_target_volatility(spec.get(), 0.5), PO_OK);
    ResultPtr result = po::testing::solve_ok(spec.get(), data.get());
    expect_no_point(result.get(), 2);
}

TEST(CapiSolve, ResultOutlivesSpecAndData)
{
    DataPtr data = sample_data();
    SpecPtr spec = po::testing::make_spec();
    ResultPtr result = po::testing::solve_ok(spec.get(), data.get());
    std::vector<double> expected = po::testing::weights_of(result.get());
    spec.reset();
    data.reset();
    EXPECT_EQ(po::testing::weights_of(result.get()), expected);
    EXPECT_EQ(po::testing::status_of(result.get()), PO_SOLVE_OPTIMAL);
}

TEST(CapiSolve, WrongVectorSizeIsInvalid)
{
    DataPtr data = sample_data();
    SpecPtr spec = po::testing::make_spec();
    ResultPtr result = po::testing::solve_ok(spec.get(), data.get());
    std::vector<double> buffer(kAssets + 1, 0.0);
    for (std::size_t n : {std::size_t{0}, kAssets - 1, kAssets + 1})
    {
        SCOPED_TRACE(n);
        EXPECT_EQ(po_result_weights(result.get(), buffer.data(), n), PO_INVALID_ARGUMENT);
        EXPECT_STREQ(po_last_error(), ("n must be 4, got " + std::to_string(n)).c_str());
        EXPECT_EQ(po_result_bound_multipliers(result.get(), buffer.data(), n), PO_INVALID_ARGUMENT);
    }
    EXPECT_EQ(buffer, std::vector<double>(kAssets + 1, 0.0));
}

TEST(CapiSolve, InvalidConstraintIsInvalid)
{
    DataPtr data = sample_data();
    SpecPtr spec = po::testing::make_spec();
    ResultPtr result = po::testing::solve_ok(spec.get(), data.get());
    double value = 7.0;
    for (int32_t constraint : {-1, 4})
    {
        EXPECT_EQ(po_result_multiplier(result.get(), constraint, &value), PO_INVALID_ARGUMENT);
    }
    EXPECT_STREQ(po_last_error(), "constraint is not one of PO_CONSTRAINT_*");
    EXPECT_EQ(value, 7.0);
}

TEST(CapiSolve, NullHandles)
{
    DataPtr data = sample_data();
    SpecPtr spec = po::testing::make_spec();
    const double values[kAssets] = {0.0, 0.0, 0.0, 0.0};
    double number = 0.0;
    double buffer[kAssets] = {};
    size_t count = 0;
    int32_t code = 0;
    po_result* out = nullptr;

    EXPECT_EQ(po_spec_set_bounds(nullptr, 0.0, 1.0), PO_INVALID_HANDLE);
    EXPECT_EQ(po_spec_set_weight_bounds(nullptr, values, values, kAssets), PO_INVALID_HANDLE);
    EXPECT_EQ(po_spec_set_target_return(nullptr, 0.1), PO_INVALID_HANDLE);
    EXPECT_EQ(po_spec_set_risk_aversion(nullptr, 2.0), PO_INVALID_HANDLE);
    EXPECT_EQ(po_spec_set_importance(nullptr, 0.5), PO_INVALID_HANDLE);
    EXPECT_EQ(po_spec_set_max_volatility(nullptr, 0.2), PO_INVALID_HANDLE);
    EXPECT_EQ(po_spec_set_leverage(nullptr, 2.0), PO_INVALID_HANDLE);
    EXPECT_EQ(po_data_set_returns(nullptr, values, values, 1), PO_INVALID_HANDLE);
    EXPECT_EQ(po_data_set_risk_free(nullptr, 0.1), PO_INVALID_HANDLE);
    EXPECT_EQ(po_data_asset_count(nullptr, &count), PO_INVALID_HANDLE);
    EXPECT_EQ(po_data_weight_count(nullptr, &count), PO_INVALID_HANDLE);
    EXPECT_EQ(po_solve(nullptr, data.get(), &out), PO_INVALID_HANDLE);
    EXPECT_EQ(po_solve(spec.get(), nullptr, &out), PO_INVALID_HANDLE);
    EXPECT_EQ(out, nullptr);
    EXPECT_EQ(po_result_status(nullptr, &code), PO_INVALID_HANDLE);
    EXPECT_EQ(po_result_weight_count(nullptr, &count), PO_INVALID_HANDLE);
    EXPECT_EQ(po_result_weights(nullptr, buffer, kAssets), PO_INVALID_HANDLE);
    EXPECT_EQ(po_result_iterations(nullptr, &code), PO_INVALID_HANDLE);
    EXPECT_EQ(po_result_multiplier(nullptr, PO_CONSTRAINT_BUDGET, &number), PO_INVALID_HANDLE);
    EXPECT_EQ(po_result_bound_multipliers(nullptr, buffer, kAssets), PO_INVALID_HANDLE);
    EXPECT_EQ(po_result_residual(nullptr, &number), PO_INVALID_HANDLE);
    EXPECT_STREQ(po_last_error(), "r must not be null");

    po_spec_free(nullptr);
    po_data_free(nullptr);
    po_result_free(nullptr);
}

TEST(CapiSolve, NullOutputPointers)
{
    DataPtr data = sample_data();
    SpecPtr spec = po::testing::make_spec();
    ResultPtr result = po::testing::solve_ok(spec.get(), data.get());

    EXPECT_EQ(po_data_create(kAssets, nullptr, nullptr, nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_data_asset_count(data.get(), nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_data_weight_count(data.get(), nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_solve(spec.get(), data.get(), nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_result_status(result.get(), nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_result_weight_count(result.get(), nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_result_weights(result.get(), nullptr, kAssets), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_result_iterations(result.get(), nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_result_multiplier(result.get(), PO_CONSTRAINT_BUDGET, nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_result_bound_multipliers(result.get(), nullptr, kAssets), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_result_residual(result.get(), nullptr), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "out must not be null");
}

TEST(CapiSolve, ParallelSolvesShareReadOnlyHandles)
{
    DataPtr data = sample_data();
    SpecPtr spec = po::testing::make_spec();
    ASSERT_EQ(po_spec_set_bounds(spec.get(), 0.0, 1.0), PO_OK);
    std::vector<double> expected = po::testing::weights_of(po::testing::solve_ok(spec.get(), data.get()).get());
    std::vector<std::vector<double>> weights(4);
    std::vector<std::thread> threads;
    for (std::vector<double>& slot : weights)
    {
        threads.emplace_back(
            [&]
            {
                po_result* raw = nullptr;
                po_solve(spec.get(), data.get(), &raw);
                ResultPtr result(raw);
                slot = po::testing::weights_of(result.get());
            });
    }
    for (std::thread& thread : threads)
    {
        thread.join();
    }
    for (const std::vector<double>& slot : weights)
    {
        EXPECT_EQ(slot, expected);
    }
}

// No QP objective reaches its iteration limit on a test-sized problem, so the rule is checked where it is made:
// a QP that stops at its limit becomes a result with NaN weights, multipliers and residual, for every shape
// of data (with cash too), which po_solve hands to the caller unchanged.
TEST(CapiSolve, NotConvergedHasNoPointForQpObjectives)
{
    constexpr int kIterations = 17;
    for (bool with_cash : {false, true})
    {
        SCOPED_TRACE(with_cash);
        DataPair pair = full_data(with_cash);
        po::qp::Solution solution;
        solution.status = po::qp::Status::MaxIterations;
        solution.iterations = kIterations;
        po::Result result = po::problems::to_result(pair.core, po::problems::ConstraintRows{}, solution);
        EXPECT_EQ(result.status, po::SolveStatus::NotConverged);
        EXPECT_EQ(result.iterations, kIterations);
        ASSERT_EQ(result.weights.size(), pair.core.weight_count());
        ASSERT_EQ(result.bound_multipliers.size(), pair.core.weight_count());
        for (double value : result.weights)
        {
            EXPECT_TRUE(std::isnan(value));
        }
        for (double value : result.bound_multipliers)
        {
            EXPECT_TRUE(std::isnan(value));
        }
        for (double value : {result.budget_multiplier, result.target_return_multiplier,
                             result.volatility_multiplier, result.leverage_multiplier, result.residual})
        {
            EXPECT_TRUE(std::isnan(value));
        }
    }
}

} // namespace
