#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <portopt/portopt.h>
#include <string>
#include <thread>
#include <vector>

#include "capi_scope.hpp"
#include "problems/solve.hpp"
#include "random_covariance.hpp"

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
    EXPECT_TRUE(std::isnan(multiplier_of(handle, PO_CONSTRAINT_LEVERAGE)));
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

TEST(CapiSolve, UnsupportedObjective)
{
    DataPtr data = sample_data();
    for (int32_t objective :
         {PO_OBJECTIVE_MAX_SHARPE, PO_OBJECTIVE_MAX_DIVERSIFICATION, PO_OBJECTIVE_RISK_PARITY})
    {
        SCOPED_TRACE(objective);
        SpecPtr spec = po::testing::make_spec(objective);
        ResultPtr existing = po::testing::solve_ok(po::testing::make_spec().get(), data.get());
        po_result* raw = existing.get();
        EXPECT_EQ(po_solve(spec.get(), data.get(), &raw), PO_UNSUPPORTED);
        EXPECT_EQ(raw, nullptr);
        EXPECT_STREQ(po_last_error(), "objective not implemented yet");
    }
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
    EXPECT_EQ(po_spec_set_asset_bounds(nullptr, values, values, kAssets), PO_INVALID_HANDLE);
    EXPECT_EQ(po_spec_set_target_return(nullptr, 0.1), PO_INVALID_HANDLE);
    EXPECT_EQ(po_spec_set_risk_aversion(nullptr, 2.0), PO_INVALID_HANDLE);
    EXPECT_EQ(po_spec_set_importance(nullptr, 0.5), PO_INVALID_HANDLE);
    EXPECT_EQ(po_spec_set_max_volatility(nullptr, 0.2), PO_INVALID_HANDLE);
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

} // namespace
