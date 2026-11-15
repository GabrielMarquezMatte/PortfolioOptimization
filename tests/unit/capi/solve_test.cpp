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
    SpecPtr spec = po::testing::make_spec(PO_OBJECTIVE_MAX_SHARPE);
    ResultPtr existing = po::testing::solve_ok(po::testing::make_spec().get(), data.get());
    po_result* raw = existing.get();
    EXPECT_EQ(po_solve(spec.get(), data.get(), &raw), PO_UNSUPPORTED);
    EXPECT_EQ(raw, nullptr);
    EXPECT_STREQ(po_last_error(), "objective not implemented yet");
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
