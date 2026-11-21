// The C++ wrapper against the shared library: it only translates, so a solve through it is the C solve bit for bit,
// and every misuse it can see (a span of the wrong size, a moved-from handle) is an error, never a read out of bounds.
#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <portopt/portopt.hpp>
#include <type_traits>
#include <utility>
#include <vector>

namespace
{

static_assert(!std::is_copy_constructible_v<po::Data>);
static_assert(std::is_nothrow_move_constructible_v<po::Data>);

constexpr std::size_t kAssets = 3;
const std::vector<double> kMu{0.08, 0.12, 0.10};
const std::vector<double> kCov{0.040, 0.006, 0.010, 0.006, 0.090, 0.012, 0.010, 0.012, 0.060};

po::Data sample_data()
{
    po::Expected<po::Data> data = po::Data::create(kAssets, kMu, kCov);
    EXPECT_TRUE(data.has_value());
    return std::move(*data);
}

struct NativeSolve
{
    std::vector<double> weights;
    std::array<double, 4> multipliers;
};

// The same problem through the C interface alone.
NativeSolve solve_through_c()
{
    po_data* data = nullptr;
    po_spec* spec = nullptr;
    po_result* result = nullptr;
    EXPECT_EQ(po_data_create(kAssets, kMu.data(), kCov.data(), &data), PO_OK);
    EXPECT_EQ(po_spec_create(PO_OBJECTIVE_MIN_VARIANCE, &spec), PO_OK);
    EXPECT_EQ(po_spec_set_bounds(spec, 0.0, 1.0), PO_OK);
    EXPECT_EQ(po_solve(spec, data, &result), PO_OK);
    NativeSolve native{std::vector<double>(kAssets), {}};
    EXPECT_EQ(po_result_weights(result, native.weights.data(), kAssets), PO_OK);
    for (std::int32_t constraint = 0; constraint < 4; ++constraint)
    {
        EXPECT_EQ(po_result_multiplier(result, constraint, &native.multipliers.at(constraint)), PO_OK);
    }
    po_result_free(result);
    po_spec_free(spec);
    po_data_free(data);
    return native;
}

TEST(Wrapper, MinVarianceMatchesTheCApiBitForBit)
{
    const po::Data data = sample_data();
    const po::Expected<po::Result> result =
        po::solve(po::Spec{.objective = po::Objective::min_variance, .bounds = po::Bounds{0.0, 1.0}}, data);
    ASSERT_TRUE(result.has_value()) << result.error().message;
    EXPECT_EQ(result->status, po::SolveStatus::optimal);
    const NativeSolve native = solve_through_c();
    ASSERT_EQ(result->weights.size(), kAssets);
    EXPECT_EQ(std::memcmp(result->weights.data(), native.weights.data(), kAssets * sizeof(double)), 0);
    const std::array<double, 4> multipliers{result->budget_multiplier, result->target_return_multiplier,
        result->volatility_multiplier, result->leverage_multiplier};
    EXPECT_EQ(std::memcmp(multipliers.data(), native.multipliers.data(), sizeof(multipliers)), 0);
}

TEST(Wrapper, DefaultSpecIsMinVariance)
{
    const po::Spec spec;
    EXPECT_EQ(spec.objective, po::Objective::min_variance);
}

TEST(Wrapper, ErrorCarriesStatusAndMessage)
{
    const po::Expected<po::Result> result =
        po::solve(po::Spec{.objective = po::Objective::min_variance, .risk_aversion = 2.0}, sample_data());
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error().status, po::Status::invalid_argument);
    EXPECT_EQ(result.error().message, "risk_aversion does not apply to MIN_VARIANCE");
}

TEST(Wrapper, SpanSizeMismatchIsInvalidArgument)
{
    const std::vector<double> short_mu{0.08, 0.12};
    const po::Expected<po::Data> created = po::Data::create(kAssets, short_mu, kCov);
    ASSERT_FALSE(created.has_value());
    EXPECT_EQ(created.error().status, po::Status::invalid_argument);
    EXPECT_EQ(created.error().message, "mu has 2 entries, expected 3");

    const po::Expected<po::Data> small_cov = po::Data::create(kAssets, kMu, std::vector<double>(4, 0.1));
    ASSERT_FALSE(small_cov.has_value());
    EXPECT_EQ(small_cov.error().message, "cov has 4 entries, expected 3 x 3 = 9");

    po::Data data = sample_data();
    const po::Expected<double> volatility = data.volatility(std::vector<double>{0.5, 0.5});
    ASSERT_FALSE(volatility.has_value());
    EXPECT_EQ(volatility.error().status, po::Status::invalid_argument);

    const std::vector<double> returns(7, 0.01);
    const po::Expected<void> series = data.set_returns(returns, std::vector<double>(2, 0.01));
    ASSERT_FALSE(series.has_value());
    EXPECT_EQ(series.error().message, "returns has 7 entries, expected 2 x 3 = 6");

    const po::Expected<po::Result> bounds = po::solve(
        po::Spec{.objective = po::Objective::min_variance,
                 .bounds = po::WeightBounds{{0.0, 0.0, 0.0}, {1.0, 1.0}}},
        data);
    ASSERT_FALSE(bounds.has_value());
    EXPECT_EQ(bounds.error().message, "the lower bounds have 3 entries and the upper bounds 2");
}

TEST(Wrapper, MovedFromDataIsInvalidHandle)
{
    po::Data first = sample_data();
    const po::Data second = std::move(first);
    const std::vector<double> weights{0.2, 0.3, 0.5};
    // NOLINTNEXTLINE(bugprone-use-after-move): the moved-from handle is what this test is about.
    const po::Expected<double> moved = first.volatility(weights);
    ASSERT_FALSE(moved.has_value());
    EXPECT_EQ(moved.error().status, po::Status::invalid_handle);
    const po::Expected<double> kept = second.volatility(weights);
    ASSERT_TRUE(kept.has_value()) << kept.error().message;
    EXPECT_GT(*kept, 0.0);
}

TEST(Wrapper, InfeasibleIsAResultNotAnError)
{
    const po::Expected<po::Result> result = po::solve(
        po::Spec{.objective = po::Objective::min_variance, .bounds = po::Bounds{0.0, 1.0}, .target_return = 0.5},
        sample_data());
    ASSERT_TRUE(result.has_value()) << result.error().message;
    EXPECT_EQ(result->status, po::SolveStatus::infeasible);
}

TEST(Wrapper, MetricsAndRiskContributions)
{
    po::Data data = sample_data();
    ASSERT_TRUE(data.set_risk_free(0.02).has_value());
    EXPECT_EQ(data.asset_count().value(), kAssets);
    EXPECT_EQ(data.weight_count().value(), kAssets + 1);
    const std::vector<double> weights{0.2, 0.3, 0.4, 0.1};
    const po::Expected<std::vector<double>> contributions = data.risk_contributions(weights);
    ASSERT_TRUE(contributions.has_value()) << contributions.error().message;
    EXPECT_NEAR((*contributions)[0] + (*contributions)[1] + (*contributions)[2], 1.0, 1e-15);
    EXPECT_EQ((*contributions)[3], 0.0);
    EXPECT_TRUE(data.sharpe_ratio(weights).has_value());
    EXPECT_TRUE(data.utility(weights, 3.0).has_value());
    EXPECT_FALSE(data.tracking_error(weights).has_value());
}

TEST(Wrapper, AbiMismatchIsUnsupported)
{
    const po::Expected<void> check = po::detail::check_abi_version(1, 0);
    ASSERT_FALSE(check.has_value());
    EXPECT_EQ(check.error().status, po::Status::unsupported);
    EXPECT_EQ(check.error().message, "the PortOpt library has ABI version 1, but these headers need version 0");
}

TEST(Wrapper, VersionIsTheLibraryVersion)
{
    EXPECT_EQ(po::version(), std::string_view(po_version()));
}

} // namespace
