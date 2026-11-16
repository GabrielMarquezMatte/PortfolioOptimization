#include "solve/requirements.hpp"

#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <expected>
#include <optional>
#include <string>
#include <vector>

#include "solve/solve.hpp"
#include "random_covariance.hpp"

namespace
{

using po::Data;
using po::Failure;
using po::FailureKind;
using po::Objective;
using po::Spec;

constexpr std::size_t kAssets = 3;
constexpr std::size_t kPeriods = 5;

Data with_returns(Data data)
{
    std::vector<double> returns = po::testing::random_vector(kPeriods * kAssets, 3, -0.05, 0.05);
    std::vector<double> benchmark = po::testing::random_vector(kPeriods, 4, -0.05, 0.05);
    EXPECT_TRUE(data.set_returns(returns, benchmark, kPeriods).has_value());
    return data;
}

Data full_data()
{
    std::vector<double> mu = {0.05, 0.08, 0.11};
    return with_returns(Data::create(kAssets, mu, po::testing::random_covariance(kAssets, 1)).value());
}

Data without_mu()
{
    return Data::create(kAssets, {}, po::testing::random_covariance(kAssets, 1)).value();
}

Data without_covariance()
{
    return Data::create(kAssets, std::vector<double>{0.05, 0.08, 0.11}, {}).value();
}

Spec spec_of(Objective objective)
{
    Spec spec;
    spec.objective = objective;
    return spec;
}

void expect_invalid(const Spec& spec, const Data& data, const std::string& message)
{
    std::optional<Failure> failure = po::problems::check_request(spec, data);
    ASSERT_TRUE(failure.has_value());
    EXPECT_EQ(failure->kind, FailureKind::InvalidArgument);
    EXPECT_EQ(failure->message, message);
}

void expect_accepted(const Spec& spec, const Data& data)
{
    std::optional<Failure> failure = po::problems::check_request(spec, data);
    EXPECT_FALSE(failure.has_value()) << failure.value_or(Failure{}).message;
}

TEST(Requirements, MinVarianceRejectsRiskAversion)
{
    Spec spec;
    spec.risk_aversion = 3.0;
    expect_invalid(spec, full_data(), "risk_aversion does not apply to MIN_VARIANCE");
}

TEST(Requirements, UtilityRejectsTargetReturn)
{
    Spec spec = spec_of(Objective::Utility);
    spec.target_return = 0.08;
    expect_invalid(spec, full_data(), "target_return does not apply to UTILITY");
}

TEST(Requirements, UtilityRequiresMu)
{
    expect_invalid(spec_of(Objective::Utility), without_mu(), "UTILITY requires mu");
}

TEST(Requirements, ImportanceRejectsMaxVolatility)
{
    Spec spec = spec_of(Objective::Importance);
    spec.max_volatility = 0.2;
    expect_invalid(spec, full_data(), "max_volatility does not apply to IMPORTANCE");
}

TEST(Requirements, KellyRejectsImportance)
{
    Spec spec = spec_of(Objective::Kelly);
    spec.importance = 0.5;
    expect_invalid(spec, full_data(), "importance does not apply to KELLY");
}

TEST(Requirements, KellyRequiresRateAboveMinusOne)
{
    Data data = full_data();
    data.set_risk_free(-1.0);
    expect_invalid(spec_of(Objective::Kelly), data, "KELLY requires risk_free > -1");
    data.set_risk_free(-0.5);
    expect_accepted(spec_of(Objective::Kelly), data);
}

TEST(Requirements, UnimplementedObjectivesStayUnsupported)
{
    Spec spec = spec_of(Objective::RiskParity);
    spec.risk_aversion = 3.0;
    std::optional<Failure> failure = po::problems::check_request(spec, Data::create(kAssets, {}, {}).value());
    ASSERT_TRUE(failure.has_value());
    EXPECT_EQ(failure->kind, FailureKind::Unsupported);
    EXPECT_EQ(failure->message, "objective not implemented yet");
}

TEST(Requirements, MaxSharpeRequiresMu)
{
    expect_invalid(spec_of(Objective::MaxSharpe), without_mu(), "MAX_SHARPE requires mu");
    expect_invalid(spec_of(Objective::MaxSharpe), without_covariance(), "MAX_SHARPE requires covariance");
    expect_accepted(spec_of(Objective::MaxSharpe), full_data());
}

TEST(Requirements, MaxDiversificationNeedsNoMu)
{
    expect_accepted(spec_of(Objective::MaxDiversification), without_mu());
    expect_invalid(spec_of(Objective::MaxDiversification), without_covariance(),
                   "MAX_DIVERSIFICATION requires covariance");
}

TEST(Requirements, MaxReturnRequiresMaxVolatility)
{
    Spec spec = spec_of(Objective::MaxReturn);
    expect_invalid(spec, full_data(), "MAX_RETURN requires max_volatility");
    spec.max_volatility = 0.2;
    expect_accepted(spec, full_data());
    expect_invalid(spec, without_mu(), "MAX_RETURN requires mu");
    expect_invalid(spec, without_covariance(), "MAX_RETURN requires covariance");
}

TEST(Requirements, MaxVolatilityOnlyForMaxReturn)
{
    constexpr std::array kOthers = {Objective::MinVariance, Objective::Utility, Objective::Importance,
                                    Objective::Kelly, Objective::TrackingError};
    constexpr std::array kNames = {"MIN_VARIANCE", "UTILITY", "IMPORTANCE", "KELLY", "TRACKING_ERROR"};
    for (std::size_t k = 0; k < kOthers.size(); ++k)
    {
        Spec spec = spec_of(kOthers[k]);
        spec.max_volatility = 0.2;
        expect_invalid(spec, full_data(), std::string("max_volatility does not apply to ") + kNames[k]);
    }
    Spec spec = spec_of(Objective::MaxReturn);
    spec.max_volatility = 0.2;
    expect_accepted(spec, full_data());
    spec.target_return = 0.08;
    expect_invalid(spec, full_data(), "target_return does not apply to MAX_RETURN");
}

TEST(Requirements, LeverageAcceptedByQpObjectives)
{
    constexpr std::array kAccepting = {Objective::MinVariance, Objective::Utility,      Objective::Importance,
                                       Objective::Kelly,       Objective::TrackingError, Objective::MaxSharpe,
                                       Objective::MaxDiversification};
    for (Objective objective : kAccepting)
    {
        SCOPED_TRACE(static_cast<int>(objective));
        Spec spec = spec_of(objective);
        spec.leverage = 1.5;
        expect_accepted(spec, full_data());
    }
    Spec spec = spec_of(Objective::MaxReturn);
    spec.max_volatility = 0.2;
    spec.leverage = 1.5;
    expect_accepted(spec, full_data());
}

TEST(Requirements, MinVarianceRequiresCovariance)
{
    expect_invalid(Spec{}, without_covariance(), "MIN_VARIANCE requires covariance");
}

TEST(Requirements, MinVarianceTargetReturnRequiresMu)
{
    Spec spec;
    expect_accepted(spec, without_mu());
    spec.target_return = 0.08;
    expect_invalid(spec, without_mu(), "MIN_VARIANCE requires mu");
}

TEST(Requirements, TrackingErrorRequiresReturns)
{
    Data data = Data::create(kAssets, std::vector<double>{0.05, 0.08, 0.11},
                             po::testing::random_covariance(kAssets, 1))
                    .value();
    expect_invalid(spec_of(Objective::TrackingError), data, "TRACKING_ERROR requires returns");
}

TEST(Requirements, TrackingErrorTargetReturnRequiresMu)
{
    Data data = with_returns(Data::create(kAssets, {}, {}).value());
    Spec spec = spec_of(Objective::TrackingError);
    expect_accepted(spec, data);
    spec.target_return = 0.08;
    expect_invalid(spec, data, "TRACKING_ERROR requires mu");
}

TEST(Requirements, MeanVarianceObjectivesRequireCovariance)
{
    expect_invalid(spec_of(Objective::Utility), without_covariance(), "UTILITY requires covariance");
    expect_invalid(spec_of(Objective::Importance), without_covariance(), "IMPORTANCE requires covariance");
    expect_invalid(spec_of(Objective::Kelly), without_covariance(), "KELLY requires covariance");
    expect_invalid(spec_of(Objective::Importance), without_mu(), "IMPORTANCE requires mu");
    expect_invalid(spec_of(Objective::Kelly), without_mu(), "KELLY requires mu");
}

struct Cell
{
    Objective objective;
    const char* objective_name;
    bool target_return;
    bool risk_aversion;
    bool importance;
    bool max_volatility;
};

Spec with_setter(Objective objective, int setter)
{
    Spec spec = spec_of(objective);
    switch (setter)
    {
    case 0:
        spec.target_return = 0.08;
        break;
    case 1:
        spec.risk_aversion = 3.0;
        break;
    case 2:
        spec.importance = 0.5;
        break;
    default:
        spec.max_volatility = 0.2;
        break;
    }
    return spec;
}

TEST(Requirements, EverySetterFollowsTheMatrix)
{
    constexpr std::array kCells = {
        Cell{Objective::MinVariance, "MIN_VARIANCE", true, false, false, false},
        Cell{Objective::MaxReturn, "MAX_RETURN", false, false, false, true},
        Cell{Objective::Utility, "UTILITY", false, true, false, false},
        Cell{Objective::Importance, "IMPORTANCE", false, false, true, false},
        Cell{Objective::Kelly, "KELLY", false, false, false, false},
        Cell{Objective::TrackingError, "TRACKING_ERROR", true, false, false, false},
        Cell{Objective::MaxSharpe, "MAX_SHARPE", false, false, false, false},
        Cell{Objective::MaxDiversification, "MAX_DIVERSIFICATION", false, false, false, false}};
    constexpr std::array kSetters = {"target_return", "risk_aversion", "importance", "max_volatility"};
    for (const Cell& cell : kCells)
    {
        std::array accepted = {cell.target_return, cell.risk_aversion, cell.importance, cell.max_volatility};
        for (int setter = 0; setter < 4; ++setter)
        {
            SCOPED_TRACE(std::string(cell.objective_name) + " " + kSetters[setter]);
            Spec spec = with_setter(cell.objective, setter);
            if (accepted[setter])
            {
                expect_accepted(spec, full_data());
                continue;
            }
            expect_invalid(spec, full_data(),
                           std::string(kSetters[setter]) + " does not apply to " + cell.objective_name);
        }
    }
}

TEST(Requirements, ImplementedCheckComesFirst)
{
    Spec spec = spec_of(Objective::RiskParity);
    spec.target_return = 0.08;
    std::optional<Failure> failure = po::problems::check_request(spec, Data::create(kAssets, {}, {}).value());
    ASSERT_TRUE(failure.has_value());
    EXPECT_EQ(failure->kind, FailureKind::Unsupported);
}

TEST(Requirements, SetterCheckComesBeforeData)
{
    Spec spec = spec_of(Objective::Utility);
    spec.importance = 0.5;
    expect_invalid(spec, Data::create(kAssets, {}, {}).value(), "importance does not apply to UTILITY");
}

TEST(Requirements, DataCheckComesBeforeBoundSizes)
{
    Spec spec = spec_of(Objective::Utility);
    spec.asset_lower.assign(kAssets + 1, 0.0);
    expect_invalid(spec, without_mu(), "UTILITY requires mu");
    expect_invalid(spec, full_data(), "asset_lower has 4 entries; expected 3");
}

TEST(Requirements, SolveReportsTheFailure)
{
    std::expected<po::Result, Failure> result = po::solve(spec_of(Objective::Utility), without_mu());
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error().kind, FailureKind::InvalidArgument);
    EXPECT_EQ(result.error().message, "UTILITY requires mu");
}

} // namespace
