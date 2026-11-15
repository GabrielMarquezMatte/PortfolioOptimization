#include "problems/solve.hpp"

#include <gtest/gtest.h>

#include <array>
#include <expected>
#include <vector>

#include "random_covariance.hpp"

namespace
{

using po::FailureKind;
using po::Objective;

void expect_failure(const std::expected<po::Result, po::Failure>& result, FailureKind kind, const char* message)
{
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error().kind, kind);
    EXPECT_EQ(result.error().message, message);
}

TEST(Solve, MissingCovarianceIsInvalid)
{
    po::Data data = po::Data::create(3, std::vector<double>{0.1, 0.2, 0.3}, {}).value();
    expect_failure(po::solve(po::Spec{}, data), FailureKind::InvalidArgument, "MIN_VARIANCE requires covariance");
}

TEST(Solve, TargetReturnWithoutMuIsInvalid)
{
    po::Data data = po::Data::create(3, {}, po::testing::random_covariance(3, 1)).value();
    po::Spec spec;
    spec.target_return = 0.1;
    expect_failure(po::solve(spec, data), FailureKind::InvalidArgument, "MIN_VARIANCE requires mu");
}

TEST(Solve, UnimplementedObjectiveIsUnsupported)
{
    constexpr std::array kObjectives = {Objective::MaxSharpe, Objective::MaxDiversification,
                                        Objective::RiskParity};
    po::Data data = po::Data::create(3, {}, {}).value();
    for (Objective objective : kObjectives)
    {
        SCOPED_TRACE(static_cast<int>(objective));
        po::Spec spec;
        spec.objective = objective;
        expect_failure(po::solve(spec, data), FailureKind::Unsupported, "objective not implemented yet");
    }
}

} // namespace
