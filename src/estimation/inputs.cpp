#include "estimation/inputs.hpp"

#include <cmath>
#include <string>

namespace po::estimation
{

std::optional<Failure> check_inputs(const Data& data, bool needs_mu)
{
    if (!data.has_covariance())
    {
        return Failure{FailureKind::InvalidArgument, "data has no covariance"};
    }
    if (needs_mu && !data.has_mu())
    {
        return Failure{FailureKind::InvalidArgument, "data has no mu"};
    }
    if (data.risk_free().has_value())
    {
        return Failure{FailureKind::InvalidArgument, "data must not have a risk-free rate"};
    }
    return std::nullopt;
}

std::optional<Failure> check_positive(double value, const char* message)
{
    if (std::isfinite(value) && value > 0.0)
    {
        return std::nullopt;
    }
    return Failure{FailureKind::InvalidArgument, message};
}

std::optional<Failure> check_views(const Data& data, const Views& views)
{
    if (views.asset_count() == data.asset_count())
    {
        return std::nullopt;
    }
    return Failure{FailureKind::InvalidArgument, "views have " + std::to_string(views.asset_count()) +
                                                     " assets; data has " + std::to_string(data.asset_count())};
}

} // namespace po::estimation
