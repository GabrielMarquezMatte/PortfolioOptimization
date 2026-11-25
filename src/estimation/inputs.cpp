#include "estimation/inputs.hpp"

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

} // namespace po::estimation
