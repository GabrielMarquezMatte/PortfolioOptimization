#include "problems/solve.hpp"

#include <optional>
#include <utility>

#include "problems/homogeneous.hpp"
#include "problems/max_return.hpp"
#include "problems/mean_variance.hpp"
#include "problems/requirements.hpp"
#include "problems/tracking_error.hpp"

namespace po
{
namespace
{

std::expected<Result, Failure> solve_request(const Spec& spec, const Data& data)
{
    switch (spec.objective)
    {
    case Objective::MinVariance:
    case Objective::Utility:
    case Objective::Importance:
    case Objective::Kelly:
        return problems::solve_mean_variance(spec, data);
    case Objective::MaxReturn:
        return problems::solve_max_return(spec, data);
    case Objective::TrackingError:
        return problems::solve_tracking_error(spec, data);
    case Objective::MaxSharpe:
        return problems::solve_max_sharpe(spec, data);
    case Objective::MaxDiversification:
        return problems::solve_max_diversification(spec, data);
    default:
        std::unreachable();
    }
}

} // namespace

std::expected<Result, Failure> solve(const Spec& spec, const Data& data)
{
    std::optional<Failure> failure = problems::check_request(spec, data);
    if (failure.has_value())
    {
        return std::unexpected(std::move(*failure));
    }
    return solve_request(spec, data);
}

} // namespace po
