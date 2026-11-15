#include "problems/solve.hpp"

#include <optional>
#include <utility>

#include "problems/mean_variance.hpp"
#include "problems/requirements.hpp"

namespace po
{
namespace
{

Result solve_request(const Spec& spec, const Data& data)
{
    switch (spec.objective)
    {
    case Objective::MinVariance:
    case Objective::Utility:
    case Objective::Importance:
    case Objective::Kelly:
        return problems::solve_mean_variance(spec, data);
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
