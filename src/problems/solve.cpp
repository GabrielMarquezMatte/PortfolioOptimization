#include "problems/solve.hpp"

#include <optional>
#include <string>
#include <utility>

#include "problems/min_variance.hpp"

namespace po
{
namespace
{

Failure invalid(std::string message)
{
    return Failure{FailureKind::InvalidArgument, std::move(message)};
}

std::optional<Failure> check_implemented(const Spec& spec)
{
    if (spec.objective != Objective::MinVariance)
    {
        return Failure{FailureKind::Unsupported, "objective not implemented yet"};
    }
    return std::nullopt;
}

std::optional<Failure> check_required_data(const Spec& spec, const Data& data)
{
    if (!data.has_covariance())
    {
        return invalid("MIN_VARIANCE requires a covariance");
    }
    if (spec.target_return.has_value() && !data.has_mu())
    {
        return invalid("target_return requires mu");
    }
    return std::nullopt;
}

std::optional<Failure> check_bound_size(const char* name, std::size_t size, std::size_t expected)
{
    if (size == 0 || size == expected)
    {
        return std::nullopt;
    }
    return invalid(std::string(name) + " has " + std::to_string(size) + " entries; expected " +
                   std::to_string(expected));
}

std::optional<Failure> check_bound_sizes(const Spec& spec, const Data& data)
{
    std::size_t expected = data.weight_count();
    return check_bound_size("asset_lower", spec.asset_lower.size(), expected).or_else([&] {
        return check_bound_size("asset_upper", spec.asset_upper.size(), expected);
    });
}

} // namespace

std::expected<Result, Failure> solve(const Spec& spec, const Data& data)
{
    std::optional<Failure> failure = check_implemented(spec)
                                         .or_else([&] { return check_required_data(spec, data); })
                                         .or_else([&] { return check_bound_sizes(spec, data); });
    if (failure.has_value())
    {
        return std::unexpected(std::move(*failure));
    }
    return problems::solve_min_variance(spec, data);
}

} // namespace po
