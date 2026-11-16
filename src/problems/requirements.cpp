#include "problems/requirements.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <utility>

namespace po::problems
{
namespace
{

// One row of the objective x setter matrix. An objective requires returns when requires_returns is set and a
// covariance otherwise, and mu whenever target_return is set. max_volatility is accepted exactly where it is
// required.
struct Requirements
{
    const char* name;
    bool accepts_target_return;
    bool accepts_risk_aversion;
    bool accepts_importance;
    bool accepts_leverage;
    bool requires_max_volatility;
    bool requires_mu;
    bool requires_returns;
};

struct Setter
{
    const char* name;
    bool present;
    bool accepted;
};

Failure invalid(std::string message)
{
    return Failure{FailureKind::InvalidArgument, std::move(message)};
}

std::optional<Requirements> requirements_of(Objective objective)
{
    switch (objective)
    {
    case Objective::MinVariance:
        return Requirements{"MIN_VARIANCE", true, false, false, true, false, false, false};
    case Objective::MaxReturn:
        return Requirements{"MAX_RETURN", false, false, false, false, true, true, false};
    case Objective::Utility:
        return Requirements{"UTILITY", false, true, false, true, false, true, false};
    case Objective::Importance:
        return Requirements{"IMPORTANCE", false, false, true, true, false, true, false};
    case Objective::Kelly:
        return Requirements{"KELLY", false, false, false, true, false, true, false};
    case Objective::TrackingError:
        return Requirements{"TRACKING_ERROR", true, false, false, true, false, false, true};
    default:
        return std::nullopt;
    }
}

std::optional<Failure> check_setters(const Spec& spec, const Requirements& requirements)
{
    std::array setters = {
        Setter{"target_return", spec.target_return.has_value(), requirements.accepts_target_return},
        Setter{"risk_aversion", spec.risk_aversion.has_value(), requirements.accepts_risk_aversion},
        Setter{"importance", spec.importance.has_value(), requirements.accepts_importance},
        Setter{"max_volatility", spec.max_volatility.has_value(), requirements.requires_max_volatility},
        Setter{"leverage", spec.leverage.has_value(), requirements.accepts_leverage},
    };
    for (const Setter& setter : setters)
    {
        if (setter.present && !setter.accepted)
        {
            return invalid(std::string(setter.name) + " does not apply to " + requirements.name);
        }
    }
    return std::nullopt;
}

std::optional<Failure> check_matrix(const Data& data, const Requirements& requirements)
{
    if (requirements.requires_returns && !data.has_returns())
    {
        return invalid(std::string(requirements.name) + " requires returns");
    }
    if (!requirements.requires_returns && !data.has_covariance())
    {
        return invalid(std::string(requirements.name) + " requires covariance");
    }
    return std::nullopt;
}

std::optional<Failure> check_mu_and_rate(const Spec& spec, const Data& data, const Requirements& requirements)
{
    bool requires_mu = requirements.requires_mu || spec.target_return.has_value();
    if (requires_mu && !data.has_mu())
    {
        return invalid(std::string(requirements.name) + " requires mu");
    }
    if (spec.objective == Objective::Kelly && !(data.risk_free().value_or(0.0) > -1.0))
    {
        return invalid("KELLY requires risk_free > -1");
    }
    if (requirements.requires_max_volatility && !spec.max_volatility.has_value())
    {
        return invalid(std::string(requirements.name) + " requires max_volatility");
    }
    return std::nullopt;
}

std::optional<Failure> check_data(const Spec& spec, const Data& data, const Requirements& requirements)
{
    return check_matrix(data, requirements).or_else([&] {
        return check_mu_and_rate(spec, data, requirements);
    });
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

std::optional<Failure> check_request(const Spec& spec, const Data& data)
{
    std::optional<Requirements> requirements = requirements_of(spec.objective);
    if (!requirements.has_value())
    {
        return Failure{FailureKind::Unsupported, "objective not implemented yet"};
    }
    return check_setters(spec, *requirements)
        .or_else([&] { return check_data(spec, data, *requirements); })
        .or_else([&] { return check_bound_sizes(spec, data); });
}

} // namespace po::problems
