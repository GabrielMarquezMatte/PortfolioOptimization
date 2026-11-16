#include "riskparity/risk_parity.hpp"

#include <cstddef>
#include <limits>
#include <span>
#include <utility>

#include "metrics/portfolio.hpp"
#include "riskparity/spinu.hpp"

namespace po::riskparity
{
namespace
{

double bound_of(const std::vector<double>& per_weight, double scalar, std::size_t weight)
{
    if (per_weight.empty())
    {
        return scalar;
    }
    return per_weight[weight];
}

bool inside_box(const Spec& spec, std::span<const double> weights)
{
    for (std::size_t i = 0; i < weights.size(); ++i)
    {
        bool above_lower = weights[i] >= bound_of(spec.asset_lower, spec.lower, i);
        bool below_upper = weights[i] <= bound_of(spec.asset_upper, spec.upper, i);
        if (!above_lower || !below_upper)
        {
            return false;
        }
    }
    return true;
}

Failure active_constraints_unsupported()
{
    return Failure{FailureKind::Unsupported, "RISK_PARITY with active constraints not implemented yet"};
}

Result no_point(std::size_t n, int iterations)
{
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    Result result;
    result.status = SolveStatus::MaxIterations;
    result.weights.assign(n, kNaN);
    result.bound_multipliers.assign(n, kNaN);
    result.budget_multiplier = kNaN;
    result.residual = kNaN;
    result.iterations = iterations;
    return result;
}

Result stationary_point(const Data& data, std::vector<double> weights, std::span<const double> budget,
                        int iterations)
{
    Result result;
    result.status = SolveStatus::Optimal;
    result.budget_multiplier = 0.0;
    result.bound_multipliers.assign(weights.size(), 0.0);
    result.residual = metrics::risk_budget_residual(data, weights, budget);
    result.iterations = iterations;
    result.weights = std::move(weights);
    return result;
}

} // namespace

std::vector<double> risk_budget_of(const Spec& spec, const Data& data)
{
    if (spec.risk_budget.empty())
    {
        return std::vector<double>(data.asset_count(), 1.0 / static_cast<double>(data.asset_count()));
    }
    return spec.risk_budget;
}

std::expected<Result, Failure> solve_risk_parity(const Spec& spec, const Data& data)
{
    if (spec.target_return.has_value() || spec.target_volatility.has_value())
    {
        return std::unexpected(active_constraints_unsupported());
    }
    std::vector<double> budget = risk_budget_of(spec, data);
    SpinuResult spinu_result = spinu(data, budget);
    if (!spinu_result.converged)
    {
        return no_point(data.asset_count(), spinu_result.iterations);
    }
    if (!inside_box(spec, spinu_result.weights))
    {
        return std::unexpected(active_constraints_unsupported());
    }
    return stationary_point(data, std::move(spinu_result.weights), budget, spinu_result.iterations);
}

} // namespace po::riskparity
