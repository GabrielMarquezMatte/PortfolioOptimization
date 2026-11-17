#pragma once

#include <algorithm>
#include <cstddef>
#include <format>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/data.hpp"
#include "core/spec.hpp"
#include "instances.hpp"
#include "metrics/portfolio.hpp"
#include "solve/solve.hpp"

namespace po::bench
{

// What a case adds to the default spec. Free is no bound at all; the others are named after the C API setters
// they use. Budget is the plain budget row, the RISK_PARITY case whose Spinu point is the answer.
enum class Variant
{
    Free,
    LongOnly,
    Box,
    Leverage,
    Cash,
    TargetReturn,
    Budget,
    TargetVolatility
};

inline constexpr double kLeverage = 1.5;
inline constexpr double kTargetReturn = 0.08;
// A volatility cap or target at this multiple of the minimum volatility of the feasible set.
inline constexpr double kVolatilityFactor = 1.5;
// The RISK_PARITY box is this multiple of 1 / n: a cap that cuts the Spinu point.
inline constexpr double kRiskParityBoxFactor = 1.5;
inline constexpr std::size_t kObservationsPerAsset = 2;

struct Family
{
    po::Objective objective;
    Variant variant;
};

struct Case
{
    std::string name;
    po::Spec spec;
    const po::Data* data;
};

inline const char* objective_name(po::Objective objective)
{
    switch (objective)
    {
    case po::Objective::MinVariance:
        return "min_variance";
    case po::Objective::MaxReturn:
        return "max_return";
    case po::Objective::MaxSharpe:
        return "max_sharpe";
    case po::Objective::Utility:
        return "utility";
    case po::Objective::Importance:
        return "importance";
    case po::Objective::Kelly:
        return "kelly";
    case po::Objective::TrackingError:
        return "tracking_error";
    case po::Objective::MaxDiversification:
        return "max_diversification";
    case po::Objective::RiskParity:
        return "risk_parity";
    }
    throw std::logic_error("unknown objective");
}

inline const char* variant_name(Variant variant)
{
    switch (variant)
    {
    case Variant::Free:
        return "free";
    case Variant::LongOnly:
        return "long_only";
    case Variant::Box:
        return "box";
    case Variant::Leverage:
        return "leverage";
    case Variant::Cash:
        return "cash";
    case Variant::TargetReturn:
        return "target_return";
    case Variant::Budget:
        return "budget";
    case Variant::TargetVolatility:
        return "target_volatility";
    }
    throw std::logic_error("unknown variant");
}

inline std::vector<Family> families()
{
    using po::Objective;
    std::vector<Family> list;
    const std::vector<Variant> general = {Variant::Free, Variant::LongOnly, Variant::Box, Variant::Leverage,
                                          Variant::Cash};
    for (Objective objective : {Objective::MinVariance, Objective::MaxReturn, Objective::MaxSharpe,
                                Objective::Utility, Objective::TrackingError, Objective::MaxDiversification})
    {
        for (Variant variant : general)
        {
            list.push_back({objective, variant});
        }
    }
    list.push_back({Objective::MinVariance, Variant::TargetReturn});
    for (Variant variant : {Variant::Budget, Variant::TargetReturn, Variant::TargetVolatility, Variant::Box})
    {
        list.push_back({Objective::RiskParity, variant});
    }
    return list;
}

// The sizes a family runs at. n = 1000 is only the long-only family of each objective (and the plain RISK_PARITY
// case): the other variants at that size take minutes per solve.
inline std::vector<std::size_t> sizes_of(const Family& family)
{
    std::vector<std::size_t> sizes = {10, 100};
    bool large = family.variant == Variant::LongOnly || family.variant == Variant::Budget;
    if (large)
    {
        sizes.push_back(1000);
    }
    return sizes;
}

inline double box_cap(const Family& family, std::size_t n)
{
    if (family.objective == po::Objective::RiskParity)
    {
        return kRiskParityBoxFactor / static_cast<double>(n);
    }
    return std::max(0.3, 2.0 / static_cast<double>(n));
}

inline void apply_variant(const Family& family, std::size_t n, po::Spec& spec, DataOptions& options)
{
    switch (family.variant)
    {
    case Variant::Free:
    case Variant::Budget:
    case Variant::TargetVolatility:
        return;
    case Variant::LongOnly:
        spec.lower = 0.0;
        return;
    case Variant::Box:
        spec.lower = 0.0;
        spec.upper = box_cap(family, n);
        return;
    case Variant::Leverage:
        spec.leverage = kLeverage;
        return;
    case Variant::Cash:
        spec.lower = 0.0;
        options.cash = true;
        return;
    case Variant::TargetReturn:
        spec.lower = 0.0;
        spec.target_return = kTargetReturn;
        return;
    }
}

// The volatility of the minimum-variance portfolio over the bounds and the leverage of spec, without cash.
inline double minimum_volatility(const po::Spec& spec, std::size_t n)
{
    po::Spec probe = spec;
    probe.objective = po::Objective::MinVariance;
    probe.target_return.reset();
    const po::Data& data = data_of(n, DataOptions{});
    po::Result result = po::solve(probe, data).value();
    if (result.status != po::SolveStatus::Optimal)
    {
        throw std::runtime_error("the minimum-variance probe of a volatility target is not optimal");
    }
    return po::metrics::portfolio_volatility(data, result.weights);
}

inline void apply_objective(const Family& family, std::size_t n, po::Spec& spec, DataOptions& options)
{
    if (family.objective == po::Objective::TrackingError)
    {
        options.observations = kObservationsPerAsset * n;
    }
    if (family.objective == po::Objective::MaxReturn)
    {
        spec.max_volatility = kVolatilityFactor * minimum_volatility(spec, n);
    }
    if (family.variant == Variant::TargetVolatility)
    {
        spec.target_volatility = kVolatilityFactor * minimum_volatility(spec, n);
    }
}

inline std::string case_name(const Family& family, std::size_t n)
{
    return std::format("{}/{}/{}", objective_name(family.objective), variant_name(family.variant), n);
}

inline Case build_case(const Family& family, std::size_t n)
{
    po::Spec spec;
    spec.objective = family.objective;
    DataOptions options;
    apply_variant(family, n, spec, options);
    apply_objective(family, n, spec, options);
    return Case{case_name(family, n), spec, &data_of(n, options)};
}

// The case of a family at a size, built once; the reference stays valid for the whole process.
inline const Case& case_of(const Family& family, std::size_t n)
{
    static std::map<std::string, Case> cache;
    std::string name = case_name(family, n);
    auto found = cache.find(name);
    if (found == cache.end())
    {
        found = cache.emplace(name, build_case(family, n)).first;
    }
    return found->second;
}

} // namespace po::bench
