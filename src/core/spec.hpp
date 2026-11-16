#pragma once

#include <limits>
#include <optional>
#include <vector>

#include "core/objective.hpp"

namespace po
{

struct Spec
{
    Objective objective = Objective::MinVariance;
    double lower = -std::numeric_limits<double>::infinity();
    double upper = std::numeric_limits<double>::infinity();
    // Non-empty replaces the scalar bound; one entry per weight, cash last.
    std::vector<double> asset_lower;
    std::vector<double> asset_upper;
    std::optional<double> target_return;
    std::optional<double> risk_aversion;
    std::optional<double> importance;
    std::optional<double> max_volatility;
    // sum_i |w_i| <= leverage over the risky weights; cash is not counted.
    std::optional<double> leverage;
    // RISK_PARITY: one entry per asset, empty means 1 / n.
    std::vector<double> risk_budget;
    std::optional<double> target_volatility;
};

} // namespace po
