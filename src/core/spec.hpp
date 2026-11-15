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
};

} // namespace po
