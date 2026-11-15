#pragma once

#include <cstdint>

namespace po
{

enum class Objective : std::int32_t
{
    MinVariance = 0,
    MaxReturn,
    MaxSharpe,
    Utility,
    Importance,
    Kelly,
    TrackingError,
    MaxDiversification,
    RiskParity
};

} // namespace po
