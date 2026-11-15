#pragma once

#include <optional>

#include "core/data.hpp"
#include "core/failure.hpp"
#include "core/spec.hpp"

namespace po::problems
{

// In order: the objective is implemented, every setter applies to it, the data it requires is present, and the
// per-weight bounds have one entry per weight.
std::optional<Failure> check_request(const Spec& spec, const Data& data);

} // namespace po::problems
