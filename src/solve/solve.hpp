#pragma once

#include <expected>

#include "core/data.hpp"
#include "core/failure.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"

namespace po
{

std::expected<Result, Failure> solve(const Spec& spec, const Data& data);

} // namespace po
