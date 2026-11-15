#pragma once

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "problems/constraints.hpp"

namespace po::problems
{

Translation min_variance_problem(const Spec& spec, const Data& data);

Result solve_min_variance(const Spec& spec, const Data& data);

} // namespace po::problems
