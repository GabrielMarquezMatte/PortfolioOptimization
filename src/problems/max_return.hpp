#pragma once

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"

namespace po::problems
{

// f = -(mu - rf 1)'w subject to sqrt(w'Sigma w) <= s, with rf = 0 without cash and s = spec.max_volatility.
// Solved along Q(t): min 0.5 w'Sigma w - t (mu - rf 1)'w over the same constraint rows, whose optimum w(t) is
// piecewise affine in t with a nondecreasing volatility. At the t* where the volatility reaches s, the user
// multipliers are the QP multipliers over t* and the volatility multiplier is -s / t*. A volatility that
// never reaches s is inactive: its multiplier is 0. s below the minimum volatility is Infeasible.
Result solve_max_return(const Spec& spec, const Data& data);

} // namespace po::problems
