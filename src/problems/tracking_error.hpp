#pragma once

#include <expected>

#include "core/data.hpp"
#include "core/failure.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "problems/constraints.hpp"

namespace po::problems
{

// f = |R~ w + rf 1 - b|^2 with R~ = R - rf 1 1' (rf = 0 without cash), so that G = 2 R~'R~ and
// g = -2 R~'(b - rf 1). The Hessian factor is sqrt(2) times the Cholesky factor of R~'R~, which is the one
// decision on the column rank of R~: Failure{InvalidArgument} when it is not positive definite.
std::expected<Translation, Failure> tracking_error_problem(const Spec& spec, const Data& data);

std::expected<Result, Failure> solve_tracking_error(const Spec& spec, const Data& data);

} // namespace po::problems
