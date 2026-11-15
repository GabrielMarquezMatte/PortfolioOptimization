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
// g = -2 R~'(b - rf 1). R~ is factored by Householder QR, R~ = Q R_qr, never through R~'R~, whose condition
// number is the square of that of R~: the Hessian factor is sqrt(2) R_qr', and the one decision on the column
// rank of R~ is on the diagonal of R_qr (linalg::has_full_column_rank). Failure{InvalidArgument} when it
// fails or when there are fewer periods than assets.
std::expected<Translation, Failure> tracking_error_problem(const Spec& spec, const Data& data);

std::expected<Result, Failure> solve_tracking_error(const Spec& spec, const Data& data);

} // namespace po::problems
