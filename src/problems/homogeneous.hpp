#pragma once

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"

namespace po::problems
{

// f = -v'w / sqrt(w'Sigma w) over 1'w = 1, the bounds of the risky weights and the leverage, with
// v = mu - rf 1 (MAX_SHARPE, rf = 0 without cash) or sigma = sqrt(diag Sigma) (MAX_DIVERSIFICATION). With
// cash the cash weight is fixed at 0, so 1'w = 1 holds over the risky weights; cash bounds that exclude 0 are
// INFEASIBLE before solving. Solved by homogenization, y = kappa w with kappa = 1'y > 0 relaxed to 1'y >= 0:
//   min y'Sigma y  s.t.  v'y = 1,  y_i - lo_i 1'y >= 0,  hi_i 1'y - y_i >= 0 (finite bounds),  1'y >= 0,
// and with leverage L 1'y - ||y||_1 >= 0 (gross_exposure {L 1, 0}); G = 2 Sigma, factor sqrt(2) L.
// A QP that is Infeasible is INFEASIBLE (no feasible w has v'w > 0). The row 1'y >= 0 in the active set at
// the optimum is UNBOUNDED: the supremum is approached as kappa -> 0 and never attained. On an unbounded cone (no
// leverage and some i != j with hi_i = +inf, lo_j = -inf) so is a kappa = 1'y within 16 times its rounding, whose
// sign the data do not decide (near 1' Sigma^-1 v = 0 without bounds); on a bounded one kappa cannot reach 0 and only
// kappa <= 0 counts. Otherwise OPTIMAL with w = y / kappa.
Result solve_max_sharpe(const Spec& spec, const Data& data);
Result solve_max_diversification(const Spec& spec, const Data& data);

} // namespace po::problems
