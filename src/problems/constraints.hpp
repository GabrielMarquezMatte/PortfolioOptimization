#pragma once

#include <cstddef>
#include <optional>
#include <vector>

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "qp/problem.hpp"

namespace po::problems
{

struct BoundRow
{
    std::size_t weight;
    std::size_t row;
    bool lower;
};

struct ConstraintRows
{
    std::optional<std::size_t> budget;
    std::optional<std::size_t> target_return;
    std::vector<BoundRow> bounds;
    // The spec set a leverage, and the QP has its gross exposure family.
    bool gross_exposure = false;
};

// The QP variables are the risky weights only; with cash, the cash weight 1 - 1'w is eliminated.
struct Translation
{
    qp::Problem problem;
    ConstraintRows rows;
};

ConstraintRows add_constraint_rows(const Spec& spec, const Data& data, qp::Problem& problem);

// The QP objective must equal multiplier_scale times the minimization-form objective of Result, up to a
// constant, so that the QP multipliers divided by multiplier_scale are the user's.
Result to_result(const Data& data, const ConstraintRows& rows, const qp::Solution& solution,
                 double multiplier_scale = 1.0);

} // namespace po::problems
