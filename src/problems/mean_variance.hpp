#pragma once

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "problems/constraints.hpp"

namespace po::problems
{

// f = (quadratic / 2) w'Sigma w - linear (mu - rf 1)'w, with rf = 0 without cash.
struct MeanVariance
{
    double quadratic;
    double linear;
};

// MIN_VARIANCE (2, 0); UTILITY (delta, 1); IMPORTANCE (2 (1 - lambda), lambda); KELLY (1 / (1 + rf), 1).
MeanVariance mean_variance_of(const Spec& spec, const Data& data);

Translation mean_variance_problem(const Spec& spec, const Data& data, MeanVariance scales);

Result solve_mean_variance(const Spec& spec, const Data& data);

} // namespace po::problems
