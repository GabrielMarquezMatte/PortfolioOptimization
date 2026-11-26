#pragma once

#include <cstddef>
#include <expected>
#include <vector>

#include "core/data.hpp"
#include "core/failure.hpp"

namespace po::estimation
{

struct BayesStein
{
    std::vector<double> mu;
    double grand_mean;
    double weight;
};

// Jorion (1986) shrinkage of data.mu() toward the minimum-variance mean. data.covariance() is the sample
// covariance with divisor t - 1 of t observations; needs t > n + 2. The weight depends on the frequency of
// the data (scaling mu and the covariance by 252 scales the quadratic form by 252), not on its unit.
std::expected<BayesStein, Failure> bayes_stein(const Data& data, std::size_t t);

} // namespace po::estimation
