#include "riskparity/objective.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <vector>

#include "core/data.hpp"
#include "fuzz/risk_parity/risk_parity_reference.hpp"
#include "metrics/portfolio.hpp"
#include "random_covariance.hpp"

namespace
{

using po::Data;
using po::riskparity::Linearization;
using po::riskparity::linearize;
using po::testing::largest_magnitude;
using po::testing::random_covariance;
using po::testing::random_vector;

constexpr double kStep = 1e-6;

Data sample_data(std::size_t n, std::uint64_t seed)
{
    return Data::create(n, random_vector(n, seed + 500, 0.03, 0.15), random_covariance(n, seed, 1e4)).value();
}

std::vector<double> normalized(std::vector<double> values)
{
    double total = std::accumulate(values.begin(), values.end(), 0.0);
    for (double& value : values)
    {
        value /= total;
    }
    return values;
}

// A point of the budget row with a short position, away from the risk parity solution.
std::vector<double> sample_weights(std::size_t n, std::uint64_t seed)
{
    return normalized(random_vector(n, seed + 600, -0.2, 1.0));
}

std::vector<double> sample_budget(std::size_t n, std::uint64_t seed)
{
    return normalized(random_vector(n, seed + 700, 0.2, 1.0));
}

std::vector<double> gradient_at(const Data& data, const std::vector<double>& budget,
                                const std::vector<double>& weights)
{
    return po::testing::risk_parity_gradient(data, budget, weights);
}

std::vector<double> moved(std::vector<double> weights, std::size_t index, double step)
{
    weights[index] += step;
    return weights;
}

TEST(RiskParityObjective, GradientMatchesIndependentReference)
{
    for (std::size_t n : {2u, 7u, 30u})
    {
        for (std::uint64_t seed = 1; seed <= 5; ++seed)
        {
            SCOPED_TRACE(testing::Message() << "n " << n << " seed " << seed);
            Data data = sample_data(n, seed);
            std::vector<double> budget = sample_budget(n, seed);
            std::vector<double> weights = sample_weights(n, seed);
            Linearization linearization = linearize(data, budget, weights, false);
            std::vector<double> reference = po::testing::risk_parity_gradient(data, budget, weights);
            double largest_error = 0.0;
            for (std::size_t j = 0; j < n; ++j)
            {
                largest_error = std::fmax(largest_error, std::fabs(linearization.gradient[j] - reference[j]));
            }
            EXPECT_LE(largest_error, 1e-12 * largest_magnitude(reference));
            EXPECT_TRUE(linearization.hessian.empty());
        }
    }
}

TEST(RiskParityObjective, ValueMatchesMetric)
{
    for (std::size_t n : {2u, 7u, 30u})
    {
        SCOPED_TRACE(testing::Message() << "n " << n);
        Data data = sample_data(n, 3);
        std::vector<double> budget = sample_budget(n, 3);
        std::vector<double> weights = sample_weights(n, 3);
        EXPECT_EQ(linearize(data, budget, weights, false).value,
                  po::metrics::risk_budget_residual(data, weights, budget));
    }
}

TEST(RiskParityObjective, JacobianByCentralDifferences)
{
    for (std::size_t n : {3u, 10u})
    {
        SCOPED_TRACE(testing::Message() << "n " << n);
        Data data = sample_data(n, 11);
        std::vector<double> budget = sample_budget(n, 11);
        std::vector<double> weights = sample_weights(n, 11);
        std::vector<double> jacobian = linearize(data, budget, weights, false).jacobian;
        ASSERT_EQ(jacobian.size(), n * n);
        double largest_error = 0.0;
        for (std::size_t j = 0; j < n; ++j)
        {
            std::vector<double> ahead = po::testing::risk_contributions(data, moved(weights, j, kStep));
            std::vector<double> behind = po::testing::risk_contributions(data, moved(weights, j, -kStep));
            for (std::size_t i = 0; i < n; ++i)
            {
                double difference = (ahead[i] - behind[i]) / (2.0 * kStep);
                largest_error = std::fmax(largest_error, std::fabs(jacobian[i * n + j] - difference));
            }
        }
        EXPECT_LE(largest_error, 1e-6 * largest_magnitude(jacobian));
        for (std::size_t i = 0; i < n; ++i)
        {
            double product = 0.0;
            double magnitude = 0.0;
            for (std::size_t j = 0; j < n; ++j)
            {
                product += jacobian[i * n + j] * weights[j];
                magnitude += std::fabs(jacobian[i * n + j] * weights[j]);
            }
            EXPECT_LE(std::fabs(product), 1e-12 * magnitude) << "row " << i;
        }
    }
}

TEST(RiskParityObjective, HessianByCentralDifferencesOfGradient)
{
    for (std::size_t n : {3u, 10u})
    {
        SCOPED_TRACE(testing::Message() << "n " << n);
        Data data = sample_data(n, 13);
        std::vector<double> budget = sample_budget(n, 13);
        std::vector<double> weights = sample_weights(n, 13);
        std::vector<double> hessian = linearize(data, budget, weights, true).hessian;
        ASSERT_EQ(hessian.size(), n * n);
        double largest_error = 0.0;
        for (std::size_t k = 0; k < n; ++k)
        {
            std::vector<double> ahead = gradient_at(data, budget, moved(weights, k, kStep));
            std::vector<double> behind = gradient_at(data, budget, moved(weights, k, -kStep));
            for (std::size_t j = 0; j < n; ++j)
            {
                double difference = (ahead[j] - behind[j]) / (2.0 * kStep);
                largest_error = std::fmax(largest_error, std::fabs(hessian[j * n + k] - difference));
            }
        }
        EXPECT_LE(largest_error, 1e-6 * largest_magnitude(hessian));
    }
}

// f(w + s) - f(w) against its second-order expansion g's + 1/2 s'Hs, whose error is of order |s|^3: a
// difference of two values of f would lose every digit below the rounding of f, about 1e-18 here.
TEST(RiskParityObjective, ValueChangeKeepsRelativeAccuracy)
{
    std::size_t n = 7;
    Data data = sample_data(n, 17);
    std::vector<double> budget = sample_budget(n, 17);
    std::vector<double> weights = sample_weights(n, 17);
    Linearization linearization = linearize(data, budget, weights, true);
    std::vector<double> direction = random_vector(n, 900, -1.0, 1.0);
    for (double length : {1e-6, 1e-9, 1e-12, 1e-15})
    {
        SCOPED_TRACE(testing::Message() << "length " << length);
        std::vector<double> step(n);
        for (std::size_t i = 0; i < n; ++i)
        {
            step[i] = length * direction[i];
        }
        double slope = po::testing::dot(linearization.gradient, step);
        double curvature = 0.0;
        for (std::size_t j = 0; j < n; ++j)
        {
            curvature += step[j] * po::testing::row_dot(linearization.hessian, j, step);
        }
        double expansion = slope + 0.5 * curvature;
        po::riskparity::Change change = po::riskparity::value_change(data, budget, weights, step);
        EXPECT_LE(std::fabs(change.value - expansion), 1e-9 * std::fabs(slope));
        double variance_change = 2.0 * po::testing::dot(po::testing::covariance_times(data, weights), step) +
                                 po::testing::dot(po::testing::covariance_times(data, step), step);
        EXPECT_NEAR(change.variance, variance_change, 1e-14 * std::fabs(variance_change));
    }
}

} // namespace
