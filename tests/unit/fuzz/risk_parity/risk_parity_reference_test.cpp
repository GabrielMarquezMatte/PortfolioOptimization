#include <gtest/gtest-spi.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <vector>

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "fuzz/risk_parity/risk_parity_reference.hpp"
#include "random_covariance.hpp"
#include "riskparity/spinu.hpp"
#include "user_kkt.hpp"

namespace
{

using po::Data;
using po::Result;
using po::SolveStatus;
using po::Spec;
using po::testing::ScanResult;
using po::testing::random_covariance;
using po::testing::random_vector;
using po::testing::risk_contributions;
using po::testing::risk_parity_gradient;
using po::testing::risk_parity_value;
using po::testing::scan_reference;
using po::testing::two_asset_risk_parity;

std::vector<double> normalized(std::vector<double> values)
{
    double total = std::accumulate(values.begin(), values.end(), 0.0);
    for (double& value : values)
    {
        value /= total;
    }
    return values;
}

Data covariance_data(std::size_t n, const std::vector<double>& mu, const std::vector<double>& covariance)
{
    return Data::create(n, mu, covariance).value();
}

std::vector<double> diagonal_covariance(const std::vector<double>& sigma)
{
    std::size_t n = sigma.size();
    std::vector<double> covariance(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        covariance[i * n + i] = sigma[i] * sigma[i];
    }
    return covariance;
}

std::vector<double> diagonal_solution(const std::vector<double>& sigma, const std::vector<double>& budget)
{
    std::vector<double> weights(sigma.size());
    for (std::size_t i = 0; i < sigma.size(); ++i)
    {
        weights[i] = std::sqrt(budget[i]) / sigma[i];
    }
    return normalized(weights);
}

double central_difference(const Data& data, const std::vector<double>& budget, std::vector<double> weights,
                          std::size_t coordinate)
{
    constexpr double kStep = 1e-6;
    double original = weights[coordinate];
    weights[coordinate] = original + kStep;
    double above = risk_parity_value(data, budget, weights);
    weights[coordinate] = original - kStep;
    double below = risk_parity_value(data, budget, weights);
    return (above - below) / (2.0 * kStep);
}

TEST(RiskParityReference, GradientMatchesCentralDifferences)
{
    for (std::size_t n : {std::size_t{5}, std::size_t{12}})
    {
        for (std::uint64_t seed = 1; seed <= 20; ++seed)
        {
            SCOPED_TRACE(seed);
            Data data = covariance_data(n, {}, random_covariance(n, seed, 100.0));
            std::vector<double> budget = normalized(random_vector(n, seed + 100, 0.2, 1.0));
            std::vector<double> weights = normalized(random_vector(n, seed + 200, 0.1, 1.0));
            std::vector<double> gradient = risk_parity_gradient(data, budget, weights);
            double scale = po::testing::largest_magnitude(gradient);
            for (std::size_t i = 0; i < n; ++i)
            {
                double numeric = central_difference(data, budget, weights, i);
                EXPECT_LE(std::fabs(numeric - gradient[i]), 1e-6 * scale);
            }
        }
    }
}

TEST(RiskParityReference, DiagonalClosedForm)
{
    for (std::uint64_t seed = 1; seed <= 20; ++seed)
    {
        SCOPED_TRACE(seed);
        std::size_t n = 3 + seed % 8;
        std::vector<double> sigma = random_vector(n, seed, 0.05, 0.6);
        std::vector<double> budget = normalized(random_vector(n, seed + 100, 0.01, 1.0));
        Data data = covariance_data(n, {}, diagonal_covariance(sigma));
        std::vector<double> weights = diagonal_solution(sigma, budget);
        EXPECT_LE(risk_parity_value(data, budget, weights), 1e-14);
    }
}

TEST(RiskParityReference, TwoAssetsClosedFormHasZeroValue)
{
    constexpr double kSigmaFirst = 0.2;
    constexpr double kSigmaSecond = 0.35;
    for (double rho : {-0.5, 0.0, 0.9})
    {
        for (double first_budget : {0.5, 0.2, 0.01})
        {
            SCOPED_TRACE(rho);
            SCOPED_TRACE(first_budget);
            double cross = rho * kSigmaFirst * kSigmaSecond;
            Data data = covariance_data(2, {},
                                        {kSigmaFirst * kSigmaFirst, cross, cross, kSigmaSecond * kSigmaSecond});
            std::vector<double> budget = {first_budget, 1.0 - first_budget};
            std::vector<double> weights = two_asset_risk_parity(kSigmaFirst, kSigmaSecond, rho, first_budget);
            ASSERT_EQ(weights.size(), 2u);
            EXPECT_GT(weights[0], 0.0);
            EXPECT_GT(weights[1], 0.0);
            EXPECT_NEAR(weights[0] + weights[1], 1.0, 1e-15);
            EXPECT_LE(risk_parity_value(data, budget, weights), 1e-14);
        }
    }
}

TEST(RiskParityReference, ContributionsSumToOne)
{
    for (std::uint64_t seed = 1; seed <= 20; ++seed)
    {
        SCOPED_TRACE(seed);
        std::size_t n = 4 + seed % 9;
        Data data = covariance_data(n, {}, random_covariance(n, seed, 1e4));
        std::vector<double> weights = normalized(random_vector(n, seed + 200, -0.3, 1.0));
        std::vector<double> contributions = risk_contributions(data, weights);
        EXPECT_NEAR(std::accumulate(contributions.begin(), contributions.end(), 0.0), 1.0, 1e-14);
    }
}

TEST(RiskParityReference, ScanFindsKnownMinimumOnTheReturnLine)
{
    for (std::uint64_t seed = 1; seed <= 10; ++seed)
    {
        SCOPED_TRACE(seed);
        std::vector<double> sigma = random_vector(3, seed, 0.1, 0.5);
        std::vector<double> budget = normalized(random_vector(3, seed + 100, 0.1, 1.0));
        std::vector<double> mu = random_vector(3, seed + 300, 0.03, 0.15);
        Data data = covariance_data(3, mu, diagonal_covariance(sigma));
        std::vector<double> known = diagonal_solution(sigma, budget);
        Spec spec;
        spec.target_return = po::testing::dot(known, mu);
        ScanResult scan = scan_reference(spec, data, budget);
        EXPECT_LE(scan.value, 1e-9);
        ASSERT_EQ(scan.weights.size(), 3u);
        EXPECT_NEAR(scan.weights[0] + scan.weights[1] + scan.weights[2], 1.0, 1e-12);
    }
}

TEST(RiskParityReference, ScanFindsKnownMinimumOnTheVolatilityEllipse)
{
    for (std::uint64_t seed = 1; seed <= 10; ++seed)
    {
        SCOPED_TRACE(seed);
        std::vector<double> sigma = random_vector(3, seed, 0.1, 0.5);
        std::vector<double> budget = normalized(random_vector(3, seed + 100, 0.1, 1.0));
        Data data = covariance_data(3, {}, diagonal_covariance(sigma));
        std::vector<double> known = diagonal_solution(sigma, budget);
        Spec spec;
        spec.target_volatility = po::testing::portfolio_volatility(data, known);
        ScanResult scan = scan_reference(spec, data, budget);
        EXPECT_LE(scan.value, 1e-9);
    }
}

TEST(RiskParityReference, ScanDiscardsPointsOutsideTheBox)
{
    std::vector<double> sigma = {0.1, 0.3, 0.5};
    std::vector<double> budget = {0.2, 0.3, 0.5};
    Data data = covariance_data(3, {0.05, 0.08, 0.11}, diagonal_covariance(sigma));
    std::vector<double> known = diagonal_solution(sigma, budget);
    Spec spec;
    spec.target_return = po::testing::dot(known, data.mu());
    spec.asset_lower = {known[0] + 0.05, -1.0, -1.0};
    ScanResult scan = scan_reference(spec, data, budget);
    ASSERT_FALSE(std::isnan(scan.value));
    EXPECT_GT(scan.value, 1e-9);
    for (double weight : scan.weights)
    {
        EXPECT_GT(weight, -1.0 - 1e-12);
    }
    EXPECT_GE(scan.weights[0], spec.asset_lower[0] - 1e-12);
}

TEST(RiskParityReference, ScanWithoutAFeasiblePointIsNaN)
{
    std::vector<double> sigma = {0.1, 0.3, 0.5};
    std::vector<double> budget = {0.2, 0.3, 0.5};
    Data data = covariance_data(3, {0.05, 0.08, 0.11}, diagonal_covariance(sigma));
    Spec spec;
    spec.target_return = 0.5;
    spec.lower = 0.0;
    spec.upper = 1.0;
    EXPECT_TRUE(std::isnan(scan_reference(spec, data, budget).value));
}

// Weights pinned at w (lower == upper) and multipliers chosen to satisfy the stationarity exactly: the KKT check
// accepts a point built to be stationary, with and without the volatility equality.
TEST(RiskParityReference, KktAcceptsAPinnedStationaryPoint)
{
    for (bool with_volatility : {false, true})
    {
        SCOPED_TRACE(with_volatility);
        std::size_t n = 6;
        Data data = covariance_data(n, random_vector(n, 7, 0.03, 0.15), random_covariance(n, 3, 100.0));
        std::vector<double> budget = normalized(random_vector(n, 11, 0.2, 1.0));
        std::vector<double> weights = normalized(random_vector(n, 13, 0.1, 1.0));
        std::vector<double> gradient = risk_parity_gradient(data, budget, weights);
        std::vector<double> product = po::testing::covariance_times(data, weights);
        double volatility = po::testing::portfolio_volatility(data, weights);
        constexpr double kVolatilityMultiplier = 0.3;
        Spec spec;
        spec.risk_budget = budget;
        spec.asset_lower = weights;
        spec.asset_upper = weights;
        Result result;
        result.status = SolveStatus::Optimal;
        result.weights = weights;
        result.budget_multiplier = 0.1;
        result.bound_multipliers.resize(n);
        if (with_volatility)
        {
            spec.target_volatility = volatility;
            result.volatility_multiplier = kVolatilityMultiplier;
        }
        for (std::size_t i = 0; i < n; ++i)
        {
            result.bound_multipliers[i] = gradient[i] - result.budget_multiplier;
            if (with_volatility)
            {
                result.bound_multipliers[i] -= kVolatilityMultiplier * product[i] / volatility;
            }
        }
        po::testing::expect_risk_parity_kkt(spec, data, result, 1e-10);
    }
}

// A point near the Spinu solution (gap about 1e-3), pinned by its bounds, with the exact stationary multipliers and
// then an error of 1e-8 of the gradient in one of them: relative to the terms it is 1e-8, far above the 1e-10
// tolerance. A scale of the whole problem (1 / gap larger) would hide it.
TEST(RiskParityReference, KktRejectsAStationarityErrorAwayFromTheFloor)
{
    std::size_t n = 6;
    Data data = covariance_data(n, random_vector(n, 7, 0.03, 0.15), random_covariance(n, 3, 100.0));
    std::vector<double> budget(n, 1.0 / static_cast<double>(n));
    std::vector<double> weights = po::riskparity::spinu(data, budget).weights;
    std::vector<double> nudge = random_vector(n, 17, -1e-3, 1e-3);
    for (std::size_t i = 0; i < n; ++i)
    {
        weights[i] *= 1.0 + nudge[i];
    }
    weights = normalized(weights);
    std::vector<double> gradient = risk_parity_gradient(data, budget, weights);
    double largest = 0.0;
    Spec spec;
    spec.asset_lower = weights;
    spec.asset_upper = weights;
    Result result;
    result.status = SolveStatus::Optimal;
    result.weights = weights;
    result.budget_multiplier = 0.0;
    result.bound_multipliers = gradient;
    for (double entry : gradient)
    {
        largest = std::fmax(largest, std::fabs(entry));
    }
    ASSERT_GT(po::testing::contribution_gap(data, budget, weights), 1e-4);
    ASSERT_GT(po::testing::unit_gap_gradient_scale(data, weights), 100.0 * largest);
    po::testing::expect_risk_parity_kkt(spec, data, result, 1e-10);
    result.bound_multipliers[0] += 1e-8 * largest;
    EXPECT_NONFATAL_FAILURE(po::testing::expect_risk_parity_kkt(spec, data, result, 1e-10), "tolerance");
}

TEST(RiskParityReference, ShapesCoverEveryCombination)
{
    std::vector<po::testing::RiskParityShape> shapes = po::testing::risk_parity_shapes(192, 5);
    ASSERT_EQ(shapes.size(), 192u);
    for (const po::testing::RiskParityShape& shape : shapes)
    {
        SCOPED_TRACE(po::testing::describe(shape));
        po::testing::GeneratedRiskParity problem = po::testing::generate(shape);
        EXPECT_EQ(problem.data.asset_count(), shape.n);
        ASSERT_EQ(problem.spec.risk_budget.size(), shape.n);
        EXPECT_NEAR(std::accumulate(problem.spec.risk_budget.begin(), problem.spec.risk_budget.end(), 0.0), 1.0,
                    1e-14);
        EXPECT_EQ(problem.spec.target_return.has_value(), shape.target_return);
        EXPECT_EQ(problem.spec.target_volatility.has_value(), shape.target_volatility);
    }
}

} // namespace
