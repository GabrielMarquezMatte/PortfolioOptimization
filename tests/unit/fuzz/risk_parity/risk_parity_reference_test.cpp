#include <gtest/gtest-spi.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "core/data.hpp"
#include "core/objective.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "fuzz/frontier_reference.hpp"
#include "fuzz/risk_parity/risk_parity_reference.hpp"
#include "fuzz/risk_parity/risk_parity_refinement.hpp"
#include "random_covariance.hpp"
#include "riskparity/objective.hpp"
#include "riskparity/spinu.hpp"
#include "solve/solve.hpp"
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
    if (!po::testing::kLongDoubleIsExtended)
    {
        GTEST_SKIP() << "long double is double here: the oracle has no extra precision";
    }
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

// expect_risk_parity_kkt fails on the result, with a failure that mentions the reason.
void expect_rejected(const Spec& spec, const Data& data, const Result& result, const std::string& reason)
{
    std::vector<std::string> failures = po::testing::risk_parity_kkt_failures(spec, data, result, 1e-10);
    auto mentions = [&reason](const std::string& failure)
    { return failure.find(reason) != std::string::npos; };
    EXPECT_TRUE(std::ranges::any_of(failures, mentions))
        << failures.size() << " failures, none about " << reason;
}

// A point held at its lower bounds (upper bounds 1 above), stationary by construction: beta = grad f -
// budget 1 - target_return mu, with the budget multiplier below every entry of the rest so that every beta is
// positive.
struct LowerBoundPoint
{
    Data data;
    Spec spec;
    Result result;
};

LowerBoundPoint lower_bound_point(bool with_target_return)
{
    std::size_t n = 6;
    Data data = covariance_data(n, random_vector(n, 7, 0.03, 0.15), random_covariance(n, 3, 100.0));
    std::vector<double> weights = normalized(random_vector(n, 13, 0.1, 1.0));
    Spec spec;
    spec.risk_budget = normalized(random_vector(n, 11, 0.2, 1.0));
    spec.asset_lower = weights;
    spec.asset_upper = weights;
    for (double& upper : spec.asset_upper)
    {
        upper += 1.0;
    }
    Result result;
    result.status = SolveStatus::Optimal;
    result.weights = weights;
    result.bound_multipliers = risk_parity_gradient(data, spec.risk_budget, weights);
    if (with_target_return)
    {
        spec.target_return = po::testing::dot(weights, data.mu());
        result.target_return_multiplier = 0.7;
        for (std::size_t i = 0; i < n; ++i)
        {
            result.bound_multipliers[i] -= result.target_return_multiplier * data.mu()[i];
        }
    }
    result.budget_multiplier = *std::ranges::min_element(result.bound_multipliers) - 0.1;
    for (double& beta : result.bound_multipliers)
    {
        beta -= result.budget_multiplier;
    }
    return LowerBoundPoint{std::move(data), std::move(spec), std::move(result)};
}

// The scale of the stationarity check at the point: the largest term of grad f or of the constraints.
double stationarity_scale(const LowerBoundPoint& point)
{
    po::testing::ObjectiveGradient gradient = po::testing::sum_of(
        po::testing::risk_parity_gradient_terms(point.data, point.spec.risk_budget, point.result.weights));
    return std::fmax(gradient.scale, po::testing::constraint_terms(point.data, point.result).scale);
}

// An error of 1e-8 of the scale of the check in one multiplier, 100 times the tolerance.
TEST(RiskParityReference, KktRejectsAPerturbedBudgetMultiplier)
{
    if (!po::testing::kLongDoubleIsExtended)
    {
        GTEST_SKIP() << "long double is double here: the oracle has no extra precision";
    }
    LowerBoundPoint point = lower_bound_point(false);
    po::testing::expect_risk_parity_kkt(point.spec, point.data, point.result, 1e-10);
    point.result.budget_multiplier += 1e-8 * stationarity_scale(point);
    expect_rejected(point.spec, point.data, point.result, "tolerance");
}

TEST(RiskParityReference, KktRejectsAPerturbedTargetReturnMultiplier)
{
    if (!po::testing::kLongDoubleIsExtended)
    {
        GTEST_SKIP() << "long double is double here: the oracle has no extra precision";
    }
    LowerBoundPoint point = lower_bound_point(true);
    po::testing::expect_risk_parity_kkt(point.spec, point.data, point.result, 1e-10);
    double largest_mean = po::testing::largest_magnitude(point.data.mu());
    point.result.target_return_multiplier += 1e-8 * stationarity_scale(point) / largest_mean;
    expect_rejected(point.spec, point.data, point.result, "tolerance");
    point = lower_bound_point(true);
    point.result.target_return_multiplier = std::numeric_limits<double>::quiet_NaN();
    expect_rejected(point.spec, point.data, point.result, "target_return_multiplier");
}

// The same stationary point with the budget multiplier above every entry, so that every beta is negative at a
// lower bound: the stationarity holds exactly and only the sign of beta is wrong.
TEST(RiskParityReference, KktRejectsAWrongSignBoundMultiplier)
{
    LowerBoundPoint point = lower_bound_point(false);
    double shift = po::testing::largest_magnitude(point.result.bound_multipliers) + 0.1;
    point.result.budget_multiplier += shift;
    for (double& beta : point.result.bound_multipliers)
    {
        beta -= shift;
        ASSERT_LT(beta, 0.0);
    }
    po::testing::ObjectiveGradient gradient = po::testing::sum_of(
        po::testing::risk_parity_gradient_terms(point.data, point.spec.risk_budget, point.result.weights));
    po::testing::expect_user_stationarity(point.data, point.result, gradient, 1e-10);
    expect_rejected(point.spec, point.data, point.result, "bound_complementarity");
}

// A weight without a finite bound has beta = 0, never NaN: with the bounds of weight 0 removed, a NaN there
// must not pass as "no multiplier".
TEST(RiskParityReference, KktRejectsANaNBoundMultiplierAtAnInfiniteBound)
{
    LowerBoundPoint point = lower_bound_point(false);
    point.spec.asset_lower[0] = -std::numeric_limits<double>::infinity();
    point.spec.asset_upper[0] = std::numeric_limits<double>::infinity();
    point.result.bound_multipliers[0] = std::numeric_limits<double>::quiet_NaN();
    expect_rejected(point.spec, point.data, point.result, "tolerance");
}

// The result of s = sigma_min: the minimum-variance weights, a volatility multiplier of -inf and every other
// multiplier NaN.
struct MinimumVolatilityPoint
{
    Data data;
    Spec spec;
    Result result;
};

MinimumVolatilityPoint minimum_volatility_point()
{
    std::size_t n = 5;
    Data data = covariance_data(n, random_vector(n, 7, 0.03, 0.15), random_covariance(n, 3, 100.0));
    Spec spec;
    spec.lower = 0.0;
    spec.upper = 1.0;
    Result minimum = po::testing::frontier_point(spec, data, std::nullopt);
    spec.target_volatility = po::testing::portfolio_volatility(data, minimum.weights);
    Result result;
    result.status = SolveStatus::Optimal;
    result.weights = minimum.weights;
    result.bound_multipliers.assign(n, std::numeric_limits<double>::quiet_NaN());
    result.volatility_multiplier = -std::numeric_limits<double>::infinity();
    return MinimumVolatilityPoint{std::move(data), std::move(spec), std::move(result)};
}

TEST(RiskParityReference, KktAcceptsTheMinimumVolatilityPoint)
{
    MinimumVolatilityPoint point = minimum_volatility_point();
    po::testing::expect_risk_parity_kkt(point.spec, point.data, point.result, 1e-10);
}

TEST(RiskParityReference, KktRejectsAnInfiniteVolatilityMultiplierAwayFromTheMinimum)
{
    MinimumVolatilityPoint point = minimum_volatility_point();
    *point.spec.target_volatility *= 1.0 + 1e-9;
    expect_rejected(point.spec, point.data, point.result, "s is not the minimum volatility");
}

TEST(RiskParityReference, KktRejectsOtherWeightsOrMultipliersAtTheMinimumVolatility)
{
    MinimumVolatilityPoint moved = minimum_volatility_point();
    moved.result.weights[0] += 1e-6;
    moved.result.weights[1] -= 1e-6;
    expect_rejected(moved.spec, moved.data, moved.result, "weight 0");
    MinimumVolatilityPoint finite = minimum_volatility_point();
    finite.result.budget_multiplier = 0.0;
    expect_rejected(finite.spec, finite.data, finite.result, "budget_multiplier");
    MinimumVolatilityPoint bound = minimum_volatility_point();
    bound.result.bound_multipliers[2] = 0.0;
    expect_rejected(bound.spec, bound.data, bound.result, "bound multiplier 2");
}

// Away from the minimum the volatility multiplier is finite: +inf is no multiplier.
TEST(RiskParityReference, KktRejectsAPositiveInfiniteVolatilityMultiplier)
{
    MinimumVolatilityPoint point = minimum_volatility_point();
    point.result.volatility_multiplier = std::numeric_limits<double>::infinity();
    point.result.budget_multiplier = 0.0;
    point.result.bound_multipliers.assign(point.result.weights.size(), 0.0);
    expect_rejected(point.spec, point.data, point.result, "volatility_multiplier");
}

TEST(RiskParityReference, HessianByDifferencesMatchesTheObjective)
{
    if (!po::testing::kLongDoubleIsExtended)
    {
        GTEST_SKIP() << "long double is double here: the oracle has no extra precision";
    }
    std::size_t n = 6;
    for (double condition : {1e2, 1e6})
    {
        SCOPED_TRACE(condition);
        Data data = covariance_data(n, {}, random_covariance(n, 3, condition));
        std::vector<double> budget = normalized(random_vector(n, 11, 0.2, 1.0));
        std::vector<double> weights = normalized(random_vector(n, 13, -0.3, 1.0));
        std::vector<long double> differences =
            po::testing::wide_hessian(data, budget, po::testing::widened(weights));
        std::vector<double> hessian = po::riskparity::linearize(data, budget, weights, true).hessian;
        double scale = po::testing::largest_magnitude(hessian);
        for (std::size_t k = 0; k < n * n; ++k)
        {
            EXPECT_NEAR(static_cast<double>(differences[k]), hessian[k], 1e-6 * scale) << "entry " << k;
        }
    }
}

// The path the references take where long double is double (MSVC): the entries of Sigma w, which cancel by
// up to cond(Sigma), by Dot2 in double agree with Dot2 in long double to the rounding of each entry.
TEST(RiskParityReference, CompensatedDotInDoubleIsAccurateAtCondition1e8)
{
    std::size_t n = 20;
    Data data = covariance_data(n, {}, random_covariance(n, 5, 1e8));
    std::vector<double> weights = normalized(random_vector(n, 17, -0.3, 1.0));
    double worst = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        std::span<const double> row = data.covariance().subspan(i * n, n);
        double plain = po::testing::compensated_dot_of<double>(n, row, weights);
        long double wide = po::testing::compensated_dot_of<long double>(n, row, weights);
        worst = std::fmax(worst, static_cast<double>(std::fabs(plain - wide) / std::fabs(wide)));
    }
    EXPECT_LE(worst, std::numeric_limits<double>::epsilon());
}

// At cond(Sigma) = 1e8, where the floor check and the refinement rely on it without another check: a coarse
// band (the two differ by up to 3.2e-7 of the largest entry here). Not because the analytic Hessian of the
// library is rounded by eps cond(Sigma): against 80 digits it is right to 1.3e-15 of its largest entry at
// every cond(Sigma) up to 1e10 at the points of the multiplier measurement (docs/guide/validation.md), so
// the gap is that of the differences; at these random weights the two were not measured apart.
TEST(RiskParityReference, HessianByDifferencesMatchesTheObjectiveAtCondition1e8)
{
    std::size_t n = 6;
    for (std::uint64_t seed : {3u, 5u, 7u})
    {
        SCOPED_TRACE(seed);
        Data data = covariance_data(n, {}, random_covariance(n, seed, 1e8));
        std::vector<double> budget = normalized(random_vector(n, 11 + seed, 0.2, 1.0));
        std::vector<double> weights = normalized(random_vector(n, 13 + seed, -0.3, 1.0));
        std::vector<long double> differences =
            po::testing::wide_hessian(data, budget, po::testing::widened(weights));
        std::vector<double> hessian = po::riskparity::linearize(data, budget, weights, true).hessian;
        double scale = po::testing::largest_magnitude(hessian);
        double worst = 0.0;
        for (std::size_t k = 0; k < n * n; ++k)
        {
            worst = std::fmax(worst, std::fabs(static_cast<double>(differences[k]) - hessian[k]) / scale);
        }
        RecordProperty("worst_" + std::to_string(seed), ::testing::PrintToString(worst));
        EXPECT_LE(worst, 1e-3);
    }
}

// A solution with a target volatility and no box at cond(Sigma) = 1e8: its stationarity is above the band of
// 1e-7 and within the rounding floor, where one ulp of the weights moves it by more than the band.
struct FloorLimitedPoint
{
    po::testing::GeneratedRiskParity problem;
    Result result;
};

constexpr double kFloorLimitedBand = 1e-7;

FloorLimitedPoint floor_limited_point()
{
    using Shape = po::testing::RiskParityShape;
    po::testing::GeneratedRiskParity problem = po::testing::generate(
        Shape{10, 1e8, 1.0, Shape::Budget::Uneven, Shape::Box::None, false, true, 7929});
    Result result = po::solve(problem.spec, problem.data).value();
    return FloorLimitedPoint{std::move(problem), std::move(result)};
}

constexpr double kSixtyFourUlps = 64.0 * std::numeric_limits<double>::epsilon();

// Every weight moved by a relative amount, in alternating directions.
void move_weights(std::vector<double>& weights, double relative)
{
    double direction = 1.0;
    for (double& weight : weights)
    {
        weight *= 1.0 + direction * relative;
        direction = -direction;
    }
}

TEST(RiskParityReference, KktAcceptsAStationarityAtItsRoundingFloor)
{
    FloorLimitedPoint point = floor_limited_point();
    ASSERT_EQ(point.result.status, SolveStatus::Optimal);
    po::testing::RiskParityStationarity stationarity = po::testing::risk_parity_stationarity(
        point.problem.spec, point.problem.data, point.result, kFloorLimitedBand);
    EXPECT_GT(stationarity.relative, kFloorLimitedBand);
    EXPECT_LE(stationarity.over_floor, po::testing::kRoundingFloorFactor);
    const Spec& spec = point.problem.spec;
    po::testing::expect_risk_parity_kkt(spec, point.problem.data, point.result, kFloorLimitedBand);
}

// The floor is not a wider band: 64 ulps of every weight are above it, while the equalities still hold at the
// band.
TEST(RiskParityReference, KktRejectsAStationarityAboveItsRoundingFloor)
{
    if (!po::testing::kLongDoubleIsExtended)
    {
        GTEST_SKIP() << "long double is double here: the oracle has no extra precision";
    }
    FloorLimitedPoint point = floor_limited_point();
    ASSERT_EQ(point.result.status, SolveStatus::Optimal);
    move_weights(point.result.weights, kSixtyFourUlps);
    std::vector<std::string> failures = po::testing::risk_parity_kkt_failures(
        point.problem.spec, point.problem.data, point.result, kFloorLimitedBand);
    ASSERT_EQ(failures.size(), 1u);
    EXPECT_NE(failures.front().find("of its rounding floor"), std::string::npos) << failures.front();
    po::testing::RiskParityStationarity stationarity = po::testing::risk_parity_stationarity(
        point.problem.spec, point.problem.data, point.result, kFloorLimitedBand);
    EXPECT_GT(stationarity.over_floor, 8.0);
}

// Where long double is double (MSVC) a result above the band cannot be checked against its floor: the check
// skips the running test, and risk_parity_kkt_failures passes that skip on instead of counting it as a
// failure. Forced here with extended = false, the skip is caught so that this test itself runs on.
TEST(RiskParityReference, FloorCheckWithoutExtendedPrecisionSkipsAndDoesNotFail)
{
    FloorLimitedPoint point = floor_limited_point();
    ASSERT_EQ(point.result.status, SolveStatus::Optimal);
    ::testing::TestPartResultArray parts;
    std::vector<std::string> failures;
    {
        ::testing::ScopedFakeTestPartResultReporter reporter(
            ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD, &parts);
        failures = po::testing::risk_parity_kkt_failures(point.problem.spec, point.problem.data, point.result,
                                                         kFloorLimitedBand, false);
    }
    EXPECT_TRUE(failures.empty());
    ASSERT_EQ(parts.size(), 1);
    EXPECT_TRUE(parts.GetTestPartResult(0).skipped());
    std::vector<std::string> extended = po::testing::risk_parity_kkt_failures(
        point.problem.spec, point.problem.data, point.result, kFloorLimitedBand, true);
    EXPECT_TRUE(extended.empty());
}

// The stationarity at its rounding floor bounds the backward error, not the distance to the solution, which
// is larger by the conditioning of the KKT system: the weights of the result are 1e-12 from the KKT point
// refined in long double, far inside the band. Weights moved by 1e-10 are measured as 1e-10 away.
TEST(RiskParityReference, ForwardErrorMeasuresTheDistanceToTheRefinedPoint)
{
    if (!po::testing::kLongDoubleIsExtended)
    {
        GTEST_SKIP() << "long double is double here: the refinement has no extra precision";
    }
    FloorLimitedPoint point = floor_limited_point();
    ASSERT_EQ(point.result.status, SolveStatus::Optimal);
    po::testing::ForwardError solved =
        po::testing::forward_error(point.problem.spec, point.problem.data, point.result);
    EXPECT_LE(solved.reference_over_floor, 1e-2);
    EXPECT_LE(solved.weights, 1e-11);
    Result moved = point.result;
    move_weights(moved.weights, 1e-10);
    po::testing::ForwardError error =
        po::testing::forward_error(point.problem.spec, point.problem.data, moved);
    EXPECT_LE(error.reference_over_floor, 1e-2);
    EXPECT_GE(error.weights, 0.5e-10);
    EXPECT_LE(error.weights, 2e-10);
}

// The KKT point refined from a result and the result with its equality multipliers, the budget multiplier
// then moved by a fraction of the largest of them.
struct NearRefined
{
    std::vector<long double> refined;
    Result result;
};

NearRefined near_refined(const Spec& spec, const Data& data, const Result& result, double fraction)
{
    NearRefined near{po::testing::refined_multipliers(spec, data, result).value(), result};
    po::testing::set_equality_multipliers(spec, near.refined, near.result);
    double largest =
        std::fmax(std::fabs(near.result.budget_multiplier), std::fabs(near.result.volatility_multiplier));
    near.result.budget_multiplier += fraction * largest;
    return near;
}

// The stationarity at the weights of a result is checked with the least-squares fit there, as the KKT gate of
// the library does, and the multipliers of the result against the KKT point refined from it. At
// s = sigma_min (1 + 1e-10) (n = 5, cond(Sigma) = 1e4, equal budget, no box) the gradients of the budget and
// of the volatility are nearly parallel and the fit is about 5e-8 from the refined multipliers, along the
// direction the stationarity does not see. Multipliers half as far from the refined ones, all of it in the
// budget multiplier, leave a stationarity far above the band and the floor, yet are closer than the fit and
// pass. At s = 1.2 sigma_min of another shape, where the fit is right to about 1e-11, multipliers 1e-6 from
// the refined ones fail, by the comparison with the reference alone.
TEST(RiskParityReference, KktUsesLeastSquaresMultipliersForTheFloor)
{
    if (!po::testing::kLongDoubleIsExtended)
    {
        GTEST_SKIP() << "long double is double here: the refinement has no extra precision";
    }
    std::size_t n = 5;
    Data data = covariance_data(n, random_vector(n, 513, 0.03, 0.15), random_covariance(n, 13, 1e4));
    Spec spec;
    spec.objective = po::Objective::RiskParity;
    spec.risk_budget.assign(n, 1.0 / static_cast<double>(n));
    Result minimum = po::testing::frontier_point(spec, data, std::nullopt);
    spec.target_volatility = po::testing::portfolio_volatility(data, minimum.weights) * (1.0 + 1e-10);
    Result solved = po::solve(spec, data).value();
    ASSERT_EQ(solved.status, SolveStatus::Optimal);
    Result fit = po::testing::least_squares_fit(spec, data, solved).value();
    EXPECT_LE(po::testing::risk_parity_stationarity(spec, data, fit, 1e-10).over_floor,
              po::testing::kRoundingFloorFactor);
    std::vector<long double> refined = po::testing::refined_multipliers(spec, data, solved).value();
    std::vector<long double> fitted = po::testing::wide_kkt_of(spec, fit).multipliers;
    double fit_distance = po::testing::multiplier_distance(fitted, refined);
    EXPECT_GT(fit_distance, 1e-8);
    NearRefined close = near_refined(spec, data, solved, 0.5 * fit_distance);
    po::testing::RiskParityStationarity own =
        po::testing::risk_parity_stationarity(spec, data, close.result, 1e-10);
    RecordProperty("fit_distance", ::testing::PrintToString(fit_distance));
    RecordProperty("own_over_floor", ::testing::PrintToString(own.over_floor));
    EXPECT_GT(own.relative, 1e-10);
    EXPECT_GT(own.over_floor, po::testing::kRoundingFloorFactor);
    std::vector<std::string> failures =
        po::testing::risk_parity_kkt_failures(spec, data, close.result, 1e-10);
    EXPECT_TRUE(failures.empty()) << failures.front();

    using Shape = po::testing::RiskParityShape;
    po::testing::GeneratedRiskParity plain = po::testing::generate(
        Shape{5, 1e4, 1.0, Shape::Budget::Equal, Shape::Box::None, false, true, 402});
    Result result = po::solve(plain.spec, plain.data).value();
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    po::testing::expect_risk_parity_kkt(plain.spec, plain.data, result, 1e-10);
    NearRefined away = near_refined(plain.spec, plain.data, result, 1e-6);
    failures = po::testing::risk_parity_kkt_failures(plain.spec, plain.data, away.result, 1e-10);
    ASSERT_EQ(failures.size(), 1u);
    EXPECT_NE(failures.front().find("refined KKT point"), std::string::npos) << failures.front();
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
