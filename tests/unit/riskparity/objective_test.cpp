#include "riskparity/objective.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <random>
#include <span>
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

// The objective sums Sigma w in long double and the metric in double, so they agree to the rounding of the
// metric.
TEST(RiskParityObjective, ValueMatchesMetric)
{
    for (std::size_t n : {2u, 7u, 30u})
    {
        SCOPED_TRACE(testing::Message() << "n " << n);
        Data data = sample_data(n, 3);
        std::vector<double> budget = sample_budget(n, 3);
        std::vector<double> weights = sample_weights(n, 3);
        double metric = po::metrics::risk_budget_residual(data, weights, budget);
        EXPECT_NEAR(linearize(data, budget, weights, false).value, metric, 1e-14 * metric);
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
    for (double length : {1e-8, 1e-10, 1e-12, 1e-15})
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
        EXPECT_LE(std::fabs(change.value - expansion), 1e-12 * std::fabs(slope));
        double variance_change = 2.0 * po::testing::dot(po::testing::covariance_times(data, weights), step) +
                                 po::testing::dot(po::testing::covariance_times(data, step), step);
        EXPECT_NEAR(change.variance, variance_change, 1e-14 * std::fabs(variance_change));
    }
}

// 2 J'J as it was computed before the symmetry: every row i of J added to all n^2 entries, each entry summed
// from +0 in increasing i.
std::vector<double> reference_gauss_newton(std::span<const double> jacobian, std::size_t n)
{
    std::vector<double> curvature(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            double scaled = 2.0 * jacobian[i * n + j];
            for (std::size_t k = 0; k < n; ++k)
            {
                curvature[j * n + k] += scaled * jacobian[i * n + k];
            }
        }
    }
    return curvature;
}

bool same_bits(const std::vector<double>& a, const std::vector<double>& b)
{
    return a.size() == b.size() && std::memcmp(a.data(), b.data(), a.size() * sizeof(double)) == 0;
}

TEST(GaussNewton, SymmetricIsBitIdenticalToTheFullLoop)
{
    for (std::size_t n : {1, 2, 3, 7, 16, 17, 64, 100})
    {
        for (double condition : {1e2, 1e8})
        {
            SCOPED_TRACE(testing::Message() << "n = " << n << ", condition = " << condition);
            std::vector<double> mu = random_vector(n, n + 600, 0.03, 0.15);
            Data data = Data::create(n, mu, random_covariance(n, n + 7, condition)).value();
            std::vector<double> weights = normalized(random_vector(n, n + 700, 0.2, 1.0));
            std::vector<double> budget(n, 1.0 / static_cast<double>(n));
            Linearization linearization = linearize(data, budget, weights, false);
            EXPECT_TRUE(same_bits(po::riskparity::gauss_newton(linearization.jacobian, n),
                                  reference_gauss_newton(linearization.jacobian, n)));
        }
    }
}

// Signed zeros, subnormals and entries up to the largest that doubles exactly; and one entry above it, where
// 2 J_ij overflows and the products of the two halves are no longer the same.
TEST(GaussNewton, SymmetricIsBitIdenticalOnEdgeEntries)
{
    constexpr double kLargest = std::numeric_limits<double>::max();
    constexpr double kSubnormal = std::numeric_limits<double>::denorm_min();
    for (std::size_t n : {1, 2, 3, 7, 16, 17, 64, 100})
    {
        SCOPED_TRACE(testing::Message() << "n = " << n);
        std::vector<double> jacobian = random_vector(n * n, n + 800, -1.0, 1.0);
        std::vector<double> edges = {0.0, -0.0, kSubnormal, -3.0 * kSubnormal, 1e-310, -1e-300, 1e150,
                                     -1e154, kLargest / 2.0, -kLargest / 2.0};
        for (std::size_t e = 0; e < jacobian.size(); e += 3)
        {
            jacobian[e] = edges[(e / 3) % edges.size()];
        }
        std::vector<double> reference = reference_gauss_newton(jacobian, n);
        EXPECT_TRUE(same_bits(po::riskparity::gauss_newton(jacobian, n), reference));
        jacobian[n * n / 2] = 0.75 * kLargest;
        reference = reference_gauss_newton(jacobian, n);
        EXPECT_TRUE(same_bits(po::riskparity::gauss_newton(jacobian, n), reference));
    }
}

// The compensated dot product as it was computed before Dekker's product: TwoProduct by std::fma.
double reference_compensated_dot(std::span<const double> a, std::span<const double> b)
{
    double sum = 0.0;
    double correction = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i)
    {
        double product = a[i] * b[i];
        double product_error = std::fma(a[i], b[i], -product);
        double total = sum + product;
        double carried = total - sum;
        double sum_error = (sum - (total - carried)) + (product - carried);
        sum = total;
        correction += sum_error + product_error;
    }
    return sum + correction;
}

bool same_bits(double a, double b)
{
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

void expect_error_matches_fma(double a, double b)
{
    double product = a * b;
    double error = po::riskparity::two_product_error(a, b, product);
    double reference = std::fma(a, b, -product);
    EXPECT_TRUE(same_bits(error, reference))
        << std::hexfloat << a << " * " << b << ": " << error << " against " << reference;
}

// sign (1 + u) 2^e for u uniform in [0, 1) and e uniform in [low, high]; below 2^-1022 a subnormal.
double random_double(std::mt19937_64& engine, int low, int high)
{
    std::uniform_int_distribution<int> exponent(low, high);
    double sign = (engine() & 1U) == 0U ? 1.0 : -1.0;
    return sign * std::ldexp(1.0 + po::testing::uniform_open(engine), exponent(engine));
}

// Each side of the range of Dekker's product: its ends 2^-485 and 2^500, the doubles next to them, products
// whose error reaches the smallest subnormal, splits that overflow (above 2^996), subnormals, signed zeros,
// infinities; then 2000 random pairs over the exponents of the range's ends and 2000 over every exponent.
TEST(CompensatedDot, DekkerMatchesFmaOnEdgeCases)
{
    constexpr double kLargest = std::numeric_limits<double>::max();
    constexpr double kInfinity = std::numeric_limits<double>::infinity();
    constexpr double kSubnormal = std::numeric_limits<double>::denorm_min();
    std::vector<double> edges = {0.0,
                                 -0.0,
                                 1.0,
                                 -1.0,
                                 0x1p-485,
                                 0x1.fffffffffffffp-486,
                                 0x1.0000000000001p-485,
                                 0x1.fffffffffffffp-485,
                                 0x1p500,
                                 0x1.fffffffffffffp499,
                                 0x1.0000000000001p500,
                                 0x1.fffffffffffffp500,
                                 0x1p996,
                                 0x1p997,
                                 0x1.8p1000,
                                 kLargest,
                                 0x1p-1022,
                                 0x1p-1000,
                                 kSubnormal,
                                 3.0 * kSubnormal,
                                 0x1.8p-1060,
                                 kInfinity,
                                 0x1.5555555555555p0,
                                 0x1.aaaaaaaaaaaabp-300,
                                 0x1.fffffffffffffp200,
                                 134217729.0};
    for (double a : edges)
    {
        for (double b : edges)
        {
            expect_error_matches_fma(a, b);
            expect_error_matches_fma(-a, b);
        }
    }
    std::mt19937_64 engine(2026);
    for (int k = 0; k < 2000; ++k)
    {
        expect_error_matches_fma(random_double(engine, -490, -480), random_double(engine, -490, -480));
        expect_error_matches_fma(random_double(engine, 495, 505), random_double(engine, -10, 10));
        expect_error_matches_fma(random_double(engine, -1074, 1023), random_double(engine, -1074, 1023));
    }
}

TEST(CompensatedDot, IsBitIdenticalToTheFmaDot)
{
    std::mt19937_64 engine(7);
    for (std::size_t n : {1, 2, 3, 7, 16, 17, 64, 100})
    {
        for (double condition : {1e2, 1e8})
        {
            SCOPED_TRACE(testing::Message() << "n = " << n << ", condition = " << condition);
            std::vector<double> covariance = random_covariance(n, n + 9, condition);
            std::vector<double> weights = random_vector(n, n + 10, -1.0, 1.0);
            for (std::size_t i = 0; i < n; ++i)
            {
                std::span<const double> row = std::span<const double>(covariance).subspan(i * n, n);
                EXPECT_TRUE(same_bits(po::riskparity::compensated_dot(row, weights),
                                      reference_compensated_dot(row, weights)));
            }
            std::vector<double> wide(n);
            for (double& entry : wide)
            {
                entry = random_double(engine, -600, 600);
            }
            wide[0] = -0.0;
            EXPECT_TRUE(same_bits(po::riskparity::compensated_dot(wide, weights),
                                  reference_compensated_dot(wide, weights)));
        }
    }
}

} // namespace
