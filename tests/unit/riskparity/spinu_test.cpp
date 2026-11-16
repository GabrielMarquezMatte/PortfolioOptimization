#include "riskparity/spinu.hpp"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <span>
#include <vector>

#include "core/data.hpp"
#include "fuzz/risk_parity/risk_parity_reference.hpp"
#include "metrics/portfolio.hpp"
#include "random_covariance.hpp"

namespace
{

using po::Data;
using po::riskparity::SpinuResult;
using po::testing::random_covariance;
using po::testing::random_vector;
using po::testing::risk_contributions;
using po::testing::two_asset_risk_parity;

Data covariance_data(std::size_t n, const std::vector<double>& covariance)
{
    return Data::create(n, {}, covariance).value();
}

double band_of(double condition)
{
    if (condition <= 1e4)
    {
        return 1e-10;
    }
    if (condition <= 1e6)
    {
        return 1e-9;
    }
    return 1e-7;
}

std::vector<double> equal_budget(std::size_t n)
{
    return std::vector<double>(n, 1.0 / static_cast<double>(n));
}

std::vector<double> uneven_budget(std::size_t n, std::uint64_t seed)
{
    std::vector<double> exponents = random_vector(n, seed + 700, -6.0, 0.0);
    std::vector<double> budget(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        budget[i] = std::pow(10.0, exponents[i]);
    }
    double total = std::accumulate(budget.begin(), budget.end(), 0.0);
    for (double& entry : budget)
    {
        entry /= total;
    }
    return budget;
}

double largest_gap(const Data& data, const std::vector<double>& weights, const std::vector<double>& budget)
{
    std::vector<double> contributions = risk_contributions(data, weights);
    double largest = 0.0;
    for (std::size_t i = 0; i < budget.size(); ++i)
    {
        largest = std::fmax(largest, std::fabs(contributions[i] - budget[i]));
    }
    return largest;
}

using Wide = long double;

// Gaussian elimination with partial pivoting in long double; the right-hand side is overwritten by the solution.
void eliminate(std::vector<Wide>& matrix, std::vector<Wide>& rhs, std::size_t n)
{
    for (std::size_t k = 0; k < n; ++k)
    {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i)
        {
            if (std::fabs(matrix[i * n + k]) > std::fabs(matrix[pivot * n + k]))
            {
                pivot = i;
            }
        }
        for (std::size_t j = 0; j < n; ++j)
        {
            std::swap(matrix[k * n + j], matrix[pivot * n + j]);
        }
        std::swap(rhs[k], rhs[pivot]);
        for (std::size_t i = k + 1; i < n; ++i)
        {
            Wide factor = matrix[i * n + k] / matrix[k * n + k];
            for (std::size_t j = k; j < n; ++j)
            {
                matrix[i * n + j] -= factor * matrix[k * n + j];
            }
            rhs[i] -= factor * rhs[k];
        }
    }
    for (std::size_t k = n; k-- > 0;)
    {
        for (std::size_t j = k + 1; j < n; ++j)
        {
            rhs[k] -= matrix[k * n + j] * rhs[j];
        }
        rhs[k] /= matrix[k * n + k];
    }
}

// One Newton step of F(y) = 1/2 y'Sigma y - b' log y in long double, at most 0.9 of the way to the boundary; the
// largest |step_i / y_i| it took.
Wide wide_newton_step(const std::vector<double>& covariance, const std::vector<double>& budget,
                      std::vector<Wide>& y)
{
    std::size_t n = y.size();
    std::vector<Wide> hessian(covariance.begin(), covariance.end());
    std::vector<Wide> step(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        std::span<const double> row(covariance.data() + i * n, n);
        Wide product = po::testing::compensated_dot_of<Wide>(n, row, y);
        step[i] = budget[i] / y[i] - product;
        hessian[i * n + i] += budget[i] / (y[i] * y[i]);
    }
    eliminate(hessian, step, n);
    Wide length = 1.0L;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (step[i] < 0.0L)
        {
            length = std::fmin(length, 0.9L * y[i] / -step[i]);
        }
    }
    Wide relative = 0.0L;
    for (std::size_t i = 0; i < n; ++i)
    {
        y[i] += length * step[i];
        relative = std::fmax(relative, std::fabs(step[i] / y[i]));
    }
    return relative;
}

// Below 1e-17 of the weights, or 8 ulps where long double is double (MSVC): the residual b / y - Sigma y is a
// compensated dot product, so the iteration converges to the rounding of Wide either way.
constexpr Wide kWideStop =
    std::numeric_limits<Wide>::digits > 53 ? 1e-17L : 8.0L * std::numeric_limits<Wide>::epsilon();

// The weights of Spinu by Newton in long double until the step is below kWideStop of the weights: a reference
// for the relative error of each weight, independent of the code under test.
std::vector<Wide> wide_spinu_weights(const std::vector<double>& covariance, const std::vector<double>& budget)
{
    std::size_t n = budget.size();
    std::vector<Wide> y(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        y[i] = std::sqrt(static_cast<Wide>(budget[i]) / covariance[i * n + i]);
    }
    int iteration = 0;
    while (iteration < 400 && wide_newton_step(covariance, budget, y) > kWideStop)
    {
        ++iteration;
    }
    Wide total = std::accumulate(y.begin(), y.end(), 0.0L);
    for (Wide& entry : y)
    {
        entry /= total;
    }
    return y;
}

Wide largest_relative_error(const std::vector<double>& weights, const std::vector<Wide>& reference)
{
    Wide largest = 0.0L;
    for (std::size_t i = 0; i < weights.size(); ++i)
    {
        largest = std::fmax(largest, std::fabs(weights[i] - reference[i]) / reference[i]);
    }
    return largest;
}

void expect_diagonal_closed_form(const std::vector<double>& sigma, const std::vector<double>& budget)
{
    std::size_t n = sigma.size();
    std::vector<double> covariance(n * n, 0.0);
    std::vector<double> expected(n);
    double total = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        covariance[i * n + i] = sigma[i] * sigma[i];
        expected[i] = std::sqrt(budget[i]) / sigma[i];
        total += expected[i];
    }
    SpinuResult result = po::riskparity::spinu(covariance_data(n, covariance), budget);
    ASSERT_TRUE(result.converged);
    for (std::size_t i = 0; i < n; ++i)
    {
        double weight = expected[i] / total;
        EXPECT_NEAR(result.weights[i], weight, 1e-12 * weight) << "asset " << i;
    }
}

TEST(Spinu, DiagonalClosedForm)
{
    std::vector<double> sigma = {0.1, 0.2, 0.35, 0.05, 0.4};
    expect_diagonal_closed_form(sigma, equal_budget(5));
    expect_diagonal_closed_form(sigma, {0.1, 0.4, 0.2, 0.05, 0.25});
    expect_diagonal_closed_form({0.01, 0.4, 0.1, 0.002, 0.3}, {0.9, 1e-6, 0.05, 0.049998, 1e-6});
}

TEST(Spinu, TwoAssetsClosedForm)
{
    for (double rho : {-0.5, 0.0, 0.3, 0.9})
    {
        for (double budget : {0.5, 0.2, 0.9})
        {
            SCOPED_TRACE(testing::Message() << "rho " << rho << " b1 " << budget);
            double sigma1 = 0.15;
            double sigma2 = 0.3;
            double cross = rho * sigma1 * sigma2;
            std::vector<double> covariance = {sigma1 * sigma1, cross, cross, sigma2 * sigma2};
            std::vector<double> budgets = {budget, 1.0 - budget};
            SpinuResult result = po::riskparity::spinu(covariance_data(2, covariance), budgets);
            std::vector<double> expected = two_asset_risk_parity(sigma1, sigma2, rho, budget);
            ASSERT_TRUE(result.converged);
            EXPECT_NEAR(result.weights[0], expected[0], 1e-12);
            EXPECT_NEAR(result.weights[1], expected[1], 1e-12);
        }
    }
}

TEST(Spinu, ContributionsEqualBudget)
{
    for (std::size_t n : {1u, 2u, 10u, 100u})
    {
        for (double condition : {1e2, 1e4, 1e6})
        {
            for (std::uint64_t seed = 1; seed <= 20; ++seed)
            {
                SCOPED_TRACE(testing::Message() << "n " << n << " cond " << condition << " seed " << seed);
                Data data = covariance_data(n, random_covariance(n, seed, condition));
                std::vector<double> budget = equal_budget(n);
                SpinuResult result = po::riskparity::spinu(data, budget);
                ASSERT_TRUE(result.converged);
                double band = band_of(condition);
                EXPECT_LE(largest_gap(data, result.weights, budget), band);
                EXPECT_LE(po::metrics::risk_budget_residual(data, result.weights, budget),
                          static_cast<double>(n) * band * band);
            }
        }
    }
}

TEST(Spinu, UnevenBudgetsConditioningSweep)
{
    for (std::size_t n : {5u, 20u})
    {
        for (double condition : {1e2, 1e4, 1e6, 1e8})
        {
            for (double scale : {1.0, 1e-4})
            {
                for (std::uint64_t seed = 1; seed <= 5; ++seed)
                {
                    SCOPED_TRACE(testing::Message()
                                 << "n " << n << " cond " << condition << " scale " << scale << " seed " << seed);
                    std::vector<double> covariance = random_covariance(n, seed, condition);
                    for (double& entry : covariance)
                    {
                        entry *= scale;
                    }
                    Data data = covariance_data(n, covariance);
                    std::vector<double> budget = uneven_budget(n, seed);
                    SpinuResult result = po::riskparity::spinu(data, budget);
                    ASSERT_TRUE(result.converged);
                    EXPECT_TRUE(std::ranges::all_of(result.weights, [](double weight) { return weight > 0.0; }));
                    EXPECT_LE(largest_gap(data, result.weights, budget), band_of(condition));
                    EXPECT_LE(largest_relative_error(result.weights, wide_spinu_weights(covariance, budget)),
                              band_of(condition));
                }
            }
        }
    }
}

} // namespace
