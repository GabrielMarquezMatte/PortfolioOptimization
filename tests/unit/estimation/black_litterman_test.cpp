#include "estimation/black_litterman.hpp"
#include "estimation/posterior.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

namespace
{

using po::Data;
using po::estimation::equilibrium;
using po::estimation::omega_he_litterman;
using po::estimation::omega_idzorek;
using po::estimation::posterior;
using po::estimation::Views;

const std::vector<double> kCovariance = {0.04, 0.01, 0.01, 0.09};

Data make_data(std::size_t n, const std::vector<double>& covariance)
{
    auto data = Data::create(n, {}, covariance);
    EXPECT_TRUE(data.has_value());
    return std::move(*data);
}

Views make_views(std::size_t n, std::vector<double> p)
{
    auto views = Views::create(1, n, p, std::vector<double>{0.02});
    EXPECT_TRUE(views.has_value());
    return std::move(*views);
}

template <typename T>
void expect_failure(const std::expected<T, po::Failure>& result, const char* message)
{
    ASSERT_FALSE(result.has_value());
    EXPECT_EQ(result.error().kind, po::FailureKind::InvalidArgument);
    EXPECT_EQ(result.error().message, message);
}

TEST(Equilibrium, ByHand)
{
    auto pi = equilibrium(make_data(2, kCovariance), std::vector<double>{0.6, 0.4}, 2.5);
    ASSERT_TRUE(pi.has_value());
    EXPECT_NEAR((*pi)[0], 0.07, 1e-16);
    EXPECT_NEAR((*pi)[1], 0.105, 1e-16);
}

TEST(Equilibrium, FloatingPointSumIsAccepted)
{
    std::vector<double> diagonal = {0.04, 0, 0, 0, 0.09, 0, 0, 0, 0.16};
    double third = 1.0 / 3.0;
    EXPECT_TRUE(equilibrium(make_data(3, diagonal), std::vector<double>{third, third, third}, 2.5).has_value());
}

TEST(Equilibrium, WeightsOffByOnePercentAreRejected)
{
    expect_failure(equilibrium(make_data(2, kCovariance), std::vector<double>{0.6, 0.39}, 2.5),
                   "w_mkt must sum to 1 within 1e-8");
}

TEST(Equilibrium, DeltaMustBePositive)
{
    Data data = make_data(2, kCovariance);
    std::vector<double> w = {0.6, 0.4};
    expect_failure(equilibrium(data, w, 0.0), "delta must be positive and finite");
    expect_failure(equilibrium(data, w, std::nan("")), "delta must be positive and finite");
}

TEST(Equilibrium, InputsAreChecked)
{
    Data data = make_data(2, kCovariance);
    expect_failure(equilibrium(data, std::vector<double>{1.0}, 2.5), "w_mkt has 1 entries; expected 2");
    expect_failure(equilibrium(data, std::vector<double>{std::nan(""), 0.4}, 2.5), "w_mkt is not finite at 0");
    auto no_covariance = Data::create(2, std::vector<double>{0.01, 0.02}, {});
    ASSERT_TRUE(no_covariance.has_value());
    expect_failure(equilibrium(*no_covariance, std::vector<double>{0.6, 0.4}, 2.5), "data has no covariance");
}

TEST(Omega, HeLitterman)
{
    auto omega = omega_he_litterman(make_data(2, kCovariance), make_views(2, {1.0, -1.0}), 0.05);
    ASSERT_TRUE(omega.has_value());
    ASSERT_EQ(omega->size(), 1u);
    EXPECT_NEAR((*omega)[0], 0.0055, 1e-18);
}

TEST(Omega, HeLittermanIsDiagonal)
{
    auto views = Views::create(2, 2, std::vector<double>{1, 0, 0, 1}, std::vector<double>{0.01, 0.02});
    ASSERT_TRUE(views.has_value());
    auto omega = omega_he_litterman(make_data(2, kCovariance), *views, 0.05);
    ASSERT_TRUE(omega.has_value());
    EXPECT_NEAR((*omega)[0], 0.05 * 0.04, 1e-18);
    EXPECT_EQ((*omega)[1], 0.0);
    EXPECT_EQ((*omega)[2], 0.0);
    EXPECT_NEAR((*omega)[3], 0.05 * 0.09, 1e-18);
}

TEST(Omega, Idzorek)
{
    Data data = make_data(2, kCovariance);
    Views views = make_views(2, {1.0, -1.0});
    auto half = omega_idzorek(data, views, 0.05, std::vector<double>{0.5});
    ASSERT_TRUE(half.has_value());
    EXPECT_NEAR((*half)[0], 0.0055, 1e-18);
    auto full = omega_idzorek(data, views, 0.05, std::vector<double>{1.0});
    ASSERT_TRUE(full.has_value());
    EXPECT_EQ((*full)[0], 0.0);
    auto quarter = omega_idzorek(data, views, 0.05, std::vector<double>{0.25});
    ASSERT_TRUE(quarter.has_value());
    EXPECT_NEAR((*quarter)[0], 0.0165, 1e-18);
}

TEST(Omega, IdzorekConfidenceIsChecked)
{
    Data data = make_data(2, kCovariance);
    Views views = make_views(2, {1.0, -1.0});
    expect_failure(omega_idzorek(data, views, 0.05, std::vector<double>{0.0}), "confidence must be in (0, 1] at 0");
    expect_failure(omega_idzorek(data, views, 0.05, std::vector<double>{1.5}), "confidence must be in (0, 1] at 0");
    expect_failure(omega_idzorek(data, views, 0.05, std::vector<double>{std::nan("")}),
                   "confidence must be in (0, 1] at 0");
    expect_failure(omega_idzorek(data, views, 0.05, std::vector<double>{0.5, 0.5}),
                   "confidence has 2 entries; expected 1");
}

TEST(Omega, AssetCountsMustAgree)
{
    expect_failure(omega_he_litterman(make_data(2, kCovariance), make_views(3, {1.0, -1.0, 0.0}), 0.05),
                   "views have 3 assets; data has 2");
}

TEST(Omega, TauMustBePositive)
{
    Data data = make_data(2, kCovariance);
    Views views = make_views(2, {1.0, -1.0});
    expect_failure(omega_he_litterman(data, views, 0.0), "tau must be positive and finite");
    expect_failure(omega_idzorek(data, views, std::nan(""), std::vector<double>{0.5}),
                   "tau must be positive and finite");
}

TEST(Omega, NeedsACovariance)
{
    Views views = make_views(2, {1.0, -1.0});
    auto no_covariance = Data::create(2, std::vector<double>{0.01, 0.02}, {});
    ASSERT_TRUE(no_covariance.has_value());
    expect_failure(omega_he_litterman(*no_covariance, views, 0.05), "data has no covariance");
}

const std::vector<double> kPi = {0.07, 0.105};

Data make_prior(std::size_t n, const std::vector<double>& pi, const std::vector<double>& covariance)
{
    auto data = Data::create(n, pi, covariance);
    EXPECT_TRUE(data.has_value());
    return std::move(*data);
}

Views make_views_with_omega(std::size_t n, const std::vector<double>& p, const std::vector<double>& q,
                            const std::vector<double>& omega)
{
    auto views = Views::create(q.size(), n, p, q);
    EXPECT_TRUE(views.has_value());
    EXPECT_TRUE(views->set_omega(omega).has_value());
    return std::move(*views);
}

void swap_rows(std::vector<double>& a, std::size_t n, std::size_t first, std::size_t second)
{
    for (std::size_t j = 0; j < n; ++j)
    {
        std::swap(a[first * n + j], a[second * n + j]);
    }
}

// Inverse of an n x n row-major matrix by Gauss-Jordan elimination with partial pivoting.
std::vector<double> gauss_jordan_inverse(std::vector<double> a, std::size_t n)
{
    std::vector<double> inverse(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        inverse[i * n + i] = 1.0;
    }
    for (std::size_t c = 0; c < n; ++c)
    {
        std::size_t best = c;
        for (std::size_t r = c + 1; r < n; ++r)
        {
            best = std::fabs(a[r * n + c]) > std::fabs(a[best * n + c]) ? r : best;
        }
        swap_rows(a, n, c, best);
        swap_rows(inverse, n, c, best);
        double pivot = a[c * n + c];
        for (std::size_t j = 0; j < n; ++j)
        {
            a[c * n + j] /= pivot;
            inverse[c * n + j] /= pivot;
        }
        for (std::size_t r = 0; r < n; ++r)
        {
            if (r == c)
            {
                continue;
            }
            double factor = a[r * n + c];
            for (std::size_t j = 0; j < n; ++j)
            {
                a[r * n + j] -= factor * a[c * n + j];
                inverse[r * n + j] -= factor * inverse[c * n + j];
            }
        }
    }
    return inverse;
}

void expect_symmetric(std::span<const double> matrix, std::size_t n)
{
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t j = 0; j < n; ++j)
        {
            EXPECT_EQ(matrix[i * n + j], matrix[j * n + i]);
        }
    }
}

TEST(BlackLitterman, HeLittermanByHand)
{
    auto result = posterior(make_prior(2, kPi, kCovariance), make_views_with_omega(2, {1.0, -1.0}, {0.02}, {0.0055}),
                            0.05, true);
    ASSERT_TRUE(result.has_value());
    ASSERT_EQ(result->mu.size(), 2u);
    EXPECT_NEAR(result->mu[0], 0.0775, 1e-16);
    EXPECT_NEAR(result->mu[1], 0.085, 1e-16);
    std::vector<double> b = {0.0015, -0.004};
    ASSERT_EQ(result->covariance.size(), 4u);
    for (std::size_t i = 0; i < 2; ++i)
    {
        for (std::size_t j = 0; j < 2; ++j)
        {
            double expected = 1.05 * kCovariance[i * 2 + j] - b[i] * b[j] / 0.011;
            EXPECT_NEAR(result->covariance[i * 2 + j], expected, 1e-16);
        }
    }
    expect_symmetric(result->covariance, 2);
}

TEST(BlackLitterman, NoViewsReturnsThePrior)
{
    auto views = Views::create(0, 2, {}, {});
    ASSERT_TRUE(views.has_value());
    auto result = posterior(make_prior(2, kPi, kCovariance), *views, 0.05, true);
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(result->mu, kPi);
    ASSERT_EQ(result->covariance.size(), 4u);
    for (std::size_t i = 0; i < 4; ++i)
    {
        EXPECT_EQ(result->covariance[i], 1.05 * kCovariance[i]);
    }
}

TEST(BlackLitterman, FullConfidenceReproducesTheView)
{
    auto two = posterior(make_prior(2, kPi, kCovariance), make_views_with_omega(2, {1.0, -1.0}, {0.02}, {0.0}),
                         0.05, true);
    ASSERT_TRUE(two.has_value());
    EXPECT_LE(std::fabs(two->mu[0] - two->mu[1] - 0.02), 1e-15);

    std::vector<double> covariance = {0.04,  0.006, 0.002,  0.001, 0.006, 0.09,  0.009, 0.003,
                                      0.002, 0.009, 0.0625, 0.004, 0.001, 0.003, 0.004, 0.0484};
    std::vector<double> p = {1.0, -1.0, 0.0, 0.0, 0.5, 0.5, -1.0, 0.0};
    std::vector<double> q = {0.02, 0.01};
    auto four = posterior(make_prior(4, {0.05, 0.07, 0.06, 0.04}, covariance),
                          make_views_with_omega(4, p, q, {0.0, 0.0, 0.0, 0.0}), 0.05, true);
    ASSERT_TRUE(four.has_value());
    for (std::size_t r = 0; r < 2; ++r)
    {
        double view = 0.0;
        for (std::size_t i = 0; i < 4; ++i)
        {
            view += p[r * 4 + i] * four->mu[i];
        }
        EXPECT_LE(std::fabs(view - q[r]), 1e-15) << r;
    }
}

// P' Omega^-1 P (n x n) for k views.
std::vector<double> view_precision(const std::vector<double>& p, const std::vector<double>& omega, std::size_t n)
{
    std::size_t k = p.size() / n;
    std::vector<double> omega_inverse = gauss_jordan_inverse(omega, k);
    std::vector<double> precision(n * n, 0.0);
    for (std::size_t i = 0; i < n * n; ++i)
    {
        for (std::size_t r = 0; r < k; ++r)
        {
            for (std::size_t s = 0; s < k; ++s)
            {
                precision[i] += p[r * n + i / n] * omega_inverse[r * k + s] * p[s * n + i % n];
            }
        }
    }
    return precision;
}

// Sigma + ((tau Sigma)^-1 + P' Omega^-1 P)^-1, the textbook form, which needs Omega > 0.
std::vector<double> inverse_form(const std::vector<double>& sigma, const std::vector<double>& p,
                                 const std::vector<double>& omega, double tau, std::size_t n)
{
    std::vector<double> scaled(n * n);
    for (std::size_t i = 0; i < n * n; ++i)
    {
        scaled[i] = tau * sigma[i];
    }
    std::vector<double> precision = gauss_jordan_inverse(scaled, n);
    std::vector<double> views = view_precision(p, omega, n);
    for (std::size_t i = 0; i < n * n; ++i)
    {
        precision[i] += views[i];
    }
    std::vector<double> expected = gauss_jordan_inverse(precision, n);
    for (std::size_t i = 0; i < n * n; ++i)
    {
        expected[i] += sigma[i];
    }
    return expected;
}

TEST(BlackLitterman, PosteriorCovarianceMatchesTheInverseForm)
{
    std::vector<double> sigma = {0.04, 0.006, 0.002, 0.006, 0.09, 0.009, 0.002, 0.009, 0.0625};
    std::vector<double> p = {1.0, 0.0, 0.0, 0.0, 1.0, -1.0};
    std::vector<double> omega = {0.001, 0.0, 0.0, 0.002};
    auto result = posterior(make_prior(3, {0.05, 0.07, 0.06}, sigma),
                            make_views_with_omega(3, p, {0.06, 0.01}, omega), 0.05, true);
    ASSERT_TRUE(result.has_value());
    std::vector<double> expected = inverse_form(sigma, p, omega, 0.05, 3);
    ASSERT_EQ(result->covariance.size(), 9u);
    for (std::size_t i = 0; i < 9; ++i)
    {
        EXPECT_LE(std::fabs(result->covariance[i] - expected[i]), 1e-14 * std::fabs(expected[i])) << i;
    }
    expect_symmetric(result->covariance, 3);
}

TEST(BlackLitterman, OmegaIsRequired)
{
    expect_failure(posterior(make_prior(2, kPi, kCovariance), make_views(2, {1.0, -1.0}), 0.05, true),
                   "omega is not set");
    expect_failure(posterior(make_data(2, kCovariance), make_views(2, {1.0, -1.0}), 0.05, true), "data has no mu");
}

TEST(BlackLitterman, NearlyRedundantViewsAreFiniteOrRejected)
{
    auto result = posterior(make_prior(2, kPi, kCovariance),
                            make_views_with_omega(2, {1.0, -1.0, 1.0, -1.0 + 1e-9}, {0.02, 0.02}, {0, 0, 0, 0}),
                            0.05, true);
    if (!result.has_value())
    {
        expect_failure(result, "P tau Sigma P' + omega is not positive definite: pivot 1 is not positive");
        return;
    }
    for (double value : result->mu)
    {
        EXPECT_TRUE(std::isfinite(value));
    }
    for (double value : result->covariance)
    {
        EXPECT_TRUE(std::isfinite(value));
    }
}

TEST(BlackLitterman, WithoutCovarianceLeavesItEmpty)
{
    Data data = make_prior(2, kPi, kCovariance);
    Views views = make_views_with_omega(2, {1.0, -1.0}, {0.02}, {0.0055});
    auto with = posterior(data, views, 0.05, true);
    auto without = posterior(data, views, 0.05, false);
    ASSERT_TRUE(with.has_value());
    ASSERT_TRUE(without.has_value());
    EXPECT_TRUE(without->covariance.empty());
    EXPECT_EQ(without->mu, with->mu);
}

} // namespace
