#include <gtest/gtest.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <portopt/portopt.h>
#include <string>
#include <vector>

#include "capi_scope.hpp"
#include "estimation/bayes_stein.hpp"
#include "estimation/black_litterman.hpp"
#include "estimation/moments.hpp"
#include "estimation/posterior.hpp"
#include "estimation/views.hpp"

namespace
{

using po::testing::DataPtr;

struct ViewsFree
{
    void operator()(po_views* handle) const
    {
        po_views_free(handle);
    }
};

using ViewsPtr = std::unique_ptr<po_views, ViewsFree>;

const std::vector<double> kCovariance = {0.04, 0.01, 0.01, 0.09};
const std::vector<double> kMu = {0.07, 0.105};
const std::vector<double> kP = {1.0, -1.0};
const std::vector<double> kQ = {0.02};
const std::vector<double> kReturns = {0.01, 0.02, 0.03, -0.01, -0.02, 0.04, 0.00, 0.01};
constexpr double kTau = 0.05;
constexpr double kPoison = -7.0;
constexpr std::size_t kHuge = std::numeric_limits<std::size_t>::max() / 2 + 1;

ViewsPtr make_views(std::size_t k, std::size_t n, const std::vector<double>& p, const std::vector<double>& q)
{
    po_views* raw = nullptr;
    EXPECT_EQ(po_views_create(k, n, p.data(), q.data(), &raw), PO_OK);
    return ViewsPtr(raw);
}

ViewsPtr make_views_with_omega()
{
    ViewsPtr views = make_views(1, 2, kP, kQ);
    const double omega[1] = {0.001};
    EXPECT_EQ(po_views_set_omega(views.get(), omega, 1), PO_OK);
    return views;
}

po::estimation::Views core_views()
{
    auto views = po::estimation::Views::create(1, 2, kP, kQ);
    EXPECT_TRUE(views.has_value());
    EXPECT_TRUE(views->set_omega(std::vector<double>{0.001}).has_value());
    return std::move(*views);
}

po::Data core_data(const std::vector<double>& mu, const std::vector<double>& covariance)
{
    auto data = po::Data::create(2, mu, covariance);
    EXPECT_TRUE(data.has_value());
    return std::move(*data);
}

TEST(CapiEstimation, EachFunctionMatchesTheCore)
{
    auto moments = po::estimation::sample_moments(kReturns, 4, 2);
    ASSERT_TRUE(moments.has_value());
    double mu[2] = {};
    double cov[4] = {};
    ASSERT_EQ(po_sample_moments(4, 2, kReturns.data(), mu, cov), PO_OK);
    EXPECT_EQ(std::vector<double>(mu, mu + 2), moments->mu);
    EXPECT_EQ(std::vector<double>(cov, cov + 4), moments->covariance);

    // t = 12 > n + 2 with the covariance of the four-period panel above.
    DataPtr sample = po::testing::make_data(2, moments->mu, moments->covariance);
    auto shrunk = po::estimation::bayes_stein(core_data(moments->mu, moments->covariance), 12);
    ASSERT_TRUE(shrunk.has_value());
    double shrunk_mu[2] = {};
    double grand_mean = 0.0;
    double weight = 0.0;
    ASSERT_EQ(po_bayes_stein(sample.get(), 12, shrunk_mu, &grand_mean, &weight), PO_OK);
    EXPECT_EQ(grand_mean, shrunk->grand_mean);
    EXPECT_EQ(weight, shrunk->weight);
    EXPECT_EQ(shrunk_mu[0], shrunk->mu[0]);
    EXPECT_EQ(shrunk_mu[1], shrunk->mu[1]);
}

TEST(CapiEstimation, BlackLittermanMatchesTheCore)
{
    DataPtr prior = po::testing::make_data(2, kMu, kCovariance);
    po::Data data = core_data(kMu, kCovariance);
    const double w_mkt[2] = {0.6, 0.4};
    double pi[2] = {};
    ASSERT_EQ(po_bl_equilibrium(prior.get(), w_mkt, 2, 2.5, pi), PO_OK);
    auto expected_pi = po::estimation::equilibrium(data, w_mkt, 2.5);
    ASSERT_TRUE(expected_pi.has_value());
    EXPECT_EQ(std::vector<double>(pi, pi + 2), *expected_pi);

    ViewsPtr views = make_views_with_omega();
    po::estimation::Views core = core_views();
    double omega[1] = {};
    ASSERT_EQ(po_bl_omega_he_litterman(prior.get(), views.get(), kTau, omega), PO_OK);
    EXPECT_EQ(omega[0], (*po::estimation::omega_he_litterman(data, core, kTau))[0]);
    const double confidence[1] = {0.5};
    ASSERT_EQ(po_bl_omega_idzorek(prior.get(), views.get(), kTau, confidence, omega), PO_OK);
    EXPECT_EQ(omega[0], (*po::estimation::omega_idzorek(data, core, kTau, confidence))[0]);

    double mu[2] = {};
    double cov[4] = {};
    ASSERT_EQ(po_bl_posterior(prior.get(), views.get(), kTau, mu, cov), PO_OK);
    auto expected = po::estimation::posterior(data, core, kTau, true);
    ASSERT_TRUE(expected.has_value());
    EXPECT_EQ(std::vector<double>(mu, mu + 2), expected->mu);
    EXPECT_EQ(std::vector<double>(cov, cov + 4), expected->covariance);
}

TEST(CapiEstimation, PosteriorWithoutCovarianceOutput)
{
    DataPtr prior = po::testing::make_data(2, kMu, kCovariance);
    ViewsPtr views = make_views_with_omega();
    double mu[2] = {};
    ASSERT_EQ(po_bl_posterior(prior.get(), views.get(), kTau, mu, nullptr), PO_OK);
    auto expected = po::estimation::posterior(core_data(kMu, kCovariance), core_views(), kTau, false);
    ASSERT_TRUE(expected.has_value());
    EXPECT_EQ(std::vector<double>(mu, mu + 2), expected->mu);
}

TEST(CapiEstimation, EmptyViewsAcceptNullPointers)
{
    po_views* raw = nullptr;
    ASSERT_EQ(po_views_create(0, 2, nullptr, nullptr, &raw), PO_OK);
    ViewsPtr views(raw);
    DataPtr prior = po::testing::make_data(2, kMu, kCovariance);
    double mu[2] = {};
    double cov[4] = {};
    ASSERT_EQ(po_bl_posterior(prior.get(), views.get(), kTau, mu, cov), PO_OK);
    EXPECT_EQ(mu[0], kMu[0]);
    EXPECT_NEAR(cov[0], (1.0 + kTau) * kCovariance[0], 1e-18);
}

TEST(CapiEstimation, EmptyViewsAcceptNullOmegaPointers)
{
    po_views* raw = nullptr;
    ASSERT_EQ(po_views_create(0, 2, nullptr, nullptr, &raw), PO_OK);
    ViewsPtr views(raw);
    DataPtr prior = po::testing::make_data(2, kMu, kCovariance);
    EXPECT_EQ(po_bl_omega_he_litterman(prior.get(), views.get(), kTau, nullptr), PO_OK);
    EXPECT_EQ(po_bl_omega_idzorek(prior.get(), views.get(), kTau, nullptr, nullptr), PO_OK);
    EXPECT_EQ(po_views_set_omega(views.get(), nullptr, 0), PO_OK);
}

TEST(CapiEstimation, NullHandles)
{
    DataPtr prior = po::testing::make_data(2, kMu, kCovariance);
    ViewsPtr views = make_views_with_omega();
    const double w_mkt[2] = {0.6, 0.4};
    const double confidence[1] = {0.5};
    double out[4] = {};
    double scalar[2] = {};
    EXPECT_EQ(po_bayes_stein(nullptr, 12, out, &scalar[0], &scalar[1]), PO_INVALID_HANDLE);
    EXPECT_STREQ(po_last_error(), "d must not be null");
    EXPECT_EQ(po_bl_equilibrium(nullptr, w_mkt, 2, 2.5, out), PO_INVALID_HANDLE);
    EXPECT_STREQ(po_last_error(), "d must not be null");
    EXPECT_EQ(po_bl_omega_he_litterman(nullptr, views.get(), kTau, out), PO_INVALID_HANDLE);
    EXPECT_STREQ(po_last_error(), "d must not be null");
    EXPECT_EQ(po_bl_omega_idzorek(nullptr, views.get(), kTau, confidence, out), PO_INVALID_HANDLE);
    EXPECT_STREQ(po_last_error(), "d must not be null");
    EXPECT_EQ(po_bl_posterior(nullptr, views.get(), kTau, out, out), PO_INVALID_HANDLE);
    EXPECT_STREQ(po_last_error(), "d must not be null");
    EXPECT_EQ(po_bl_omega_he_litterman(prior.get(), nullptr, kTau, out), PO_INVALID_HANDLE);
    EXPECT_STREQ(po_last_error(), "v must not be null");
    EXPECT_EQ(po_bl_omega_idzorek(prior.get(), nullptr, kTau, confidence, out), PO_INVALID_HANDLE);
    EXPECT_STREQ(po_last_error(), "v must not be null");
    EXPECT_EQ(po_bl_posterior(prior.get(), nullptr, kTau, out, out), PO_INVALID_HANDLE);
    EXPECT_STREQ(po_last_error(), "v must not be null");
    EXPECT_EQ(po_views_set_omega(nullptr, out, 1), PO_INVALID_HANDLE);
    EXPECT_STREQ(po_last_error(), "v must not be null");
}

TEST(CapiEstimation, NullArguments)
{
    DataPtr prior = po::testing::make_data(2, kMu, kCovariance);
    ViewsPtr views = make_views_with_omega();
    const double w_mkt[2] = {0.6, 0.4};
    const double confidence[1] = {0.5};
    double out[4] = {};
    double grand_mean = 0.0;
    double weight = 0.0;
    po_views* created = nullptr;
    EXPECT_EQ(po_sample_moments(4, 2, nullptr, out, out), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_sample_moments(4, 2, kReturns.data(), nullptr, out), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_sample_moments(4, 2, kReturns.data(), out, nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_bayes_stein(prior.get(), 12, nullptr, &grand_mean, &weight), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_bayes_stein(prior.get(), 12, out, nullptr, &weight), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_bayes_stein(prior.get(), 12, out, &grand_mean, nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_bl_equilibrium(prior.get(), nullptr, 2, 2.5, out), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_bl_equilibrium(prior.get(), w_mkt, 2, 2.5, nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_views_create(1, 2, nullptr, kQ.data(), &created), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_views_create(1, 2, kP.data(), nullptr, &created), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_views_create(1, 2, kP.data(), kQ.data(), nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_views_set_omega(views.get(), nullptr, 1), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_bl_omega_he_litterman(prior.get(), views.get(), kTau, nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_bl_omega_idzorek(prior.get(), views.get(), kTau, nullptr, out), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_bl_omega_idzorek(prior.get(), views.get(), kTau, confidence, nullptr), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_bl_posterior(prior.get(), views.get(), kTau, nullptr, out), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "out_mu must not be null");
    EXPECT_EQ(created, nullptr);
}

TEST(CapiEstimation, OutputUntouchedOnFailure)
{
    DataPtr prior = po::testing::make_data(2, kMu, kCovariance);
    ViewsPtr views = make_views_with_omega();
    const double w_mkt[2] = {0.6, 0.4};
    const double confidence[1] = {0.5};
    double a[4] = {kPoison, kPoison, kPoison, kPoison};
    double b = kPoison;
    double c = kPoison;
    EXPECT_EQ(po_sample_moments(1, 2, kReturns.data(), a, a + 2), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_bayes_stein(prior.get(), 4, a, &b, &c), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_bl_equilibrium(prior.get(), w_mkt, 2, 0.0, a), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_bl_omega_he_litterman(prior.get(), views.get(), 0.0, a), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_bl_omega_idzorek(prior.get(), views.get(), 0.0, confidence, a), PO_INVALID_ARGUMENT);
    EXPECT_EQ(po_bl_posterior(prior.get(), views.get(), 0.0, a, a + 2), PO_INVALID_ARGUMENT);
    for (double value : a)
    {
        EXPECT_EQ(value, kPoison);
    }
    EXPECT_EQ(b, kPoison);
    EXPECT_EQ(c, kPoison);
}

TEST(CapiEstimation, RiskFreeDataIsRefused)
{
    DataPtr prior = po::testing::make_data(2, kMu, kCovariance);
    ASSERT_EQ(po_data_set_risk_free(prior.get(), 0.02), PO_OK);
    const double w_mkt[2] = {0.6, 0.4};
    double out[2] = {kPoison, kPoison};
    EXPECT_EQ(po_bl_equilibrium(prior.get(), w_mkt, 2, 2.5, out), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "data must not have a risk-free rate");
    EXPECT_EQ(out[0], kPoison);
}

TEST(CapiEstimation, OmegaFailureKeepsThePreviousOne)
{
    ViewsPtr views = make_views(1, 2, kP, kQ);
    const double negative[1] = {-1.0};
    EXPECT_EQ(po_views_set_omega(views.get(), negative, 1), PO_INVALID_ARGUMENT);
    DataPtr prior = po::testing::make_data(2, kMu, kCovariance);
    double mu[2] = {};
    EXPECT_EQ(po_bl_posterior(prior.get(), views.get(), kTau, mu, nullptr), PO_INVALID_ARGUMENT);
    const double wrong_k[1] = {0.001};
    EXPECT_EQ(po_views_set_omega(views.get(), wrong_k, 2), PO_INVALID_ARGUMENT);
}

TEST(CapiEstimation, FailedCreateLeavesNull)
{
    const std::vector<double> redundant = {1.0, -1.0, 2.0, -2.0};
    const std::vector<double> q = {0.01, 0.02};
    char sentinel = 0;
    po_views* raw = reinterpret_cast<po_views*>(&sentinel);
    EXPECT_EQ(po_views_create(2, 2, redundant.data(), q.data(), &raw), PO_INVALID_ARGUMENT);
    EXPECT_EQ(raw, nullptr);
    EXPECT_NE(std::string(po_last_error()).find("redundant"), std::string::npos);
    po_views_free(raw);
}

const double* const kUnread = reinterpret_cast<const double*>(std::uintptr_t{8});

TEST(CapiEstimation, ViewsSizeOverflowIsRefusedBeforeReading)
{
    po_views* raw = nullptr;
    EXPECT_EQ(po_views_create(kHuge, kHuge, kUnread, kUnread, &raw), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "views size is too large");
    EXPECT_EQ(raw, nullptr);
}

TEST(CapiEstimation, OmegaSizeOverflowIsRefusedBeforeReading)
{
    ViewsPtr views = make_views_with_omega();
    EXPECT_EQ(po_views_set_omega(views.get(), kUnread, kHuge), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "omega size is too large");
}

TEST(CapiEstimation, ReturnsSizeOverflowIsRefusedBeforeReading)
{
    double out[2] = {};
    EXPECT_EQ(po_sample_moments(kHuge, 2, kUnread, out, out), PO_INVALID_ARGUMENT);
    EXPECT_STREQ(po_last_error(), "returns size is too large");
}

TEST(CapiEstimation, PivotMessage)
{
    const std::vector<double> covariance = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0};
    po_data* raw = nullptr;
    EXPECT_EQ(po_data_create(3, nullptr, covariance.data(), &raw), PO_INVALID_ARGUMENT);
    EXPECT_NE(std::string(po_last_error()).find("pivot 2"), std::string::npos);
}

} // namespace
