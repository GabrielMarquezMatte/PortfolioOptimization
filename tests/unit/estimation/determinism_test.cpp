// The committed bits of the estimation utilities: one fixed case run from a literal, never generated, panel
// and compared bit for bit (std::bit_cast to std::uint64_t) with hex-float values recorded here. They are the
// guard that the sample moments, the Bayes-Stein shrinkage, the equilibrium returns, the Idzorek Omega and the
// Black-Litterman posterior return the same bits on every platform and compiler of the CI matrix.
//
// A legitimate change of bits (a different summation order, a new rounding in src/) rewrites these values in
// the same commit. A platform that disagrees is a defect to be fixed at its cause; the values are never
// loosened and never made conditional. On a mismatch the failure prints the whole recorded form of the actual
// results, in hex-float without the 0x prefix.

#include <gtest/gtest.h>

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <format>
#include <span>
#include <string>
#include <vector>

#include "core/data.hpp"
#include "estimation/bayes_stein.hpp"
#include "estimation/black_litterman.hpp"
#include "estimation/moments.hpp"
#include "estimation/posterior.hpp"
#include "estimation/views.hpp"

namespace
{

using po::Data;
using namespace po::estimation;

constexpr std::size_t kT = 8;
constexpr std::size_t kN = 3;

// Eight periods of three asset returns; the columns are not collinear, so the covariance is positive definite.
constexpr std::array<double, kT * kN> kReturns = {0.012,  -0.008, 0.021,  //
                                                  -0.015, 0.011,  -0.030, //
                                                  0.024,  0.019,  0.040,  //
                                                  0.006,  -0.022, 0.013,  //
                                                  -0.010, 0.016,  -0.018, //
                                                  0.017,  0.003,  0.027,  //
                                                  -0.004, -0.013, 0.008,  //
                                                  0.009,  0.025,  -0.011};
constexpr std::array<double, kN> kMarketWeights = {0.5, 0.3, 0.2};
constexpr double kDelta = 2.5;
constexpr double kTau = 0.05;
constexpr std::array<double, 6> kP = {1, -1, 0, 0, 1, -1};
constexpr std::array<double, 2> kQ = {0.01, 0.005};
constexpr std::array<double, 2> kConfidence = {0.6, 0.9};

struct Golden
{
    std::vector<double> sample_mu;
    std::vector<double> sample_covariance;
    std::vector<double> shrunk_mu;
    double grand_mean;
    double weight;
    std::vector<double> pi;
    std::vector<double> omega;
    std::vector<double> posterior_mu;
    std::vector<double> posterior_covariance;
};

std::uint64_t bits(double value)
{
    return std::bit_cast<std::uint64_t>(value);
}

std::string record_of(const std::vector<double>& values)
{
    std::string text = "{";
    for (double value : values)
    {
        text += std::format("{:a}, ", value);
    }
    return text + "}";
}

std::string record_of(const Golden& g)
{
    return std::format("sample_mu {}\nsample_covariance {}\nshrunk_mu {}\ngrand_mean {:a}, weight {:a}\npi {}\n"
                       "omega {}\nposterior_mu {}\nposterior_covariance {}",
                       record_of(g.sample_mu), record_of(g.sample_covariance), record_of(g.shrunk_mu),
                       g.grand_mean, g.weight, record_of(g.pi), record_of(g.omega), record_of(g.posterior_mu),
                       record_of(g.posterior_covariance));
}

bool same_bits(const std::vector<double>& actual, const std::vector<double>& expected)
{
    if (actual.size() != expected.size())
    {
        return false;
    }
    for (std::size_t i = 0; i < actual.size(); ++i)
    {
        if (bits(actual[i]) != bits(expected[i]))
        {
            return false;
        }
    }
    return true;
}

bool same_bits(const Golden& a, const Golden& b)
{
    return same_bits(a.sample_mu, b.sample_mu) && same_bits(a.sample_covariance, b.sample_covariance) &&
           same_bits(a.shrunk_mu, b.shrunk_mu) && bits(a.grand_mean) == bits(b.grand_mean) &&
           bits(a.weight) == bits(b.weight) && same_bits(a.pi, b.pi) && same_bits(a.omega, b.omega) &&
           same_bits(a.posterior_mu, b.posterior_mu) &&
           same_bits(a.posterior_covariance, b.posterior_covariance);
}

// Runs the whole chain on the fixed panel; every step must succeed.
Golden run_chain()
{
    Golden g;
    Moments moments = sample_moments(kReturns, kT, kN).value();
    g.sample_mu = moments.mu;
    g.sample_covariance = moments.covariance;

    Data sample = Data::create(kN, moments.mu, moments.covariance).value();
    BayesStein shrunk = bayes_stein(sample, kT).value();
    g.shrunk_mu = shrunk.mu;
    g.grand_mean = shrunk.grand_mean;
    g.weight = shrunk.weight;

    g.pi = equilibrium(sample, kMarketWeights, kDelta).value();
    Data prior = Data::create(kN, g.pi, moments.covariance).value();
    Views views = Views::create(2, kN, kP, kQ).value();
    g.omega = omega_idzorek(prior, views, kTau, kConfidence).value();
    EXPECT_TRUE(views.set_omega(g.omega).has_value());

    Posterior post = posterior(prior, views, kTau, true).value();
    g.posterior_mu = post.mu;
    g.posterior_covariance = post.covariance;
    return g;
}

TEST(EstimationDeterminism, TheChainHasTheCommittedBits)
{
    const Golden golden{
        .sample_mu = {0x1.3f7ced916872bp-8, 0x1.fbe76c8b43958p-9, 0x1.999999999999ap-8},
        .sample_covariance = {
            0x1.7e8b01ac6b0ebp-13, 0x1.7ee14a13a591ap-17, 0x1.2c3ad2bfee795p-12, 0x1.7ee14a13a591ap-17,
            0x1.26ee43f72bde9p-12, -0x1.e19891f541bc1p-14, 0x1.2c3ad2bfee795p-12, -0x1.e19891f541bc1p-14,
            0x1.2e40852b4d8bap-11},
        .shrunk_mu = {0x1.1a07bfeb15ecep-8, 0x1.197073f70dd84p-8, 0x1.1ad7c85aa1093p-8},
        .grand_mean = 0x1.19b07cd8afbdep-8,
        .weight = 0x1.fb61fe9d80df1p-1,
        .pi = {0x1.8e2d922830074p-12, 0x1.5fe8db3efa4e3p-13, 0x1.259eb89e9da8p-11},
        .omega = {0x1.ed1722513ab86p-17, 0x0p+0, 0x0p+0, 0x1.956a856f13d97p-18},
        .posterior_mu = {-0x1.a3eedf3a17758p-13, 0x1.6e02590c98163p-13, -0x1.f6b265de970cp-10},
        .posterior_covariance = {
            0x1.8b0393901c085p-13, 0x1.16baf0f3c03d2p-16, 0x1.32f2ddd37e0fdp-12, 0x1.16baf0f3c03d2p-16,
            0x1.2e8258eb9b176p-12, -0x1.c9a6739c52199p-14, 0x1.32f2ddd37e0fdp-12, -0x1.c9a6739c52199p-14,
            0x1.32f9515af60abp-11},
    };
    Golden actual = run_chain();
    EXPECT_TRUE(same_bits(actual, golden)) << "actual result:\n" << record_of(actual);
}

} // namespace
