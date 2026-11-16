#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <vector>

#include "fuzz/random_problems.hpp"
#include "fuzz/ratio_reference.hpp"
#include "kkt.hpp"
#include "linalg/cholesky.hpp"
#include "solve/solve.hpp"
#include "random_covariance.hpp"
#include "user_kkt.hpp"

namespace
{

using po::Data;
using po::Objective;
using po::Result;
using po::SolveStatus;
using po::Spec;
using po::testing::largest_magnitude;
using po::testing::ratio_of;
using po::testing::ratio_vector;

constexpr double kRiskFree = 0.02;
constexpr double kRiskFreeAboveMeans = 0.2;
// The closed forms, relative to the largest weight.
constexpr double kClosedForm = 1e-12;

Data random_data(std::size_t n, std::uint64_t seed)
{
    std::vector<double> mu = po::testing::random_vector(n, seed + 500, 0.03, 0.15);
    return Data::create(n, mu, po::testing::random_covariance(n, seed)).value();
}

std::vector<double> covariance_solve(const Data& data, std::vector<double> rhs)
{
    std::size_t n = data.asset_count();
    std::vector<double> factor = po::linalg::cholesky(data.covariance(), n).value();
    po::linalg::solve_lower(factor, n, rhs);
    po::linalg::solve_lower_transposed(factor, n, rhs);
    return rhs;
}

double sum(std::span<const double> values)
{
    double total = 0.0;
    for (double value : values)
    {
        total += value;
    }
    return total;
}

Spec ratio_spec(Objective objective, double lower, double upper)
{
    Spec spec;
    spec.objective = objective;
    spec.lower = lower;
    spec.upper = upper;
    return spec;
}

// No bounds, the default.
Spec unbounded_spec(Objective objective)
{
    Spec spec;
    spec.objective = objective;
    return spec;
}

// OPTIMAL with KKT in the user's variables for f = -v'w / vol(w), the gradient written here.
Result solve_checked(const Spec& spec, const Data& data)
{
    Result result = po::solve(spec, data).value();
    std::vector<double> v = ratio_vector(data, spec.objective);
    po::testing::expect_user_kkt(po::testing::with_cash_fixed(spec, data), data, result,
                                 po::testing::ratio_gradient(data, v, result));
    return result;
}

// The closed form without bounds: w = Sigma^-1 v / (1' Sigma^-1 v) over the risky weights.
void expect_closed_form(const Data& data, const Result& result, std::span<const double> v)
{
    std::vector<double> direction = covariance_solve(data, std::vector<double>(v.begin(), v.end()));
    double total = sum(direction);
    ASSERT_GT(total, 0.0);
    for (double& entry : direction)
    {
        entry /= total;
    }
    double scale = largest_magnitude(direction);
    for (std::size_t i = 0; i < direction.size(); ++i)
    {
        EXPECT_NEAR(result.weights[i], direction[i], kClosedForm * scale);
    }
}

void expect_no_point(const Result& result, SolveStatus status)
{
    EXPECT_EQ(result.status, status);
    for (double weight : result.weights)
    {
        EXPECT_TRUE(std::isnan(weight));
    }
    for (double multiplier : result.bound_multipliers)
    {
        EXPECT_TRUE(std::isnan(multiplier));
    }
    EXPECT_TRUE(std::isnan(result.budget_multiplier));
    EXPECT_TRUE(std::isnan(result.leverage_multiplier));
    EXPECT_TRUE(std::isnan(result.residual));
}

void expect_cash_at_zero(const Data& data, const Result& result)
{
    ASSERT_EQ(result.weights.size(), data.asset_count() + 1);
    EXPECT_EQ(result.weights.back(), 0.0);
    EXPECT_TRUE(std::isnan(result.budget_multiplier));
    EXPECT_TRUE(std::isfinite(result.bound_multipliers.back()));
}

TEST(MaxSharpe, UnconstrainedIsTangency)
{
    for (std::uint64_t seed = 1; seed <= 5; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = random_data(6, seed);
        Result result = solve_checked(unbounded_spec(Objective::MaxSharpe), data);
        expect_closed_form(data, result, data.mu());
    }
}

// With rf the tangency of mu - rf 1, the cash weight exactly 0 and the budget term in beta_0: 0 at the
// tangency, where w'grad f = 0 (f is homogeneous of degree 0) leaves no budget term, and not 0 under an
// active cap.
TEST(MaxSharpe, CashWeightZero)
{
    Data data = random_data(6, 7);
    data.set_risk_free(kRiskFree);
    Result tangency = solve_checked(unbounded_spec(Objective::MaxSharpe), data);
    expect_cash_at_zero(data, tangency);
    expect_closed_form(data, tangency, po::testing::excess_mu(data));
    double multiplier_scale = largest_magnitude(tangency.bound_multipliers);
    EXPECT_NEAR(tangency.bound_multipliers.back(), 0.0, 1e-12 * multiplier_scale);
    Result capped = solve_checked(ratio_spec(Objective::MaxSharpe, 0.0, 0.3), data);
    expect_cash_at_zero(data, capped);
    EXPECT_NE(capped.bound_multipliers.back(), 0.0);
    Result excluded = po::solve(ratio_spec(Objective::MaxSharpe, 0.05, 0.3), data).value();
    expect_no_point(excluded, SolveStatus::Infeasible);
}

TEST(MaxSharpe, RiskFreeAboveAllMeansLongOnly)
{
    Data data = random_data(6, 8);
    Spec long_only = ratio_spec(Objective::MaxSharpe, 0.0, 1.0);
    std::vector<double> mu(data.mu().begin(), data.mu().end());
    Data negative = Data::create(6, po::testing::negated(mu), data.covariance()).value();
    expect_no_point(po::solve(long_only, negative).value(), SolveStatus::Infeasible);
    data.set_risk_free(kRiskFreeAboveMeans);
    expect_no_point(po::solve(long_only, data).value(), SolveStatus::Infeasible);
    long_only.leverage = 2.0;
    expect_no_point(po::solve(long_only, data).value(), SolveStatus::Infeasible);
}

// Without bounds and 1' Sigma^-1 mu~ < 0, the tangency direction has a negative sum: along the efficient
// branch the ratio increases toward its supremum, the slope of the frontier's asymptote, reached only as
// 1'y -> 0.
TEST(MaxSharpe, UnboundedWithShorts)
{
    Data data = random_data(6, 9);
    data.set_risk_free(kRiskFreeAboveMeans);
    ASSERT_LT(sum(covariance_solve(data, po::testing::excess_mu(data))), 0.0);
    Result result = po::solve(unbounded_spec(Objective::MaxSharpe), data).value();
    ASSERT_EQ(result.weights.size(), 7u);
    expect_no_point(result, SolveStatus::Unbounded);
    Data negative = Data::create(6, po::testing::excess_mu(data), data.covariance()).value();
    Result without_cash = po::solve(unbounded_spec(Objective::MaxSharpe), negative).value();
    expect_no_point(without_cash, SolveStatus::Unbounded);
}

// Long-only, a leverage with shorts, and both with cash: every case OPTIMAL with KKT, the bounds and the
// leverage binding in some.
TEST(MaxSharpe, LongOnlyAndLeverageKkt)
{
    std::size_t active_bounds = 0;
    std::size_t binding_leverage = 0;
    for (std::uint64_t seed = 1; seed <= 10; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = random_data(8, seed);
        if (seed % 2 == 0)
        {
            data.set_risk_free(kRiskFree);
        }
        Result long_only = solve_checked(ratio_spec(Objective::MaxSharpe, 0.0, 1.0), data);
        Spec leveraged = ratio_spec(Objective::MaxSharpe, -0.5, 1.0);
        leveraged.leverage = 1.4;
        Result with_leverage = solve_checked(leveraged, data);
        for (std::size_t i = 0; i < data.asset_count(); ++i)
        {
            active_bounds += long_only.bound_multipliers[i] != 0.0 ? 1 : 0;
        }
        EXPECT_LE(with_leverage.leverage_multiplier, 0.0);
        binding_leverage += with_leverage.leverage_multiplier < 0.0 ? 1 : 0;
    }
    EXPECT_GE(active_bounds, 10u);
    EXPECT_GE(binding_leverage, 5u);
}

// A random point of the long-only box with upper 0.4: uniform weights normalized, drawn again until they fit.
std::vector<double> random_feasible(std::size_t n, std::mt19937_64& engine)
{
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::vector<double> weights(n);
    for (int draw = 0; draw < 1000; ++draw)
    {
        for (double& weight : weights)
        {
            weight = uniform(engine);
        }
        double total = sum(weights);
        bool fits = true;
        for (double& weight : weights)
        {
            weight /= total;
            fits = fits && weight <= 0.4;
        }
        if (fits)
        {
            return weights;
        }
    }
    ADD_FAILURE() << "no feasible draw";
    return weights;
}

// 10 000 feasible portfolios, the convex combinations (1 - t) w* + t u of the optimum with random feasible u
// at t log-uniform in [1e-6, 1], so that many lie close to the optimum: none has a larger Sharpe ratio.
TEST(MaxSharpe, NoSharperFeasiblePortfolio)
{
    constexpr std::size_t kAssets = 6;
    Data data = random_data(kAssets, 21);
    Result result = solve_checked(ratio_spec(Objective::MaxSharpe, 0.0, 0.4), data);
    double best = ratio_of(data, data.mu(), result.weights);
    std::mt19937_64 engine(2026);
    std::uniform_real_distribution<double> exponent(-6.0, 0.0);
    for (int sample = 0; sample < 10000; ++sample)
    {
        std::vector<double> point = random_feasible(kAssets, engine);
        double t = std::pow(10.0, exponent(engine));
        for (std::size_t i = 0; i < kAssets; ++i)
        {
            point[i] = (1.0 - t) * result.weights[i] + t * point[i];
        }
        double ratio = ratio_of(data, data.mu(), point);
        ASSERT_LE(ratio, best * (1.0 + 1e-12)) << "sample " << sample << " t " << t;
    }
}

// The boundary 1' Sigma^-1 mu~ = 0, exact in the data: mu = Sigma y0 + rf 1 with 1'y0 = 0, computed in integers, so
// that Sigma^-1 mu~ = y0 and rf is exactly the return of the global minimum-variance portfolio. The covariance is
// random_covariance(n, seed, condition) scaled to a largest eigenvalue of about 0.04 (times 2^-14 at the daily
// scale) and rounded to integers times 2^exponent, 2^-36 of its largest entry: the rounding is far below the
// smallest eigenvalue, so the condition number stays that of random_covariance.
struct Boundary
{
    Data data;
    std::vector<double> direction;
};

constexpr int kGridBits = 36;
constexpr std::int64_t kExact = std::int64_t{1} << 53;

// Integers in [-8, 8] summing to 0 (the last one too), times 2^shift.
std::vector<std::int64_t> zero_sum_direction(std::size_t n, std::uint64_t seed, int shift)
{
    std::mt19937_64 engine(seed);
    std::uniform_int_distribution<std::int64_t> draw(-8, 8);
    std::vector<std::int64_t> direction(n);
    std::int64_t total = 100;
    while (total < -8 || total > 8)
    {
        total = 0;
        for (std::size_t i = 0; i + 1 < n; ++i)
        {
            direction[i] = draw(engine);
            total += direction[i];
        }
    }
    direction[n - 1] = -total;
    for (std::int64_t& entry : direction)
    {
        entry *= std::int64_t{1} << shift;
    }
    return direction;
}

// Sigma as integers times 2^exponent.
struct DyadicCovariance
{
    std::vector<std::int64_t> entries;
    int exponent;
};

DyadicCovariance dyadic_covariance(std::size_t n, std::uint64_t seed, double condition, int scale_exponent)
{
    std::vector<double> covariance = po::testing::random_covariance(n, seed, condition);
    int grid = std::ilogb(largest_magnitude(covariance)) - kGridBits;
    DyadicCovariance dyadic{std::vector<std::int64_t>(n * n), grid - std::ilogb(condition) + scale_exponent};
    for (std::size_t k = 0; k < n * n; ++k)
    {
        dyadic.entries[k] = std::llround(std::ldexp(covariance[k], -grid));
    }
    return dyadic;
}

// rf is 2^-6 at the annual scale and 2^-13 at the daily one (2^-14 on Sigma, its square root on the means).
Boundary boundary_problem(std::size_t n, std::uint64_t seed, double condition, bool daily)
{
    int scale_exponent = daily ? -14 : 0;
    DyadicCovariance dyadic = dyadic_covariance(n, seed, condition, scale_exponent);
    std::vector<std::int64_t> direction = zero_sum_direction(n, seed + 100, daily ? 7 : 0);
    double rate = std::ldexp(1.0, daily ? -13 : -6);
    std::int64_t rate_units = std::llround(std::ldexp(rate, -dyadic.exponent));
    std::vector<double> covariance(n * n);
    std::vector<double> mu(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        std::int64_t mean = rate_units;
        for (std::size_t j = 0; j < n; ++j)
        {
            covariance[i * n + j] = std::ldexp(static_cast<double>(dyadic.entries[i * n + j]), dyadic.exponent);
            mean += dyadic.entries[i * n + j] * direction[j];
        }
        EXPECT_LT(std::llabs(mean), kExact);
        mu[i] = std::ldexp(static_cast<double>(mean), dyadic.exponent);
    }
    Boundary boundary{Data::create(n, mu, covariance).value(),
                      std::vector<double>(direction.begin(), direction.end())};
    boundary.data.set_risk_free(rate);
    return boundary;
}

// The lower factor of Sigma in long double, for a reference that does not share the rounding of the library.
std::vector<long double> extended_factor(const Data& data)
{
    std::size_t n = data.asset_count();
    std::vector<long double> factor(data.covariance().begin(), data.covariance().end());
    for (std::size_t j = 0; j < n; ++j)
    {
        for (std::size_t k = 0; k < j; ++k)
        {
            factor[j * n + j] -= factor[j * n + k] * factor[j * n + k];
        }
        factor[j * n + j] = std::sqrt(factor[j * n + j]);
        for (std::size_t i = j + 1; i < n; ++i)
        {
            for (std::size_t k = 0; k < j; ++k)
            {
                factor[i * n + j] -= factor[i * n + k] * factor[j * n + k];
            }
            factor[i * n + j] /= factor[j * n + j];
        }
    }
    return factor;
}

// Sigma^-1 1 in long double.
std::vector<long double> extended_inverse_times_ones(const Data& data)
{
    std::size_t n = data.asset_count();
    std::vector<long double> factor = extended_factor(data);
    std::vector<long double> solution(n, 1.0L);
    for (std::size_t i = 0; i < n; ++i)
    {
        for (std::size_t k = 0; k < i; ++k)
        {
            solution[i] -= factor[i * n + k] * solution[k];
        }
        solution[i] /= factor[i * n + i];
    }
    for (std::size_t i = n; i-- > 0;)
    {
        for (std::size_t k = i + 1; k < n; ++k)
        {
            solution[i] -= factor[k * n + i] * solution[k];
        }
        solution[i] /= factor[i * n + i];
    }
    return solution;
}

// With rf moved down by d > 0, Sigma^-1 mu~ = y0 + d Sigma^-1 1 exactly, so the weights are
// (y0 + d z) / (d 1'z), z = Sigma^-1 1, with no cancellation since 1'y0 = 0 exactly.
void expect_boundary_weights(const Boundary& boundary, const Result& result, double distance, double tolerance)
{
    std::vector<long double> z = extended_inverse_times_ones(boundary.data);
    long double total = 0.0L;
    for (long double entry : z)
    {
        total += entry;
    }
    std::vector<double> expected(z.size());
    for (std::size_t i = 0; i < z.size(); ++i)
    {
        expected[i] = static_cast<double>((boundary.direction[i] + distance * z[i]) / (distance * total));
    }
    double scale = largest_magnitude(expected);
    for (std::size_t i = 0; i < z.size(); ++i)
    {
        EXPECT_NEAR(result.weights[i], expected[i], tolerance * scale);
    }
}

// rf moved by ulps units in the last place, up for ulps > 0.
double moved_rate(double rate, int ulps)
{
    double direction = ulps > 0 ? 1.0 : -1.0;
    for (int step = 0; step < std::abs(ulps); ++step)
    {
        rate = std::nextafter(rate, direction);
    }
    return rate;
}

// rf at the return of the global minimum-variance portfolio, and k ulps away, over cond(Sigma) 1e2 to 1e8 at both
// scales. Without bounds the supremum is attained only when 1' Sigma^-1 mu~ > 0: at rf and above it the result
// must be UNBOUNDED; below it, where the exact weights are of order 1 / (k eps), UNBOUNDED (kappa within its noise)
// or OPTIMAL with the weights of the extended-precision reference within the band.
TEST(MaxSharpe, RiskFreeAtMinimumVarianceReturn)
{
    constexpr std::array kConditions = {1e2, 1e4, 1e6, 1e8};
    constexpr std::array kUlps = {0, 1, -1, 3, -3, 16, -16};
    for (double condition : kConditions)
    {
        for (std::uint64_t seed = 1; seed <= 4; ++seed)
        {
            SCOPED_TRACE(testing::Message() << "cond " << condition << " seed " << seed);
            Boundary boundary = boundary_problem(3 + 3 * seed, seed, condition, seed % 2 == 0);
            double rate = *boundary.data.risk_free();
            for (int ulps : kUlps)
            {
                SCOPED_TRACE(ulps);
                Boundary moved = boundary;
                moved.data.set_risk_free(moved_rate(rate, ulps));
                Result result = po::solve(unbounded_spec(Objective::MaxSharpe), moved.data).value();
                if (ulps >= 0 || result.status != SolveStatus::Optimal)
                {
                    expect_no_point(result, SolveStatus::Unbounded);
                    continue;
                }
                expect_boundary_weights(moved, result, rate - *moved.data.risk_free(),
                                        po::testing::condition_tolerance(condition));
            }
        }
    }
}

// The boundary with some weights bounded: F, the first nf assets, without bounds, and B, the last 3, in the box
// [0, 1]. Two free assets keep the homogeneous cone unbounded, so the noise of kappa decides. mu_F is that of
// boundary_problem over F (direction y0 on F, 0 on B) and mu~_B = (Sigma y0)_B - rf, so that at kappa = 0 the
// lower rows of B are active with multiplier nu rf > 0 and what is left is the boundary
// 1' Sigma_FF^-1 mu~_F = 0 of the problem over F.
struct PartialBoundary
{
    Data data;
    Data free_block;
    std::vector<double> direction;
};

constexpr std::size_t kBoundedAssets = 3;

// The leading nf x nf block of Sigma, with the means of F.
Data free_block(const Data& data, std::size_t nf)
{
    std::size_t n = data.asset_count();
    std::vector<double> block(nf * nf);
    for (std::size_t k = 0; k < nf * nf; ++k)
    {
        block[k] = data.covariance()[(k / nf) * n + k % nf];
    }
    std::vector<double> mu(data.mu().begin(), data.mu().begin() + static_cast<std::ptrdiff_t>(nf));
    return Data::create(nf, mu, block).value();
}

PartialBoundary partial_boundary_problem(std::size_t nf, std::uint64_t seed, double condition, bool daily)
{
    std::size_t n = nf + kBoundedAssets;
    DyadicCovariance dyadic = dyadic_covariance(n, seed, condition, daily ? -14 : 0);
    std::vector<std::int64_t> direction = zero_sum_direction(nf, seed + 100, daily ? 7 : 0);
    direction.resize(n, 0);
    double rate = std::ldexp(1.0, daily ? -13 : -6);
    std::int64_t rate_units = std::llround(std::ldexp(rate, -dyadic.exponent));
    std::vector<double> covariance(n * n);
    std::vector<double> mu(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        std::int64_t mean = i < nf ? rate_units : 0;
        for (std::size_t j = 0; j < n; ++j)
        {
            covariance[i * n + j] = std::ldexp(static_cast<double>(dyadic.entries[i * n + j]), dyadic.exponent);
            mean += dyadic.entries[i * n + j] * direction[j];
        }
        EXPECT_LT(std::llabs(mean), kExact);
        mu[i] = std::ldexp(static_cast<double>(mean), dyadic.exponent);
    }
    Data data = Data::create(n, mu, covariance).value();
    data.set_risk_free(rate);
    return PartialBoundary{data, free_block(data, nf), std::vector<double>(direction.begin(), direction.end())};
}

// The free weights unbounded, those of B in [0, 1], cash unbounded.
Spec partial_spec(std::size_t n)
{
    constexpr double kInfinity = std::numeric_limits<double>::infinity();
    Spec spec = unbounded_spec(Objective::MaxSharpe);
    spec.asset_lower.assign(n + 1, -kInfinity);
    spec.asset_upper.assign(n + 1, kInfinity);
    for (std::size_t i = n - kBoundedAssets; i < n; ++i)
    {
        spec.asset_lower[i] = 0.0;
        spec.asset_upper[i] = 1.0;
    }
    return spec;
}

// The exact optimum at rf moved down by d > 0, in long double: y = (y0 + d z_F) / (y0' Sigma_FF y0 + d^2 1'z_F) on F
// and 0 on B, with z_F = Sigma_FF^-1 1 (v_F = Sigma_FF y0 + d 1 and 1'y0 = 0 leave no cancellation).
std::vector<double> partial_exact_y(const PartialBoundary& boundary, double distance)
{
    std::size_t nf = boundary.free_block.asset_count();
    std::vector<long double> z = extended_inverse_times_ones(boundary.free_block);
    long double curvature = 0.0L;
    long double total = 0.0L;
    for (std::size_t i = 0; i < nf; ++i)
    {
        total += z[i];
        for (std::size_t j = 0; j < nf; ++j)
        {
            curvature += boundary.direction[i] *
                         static_cast<long double>(boundary.free_block.covariance()[i * nf + j]) *
                         boundary.direction[j];
        }
    }
    long double denominator = curvature + static_cast<long double>(distance) * distance * total;
    std::vector<double> y(boundary.data.asset_count(), 0.0);
    for (std::size_t i = 0; i < nf; ++i)
    {
        y[i] = static_cast<double>((boundary.direction[i] + distance * z[i]) / denominator);
    }
    return y;
}

// The noise model of src/problems/homogeneous.cpp evaluated at a given y, written again here:
// n eps (|z|'|Sigma||y| + sigma_y^2 |z|'m), z = Sigma^-1 1 of the whole covariance, m = |mu| + |rf|.
double kappa_noise_at(const Data& data, std::span<const double> y)
{
    std::size_t n = data.asset_count();
    std::vector<double> z = covariance_solve(data, std::vector<double>(n, 1.0));
    double variance = 0.0;
    double curvature = 0.0;
    double rounding = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        double magnitude = 0.0;
        double product = 0.0;
        for (std::size_t j = 0; j < n; ++j)
        {
            magnitude += std::fabs(data.covariance()[i * n + j] * y[j]);
            product += data.covariance()[i * n + j] * y[j];
        }
        variance += y[i] * product;
        curvature += std::fabs(z[i]) * magnitude;
        rounding += std::fabs(z[i]) * (std::fabs(data.mu()[i]) + std::fabs(*data.risk_free()));
    }
    return static_cast<double>(n) * std::numeric_limits<double>::epsilon() * (curvature + variance * rounding);
}

struct PartialCounts
{
    std::size_t optimal = 0;
    std::size_t unbounded = 0;
};

// An OPTIMAL point: the budget within the band relative to the largest weight (the weights near the boundary are
// far from order 1), the bounds within the band, and the lower rows of B active.
void expect_partial_point(const Spec& spec, const Data& data, const Result& result, double tolerance)
{
    double total = 0.0;
    for (double weight : result.weights)
    {
        total += weight;
    }
    EXPECT_NEAR(total, 1.0, tolerance * largest_magnitude(result.weights));
    po::testing::expect_user_bounds(spec, result, tolerance);
    for (std::size_t i = data.asset_count() - kBoundedAssets; i < data.asset_count(); ++i)
    {
        EXPECT_GT(result.bound_multipliers[i], 0.0);
    }
}

// Twice the threshold of 16 kappa_noise. The measured error of the computed kappa is at most 0.040 kappa_noise
// (docs/guide/validation.md), so an exact kappa above this cannot fall below the threshold: such a problem is
// decided by the data and must be OPTIMAL.
constexpr double kDecidedKappa = 32.0;

void check_partial_rate(const PartialBoundary& boundary, double moved, double tolerance, PartialCounts& counts)
{
    PartialBoundary problem = boundary;
    problem.data.set_risk_free(moved);
    Spec spec = partial_spec(problem.data.asset_count());
    Result result = po::solve(spec, problem.data).value();
    counts.optimal += result.status == SolveStatus::Optimal ? 1 : 0;
    counts.unbounded += result.status == SolveStatus::Unbounded ? 1 : 0;
    double distance = *boundary.data.risk_free() - moved;
    if (distance <= 0.0)
    {
        expect_no_point(result, SolveStatus::Unbounded);
        return;
    }
    std::vector<double> y = partial_exact_y(boundary, distance);
    double kappa = sum(y);
    if (kappa > kDecidedKappa * kappa_noise_at(problem.data, y))
    {
        ASSERT_EQ(result.status, SolveStatus::Optimal) << "kappa " << kappa;
    }
    if (result.status == SolveStatus::Optimal)
    {
        expect_partial_point(po::testing::with_cash_fixed(spec, problem.data), problem.data, result, tolerance);
        return;
    }
    expect_no_point(result, SolveStatus::Unbounded);
}

// rf at the return of the global minimum-variance portfolio of F, k ulps away (k in {0, +-1, +-3, +-16}), and
// moved down by 3 to 3e4 times cond(Sigma) eps relative, over cond(Sigma) 1e2 to 1e8 at both scales. No exact
// kappa <= 0 is OPTIMAL, and no exact kappa above twice the threshold is UNBOUNDED. Of the 192 solves, 70 are
// OPTIMAL (exact kappa from 19 kappa_noise) and 122 UNBOUNDED, 58 of them with an exact kappa > 0 of at most
// 8.0 kappa_noise (docs/guide/validation.md).
TEST(MaxSharpe, PartiallyBoundedAtMinimumVarianceReturn)
{
    constexpr std::array kConditions = {1e2, 1e4, 1e6, 1e8};
    constexpr std::array kUlps = {0, 1, -1, 3, -3, 16, -16};
    constexpr std::array kMoves = {3.0, 30.0, 300.0, 3e3, 3e4};
    PartialCounts counts;
    for (double condition : kConditions)
    {
        for (std::uint64_t seed = 1; seed <= 4; ++seed)
        {
            SCOPED_TRACE(testing::Message() << "cond " << condition << " seed " << seed);
            PartialBoundary boundary = partial_boundary_problem(3 + 3 * seed, seed, condition, seed % 2 == 0);
            double rate = *boundary.data.risk_free();
            std::vector<double> rates;
            for (int ulps : kUlps)
            {
                rates.push_back(moved_rate(rate, ulps));
            }
            for (double move : kMoves)
            {
                rates.push_back(rate * (1.0 - move * condition * std::numeric_limits<double>::epsilon()));
            }
            for (double moved : rates)
            {
                SCOPED_TRACE(moved);
                check_partial_rate(boundary, moved, po::testing::condition_tolerance(condition), counts);
            }
        }
    }
    EXPECT_GE(counts.optimal, 60u);
    EXPECT_GE(counts.unbounded, 100u);
}

// The diversification ratio at 1' Sigma^-1 sigma = 0, exact in the data. With Sigma = D C D, D = diag(sigma),
// Sigma^-1 sigma = D^-1 C^-1 1. The correlation C = [[1, q e, e], [q e, 1, e], [e, e, 1]], e = (1 + q) / 2, has
// C^-1 1 = (1, 1, -q), so 1' Sigma^-1 sigma = 2 / s - q / sigma_3, which is 0 at sigma = s (1, 1, q / 2). With
// q = 1 - 2^-k every entry is a dyadic rational of at most 2k + 1 bits, exact in double, and sqrt(Sigma_33) is
// sigma_3 exactly. C is positive definite for q < 1 (eigenvalues 1 - q e and those of
// [[1 + q e, sqrt(2) e], [sqrt(2) e, 1]], of determinant (1 - q) / 2); k = 1, 8, 15 and 21 give cond(Sigma) of
// 1.3e2, 1.0e4, 1.3e6 and 8.5e7. Copies of the block on the diagonal keep the sum at 0.
Data zero_diversification_problem(int k, double s, std::size_t blocks)
{
    double q = 1.0 - std::ldexp(1.0, -k);
    double e = 0.5 * (1.0 + q);
    std::array<double, 3> sigma = {s, s, 0.5 * q * s};
    std::array<std::array<double, 3>, 3> correlation = {{{1.0, q * e, e}, {q * e, 1.0, e}, {e, e, 1.0}}};
    std::size_t n = 3 * blocks;
    std::vector<double> covariance(n * n, 0.0);
    for (std::size_t entry = 0; entry < 9 * blocks; ++entry)
    {
        std::size_t block = entry / 9;
        std::size_t i = (entry % 9) / 3;
        std::size_t j = entry % 3;
        covariance[(3 * block + i) * n + 3 * block + j] = sigma[i] * sigma[j] * correlation[i][j];
    }
    return Data::create(n, {}, covariance).value();
}

// Without bounds the supremum is not attained at 1' Sigma^-1 sigma = 0: UNBOUNDED, over the four condition numbers,
// both scales and one to three blocks.
TEST(MaxDiversification, InverseCovarianceTimesVolatilitiesSumsToZero)
{
    constexpr std::array kBits = {1, 8, 15, 21};
    constexpr std::array kScales = {0.25, 0.25 / 128.0};
    for (int k : kBits)
    {
        for (double s : kScales)
        {
            for (std::size_t blocks = 1; blocks <= 3; ++blocks)
            {
                SCOPED_TRACE(testing::Message() << "k " << k << " s " << s << " blocks " << blocks);
                Data data = zero_diversification_problem(k, s, blocks);
                Result result = po::solve(unbounded_spec(Objective::MaxDiversification), data).value();
                expect_no_point(result, SolveStatus::Unbounded);
            }
        }
    }
}

TEST(MaxDiversification, UnconstrainedClosedForm)
{
    for (std::uint64_t seed = 1; seed <= 5; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = random_data(6, seed);
        Result result = solve_checked(unbounded_spec(Objective::MaxDiversification), data);
        expect_closed_form(data, result, ratio_vector(data, Objective::MaxDiversification));
    }
}

TEST(MaxDiversification, LongOnlyKkt)
{
    std::size_t active_bounds = 0;
    for (std::uint64_t seed = 1; seed <= 10; ++seed)
    {
        SCOPED_TRACE(seed);
        Data data = Data::create(8, {}, po::testing::random_covariance(8, seed, 1e4)).value();
        Result result = solve_checked(ratio_spec(Objective::MaxDiversification, 0.0, 0.3), data);
        for (double beta : result.bound_multipliers)
        {
            active_bounds += beta != 0.0 ? 1 : 0;
        }
    }
    EXPECT_GE(active_bounds, 10u);
}

// The ratio does not depend on cash: with cash fixed at 0 the risky weights are those without cash, bit for
// bit, the budget multiplier becomes beta_0 = -lambda_budget and mu is not needed.
TEST(MaxDiversification, CashWeightZero)
{
    Data data = Data::create(6, {}, po::testing::random_covariance(6, 4, 1e4)).value();
    Spec spec = ratio_spec(Objective::MaxDiversification, 0.0, 0.3);
    Result without_cash = solve_checked(spec, data);
    data.set_risk_free(kRiskFree);
    Result with_cash = solve_checked(spec, data);
    expect_cash_at_zero(data, with_cash);
    for (std::size_t i = 0; i < data.asset_count(); ++i)
    {
        EXPECT_EQ(with_cash.weights[i], without_cash.weights[i]);
        EXPECT_EQ(with_cash.bound_multipliers[i], without_cash.bound_multipliers[i]);
    }
    EXPECT_EQ(with_cash.bound_multipliers.back(), -without_cash.budget_multiplier);
    EXPECT_NE(with_cash.bound_multipliers.back(), 0.0);
}

// The covariance of random_covariance(n, seed, condition) rescaled to eigenvalues in [scale / condition, scale],
// and means N(0.05 scale, 0.3 sqrt(scale)) of either sign.
Data large_condition_data(std::size_t n, std::uint64_t seed, double condition, double scale)
{
    std::vector<double> covariance = po::testing::random_covariance(n, seed, condition);
    for (double& entry : covariance)
    {
        entry *= scale / (0.04 * condition);
    }
    std::mt19937_64 engine(seed + 900);
    std::vector<double> mu(n);
    for (double& mean : mu)
    {
        mean = 0.3 * std::sqrt(scale) * po::testing::standard_normal(engine) + 0.05 * scale;
    }
    return Data::create(n, mu, covariance).value();
}

Spec leveraged_spec(Objective objective, double lower, double upper, double leverage)
{
    Spec spec = ratio_spec(objective, lower, upper);
    spec.leverage = leverage;
    return spec;
}

// Above cond(Sigma) = 1e8 the stationarity is not checked: where many weights are interior its residual is about
// 0.1 cond(Sigma) eps relative, the limit of a backward-stable solve in double (docs/guide/validation.md). The
// status, the budget, the bounds and the leverage are.
void expect_feasible_optimum(const Spec& spec, const Data& data, const Result& result, double tolerance)
{
    ASSERT_EQ(result.status, SolveStatus::Optimal);
    Spec problem = po::testing::with_cash_fixed(spec, data);
    po::testing::expect_user_equalities(problem, data, result, tolerance);
    po::testing::expect_user_bounds(problem, result, tolerance);
    po::testing::expect_user_leverage(problem, data, result, tolerance);
}

// Long-only passes the whole KKT check at cond(Sigma) = 1e12 with a margin of 74 below the band (the worst
// stationarity 1.4e-9 at n = 150, 5.4e-9 at n = 50); at 1e13 the margin is 2.6, so it is not checked there.
constexpr double kLongOnlyKktCondition = 1e12;

void check_bounded_cones(Objective objective, const Data& data, double condition)
{
    constexpr double kInfinity = std::numeric_limits<double>::infinity();
    std::array specs = {ratio_spec(objective, -1.0, 1.0), ratio_spec(objective, 0.0, 1.0),
                        ratio_spec(objective, 0.0, kInfinity), leveraged_spec(objective, -kInfinity, kInfinity, 2.0),
                        leveraged_spec(objective, -0.5, 1.0, 1.5)};
    double tolerance = po::testing::condition_tolerance(condition);
    std::vector<double> v = ratio_vector(data, objective);
    for (std::size_t k = 0; k < specs.size(); ++k)
    {
        SCOPED_TRACE(k);
        Result result = po::solve(specs[k], data).value();
        expect_feasible_optimum(specs[k], data, result, tolerance);
        bool long_only = specs[k].lower == 0.0 && !specs[k].leverage.has_value();
        if (long_only && condition == kLongOnlyKktCondition)
        {
            po::testing::expect_user_kkt(specs[k], data, result, po::testing::ratio_gradient(data, v, result),
                                         tolerance);
        }
    }
}

// Bounded feasible sets, where UNBOUNDED is impossible (a leverage, or no pair i != j with upper_i = +inf and
// lower_j = -inf): a box, long-only with and without a cap, a leverage alone and a leverage with a box, n = 150,
// cond(Sigma) 1e12 and 1e14, both scales. Each is OPTIMAL; before the noise threshold was limited to unbounded
// cones, the box at 1e12 and every case at 1e14 were UNBOUNDED.
void expect_bounded_cones_optimal(Objective objective)
{
    constexpr std::array kConditions = {1e12, 1e14};
    constexpr std::array kScales = {1.0, 1.0 / 252.0};
    for (double condition : kConditions)
    {
        for (std::uint64_t seed = 1; seed <= 3; ++seed)
        {
            for (double scale : kScales)
            {
                SCOPED_TRACE(testing::Message() << "cond " << condition << " seed " << seed << " scale " << scale);
                check_bounded_cones(objective, large_condition_data(150, seed, condition, scale), condition);
            }
        }
    }
}

TEST(MaxSharpe, BoundedConeAtLargeCondition)
{
    expect_bounded_cones_optimal(Objective::MaxSharpe);
}

TEST(MaxDiversification, BoundedConeAtLargeCondition)
{
    expect_bounded_cones_optimal(Objective::MaxDiversification);
}

} // namespace
