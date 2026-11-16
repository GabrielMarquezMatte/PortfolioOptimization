#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <span>
#include <vector>

#include "fuzz/random_problems.hpp"
#include "fuzz/ratio_reference.hpp"
#include "kkt.hpp"
#include "linalg/cholesky.hpp"
#include "problems/solve.hpp"
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

} // namespace
