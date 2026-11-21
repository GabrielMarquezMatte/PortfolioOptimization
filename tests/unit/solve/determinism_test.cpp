// The committed solution bits of the library: one case per objective, solved from literal, never generated,
// inputs and compared bit for bit (std::bit_cast to std::uint64_t) with hex-float values recorded here. They
// are the guard that a solve returns the same bits on every platform and compiler of the CI matrix.
//
// A legitimate change of bits (a different pivot order, a new rounding in src/) rewrites these values in the
// same commit, with the output of `quick.py compare-solutions` in the commit body. A platform that disagrees
// is a defect to be fixed at its cause (spec portable-and-max-return, section 6.4); the values are never
// loosened and never made conditional. On a mismatch the failure prints the whole recorded form of the actual
// result, in hex-float without the 0x prefix.

#include <gtest/gtest.h>

#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <format>
#include <limits>
#include <span>
#include <string>
#include <vector>

#include "solve/solve.hpp"

namespace
{

using po::Data;
using po::Objective;
using po::Result;
using po::SolveStatus;
using po::Spec;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// n = 4 assets: volatilities 0.15, 0.20, 0.25, 0.12 and correlations between 0.10 and 0.30.
constexpr std::array<double, 4> kMu = {0.071, 0.094, 0.118, 0.052};
constexpr std::array<double, 16> kCovariance = {0.0225, 0.0066, 0.0090, 0.0027, //
                                                0.0066, 0.0400, 0.0150, 0.0036, //
                                                0.0090, 0.0150, 0.0625, 0.0030, //
                                                0.0027, 0.0036, 0.0030, 0.0144};

// Eight periods of four asset returns and the benchmark return of each period.
constexpr std::array<double, 32> kReturns = {0.012,  -0.008, 0.021,  0.004,  //
                                             -0.015, 0.011,  -0.030, 0.007,  //
                                             0.024,  0.019,  0.040,  -0.002, //
                                             0.006,  -0.022, 0.013,  0.009,  //
                                             -0.010, 0.016,  -0.018, 0.001,  //
                                             0.017,  0.003,  0.027,  0.012,  //
                                             -0.004, -0.013, 0.008,  -0.006, //
                                             0.009,  0.025,  0.011,  0.005};
constexpr std::array<double, 8> kBenchmark = {0.010, -0.009, 0.020, -0.001, 0.003, 0.012, -0.007, 0.014};

constexpr double kRiskFree = 0.02;

struct Golden
{
    SolveStatus status;
    int iterations;
    std::vector<double> weights;
    double budget;
    double target_return;
    double volatility;
    double leverage;
    std::vector<double> bounds;
    double residual;
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

std::string record_of(const Result& result)
{
    return std::format(
        "status {}, iterations {}\nweights {}\nbudget {:a}, target_return {:a}, volatility {:a}, leverage {:a}\n"
        "bounds {}\nresidual {:a}",
        static_cast<int>(result.status), result.iterations, record_of(result.weights),
        result.budget_multiplier, result.target_return_multiplier, result.volatility_multiplier,
        result.leverage_multiplier, record_of(result.bound_multipliers), result.residual);
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

// The five scalars of a result in one vector, compared by the same bit rule as the vectors.
std::vector<double> scalars_of(const Result& result)
{
    return {result.budget_multiplier, result.target_return_multiplier, result.volatility_multiplier,
            result.leverage_multiplier, result.residual};
}

std::vector<double> scalars_of(const Golden& golden)
{
    return {golden.budget, golden.target_return, golden.volatility, golden.leverage, golden.residual};
}

void expect_golden(const Spec& spec, const Data& data, const Golden& golden)
{
    std::expected<Result, po::Failure> solved = po::solve(spec, data);
    ASSERT_TRUE(solved.has_value()) << solved.error().message;
    const Result& result = solved.value();
    bool same = result.status == golden.status && result.iterations == golden.iterations &&
                same_bits(result.weights, golden.weights) &&
                same_bits(result.bound_multipliers, golden.bounds) &&
                same_bits(scalars_of(result), scalars_of(golden));
    EXPECT_TRUE(same) << "actual result:\n" << record_of(result);
}

Data covariance_data(bool with_mu, bool with_cash = false)
{
    std::span<const double> mu;
    if (with_mu)
    {
        mu = kMu;
    }
    Data data = Data::create(4, mu, kCovariance).value();
    if (with_cash)
    {
        data.set_risk_free(kRiskFree);
    }
    return data;
}

TEST(Determinism, MinVarianceHasTheCommittedBits)
{
    // Target return and an active leverage, with shorts allowed.
    Spec spec;
    spec.lower = -0.4;
    spec.upper = 0.9;
    spec.target_return = 0.11;
    spec.leverage = 1.2;
    const Golden golden{
        .status = SolveStatus::Optimal,
        .iterations = 3,
        .weights = {0x1.747db7cf9595dp-4, 0x1.b8916246a3088p-2, 0x1.285acb15eefc3p-1, -0x1.999999999999p-4},
        .budget = -0x1.34649a37b382cp-4,
        .target_return = 0x1.6b024addb1e51p+0,
        .volatility = kNaN,
        .leverage = -0x1.77f74ffe7e7f1p-8,
        .bounds = {0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0},
        .residual = 0x0p+0,
    };
    expect_golden(spec, covariance_data(true), golden);
}

TEST(Determinism, MaxReturnHasTheCommittedBits)
{
    // A volatility cap that binds, with cash.
    Spec spec;
    spec.objective = Objective::MaxReturn;
    spec.lower = 0.0;
    spec.upper = 0.24;
    spec.max_volatility = 0.09;
    const Golden golden{
        .status = SolveStatus::Optimal,
        .iterations = 32,
        .weights = {0x1.c7e6e4c0eb887p-3, 0x1.751fe216175f2p-3, 0x1.5ded21c94b4f4p-3,
                    0x1.eb851eb851eb8p-3, 0x1.7986f8a75fddcp-3},
        .budget = kNaN,
        .target_return = kNaN,
        .volatility = -0x1.17dfe8348e99bp-1,
        .leverage = kNaN,
        .bounds = {0x0p+0, 0x0p+0, 0x0p+0, -0x1.12b1688d42f57p-12, 0x0p+0},
        .residual = 0x0p+0,
    };
    expect_golden(spec, covariance_data(true, true), golden);
}

TEST(Determinism, UtilityHasTheCommittedBits)
{
    Spec spec;
    spec.objective = Objective::Utility;
    spec.risk_aversion = 4.0;
    spec.lower = 0.0;
    spec.weight_upper = {0.25, 0.35, 0.3, 0.4, 1.0};
    const Golden golden{
        .status = SolveStatus::Optimal,
        .iterations = 2,
        .weights = {0x1p-2, 0x1.02cc304108c63p-2, 0x1.ffa51d4393cffp-3, 0x1.fac2823a5aa38p-3, 0x1.8p-54},
        .budget = kNaN,
        .target_return = kNaN,
        .volatility = kNaN,
        .leverage = kNaN,
        .bounds = {-0x1.cb2dd787ec023p-10, 0x0p+0, 0x0p+0, 0x0p+0, 0x1.139531f6eec23p-7},
        .residual = 0x0p+0,
    };
    expect_golden(spec, covariance_data(true, true), golden);
}

TEST(Determinism, ImportanceHasTheCommittedBits)
{
    Spec spec;
    spec.objective = Objective::Importance;
    spec.importance = 0.4;
    spec.lower = 0.0;
    spec.upper = 0.3;
    const Golden golden{
        .status = SolveStatus::Optimal,
        .iterations = 2,
        .weights = {0x1.0f3114e9cad35p-2, 0x1.312e38139ba01p-2, 0x1.3333333333333p-2, 0x1.18daff9eccb2ep-3},
        .budget = -0x1.f23a8dbc95e61p-7,
        .target_return = kNaN,
        .volatility = kNaN,
        .leverage = kNaN,
        .bounds = {0x0p+0, 0x0p+0, -0x1.974abd529e373p-11, 0x0p+0},
        .residual = 0x0p+0,
    };
    expect_golden(spec, covariance_data(true), golden);
}

TEST(Determinism, KellyHasTheCommittedBits)
{
    Spec spec;
    spec.objective = Objective::Kelly;
    spec.lower = -0.2;
    spec.upper = 0.7;
    const Golden golden{
        .status = SolveStatus::Optimal,
        .iterations = 3,
        .weights = {0x1.475900ba1b354p-4, 0x1.3d7b464f22ffcp-1, 0x1.6666666666665p-1,
                    -0x1.999999999999ap-3, -0x1.9999999999992p-3},
        .budget = kNaN,
        .target_return = kNaN,
        .volatility = kNaN,
        .leverage = kNaN,
        .bounds = {0x0p+0, 0x0p+0, -0x1.9c8170f855533p-8, 0x1.2de517e0bb3d2p-7, 0x1.443886dfd096fp-5},
        .residual = 0x0p+0,
    };
    expect_golden(spec, covariance_data(true, true), golden);
}

TEST(Determinism, TrackingErrorHasTheCommittedBits)
{
    Data data = covariance_data(true);
    ASSERT_TRUE(data.set_returns(kReturns, kBenchmark, kBenchmark.size()).has_value());
    Spec spec;
    spec.objective = Objective::TrackingError;
    spec.lower = 0.0;
    spec.upper = 0.6;
    spec.target_return = 0.09;
    const Golden golden{
        .status = SolveStatus::Optimal,
        .iterations = 2,
        .weights = {0x1.19052e9e84b45p-3, 0x1.72d5c2c501292p-2, 0x1.3923d49e3d9d7p-2, 0x1.8f07a29afdbe8p-3},
        .budget = -0x1.29de591dff227p-14,
        .target_return = 0x1.8590889f3e9e6p-10,
        .volatility = kNaN,
        .leverage = kNaN,
        .bounds = {0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0},
        .residual = 0x0p+0,
    };
    expect_golden(spec, data, golden);
}

TEST(Determinism, MaxSharpeHasTheCommittedBits)
{
    Spec spec;
    spec.objective = Objective::MaxSharpe;
    spec.lower = 0.0;
    spec.upper = 0.28;
    const Golden golden{
        .status = SolveStatus::Optimal,
        .iterations = 3,
        .weights = {0x1.1eb851eb851ecp-2, 0x1.d2c1a6751c8b5p-3, 0x1.b25d11dccef99p-3,
                    0x1.1eb851eb851ecp-2, 0x0p+0},
        .budget = kNaN,
        .target_return = kNaN,
        .volatility = kNaN,
        .leverage = kNaN,
        .bounds = {-0x1.576a1aa417aa4p-12, 0x0p+0, 0x0p+0, -0x1.0a539aa0b98cfp-6, -0x1.304ba73fd7e93p-8},
        .residual = 0x0p+0,
    };
    expect_golden(spec, covariance_data(true, true), golden);
}

TEST(Determinism, MaxDiversificationHasTheCommittedBits)
{
    Spec spec;
    spec.objective = Objective::MaxDiversification;
    spec.lower = -0.2;
    spec.upper = 0.4;
    spec.leverage = 1.1;
    const Golden golden{
        .status = SolveStatus::Optimal,
        .iterations = 2,
        .weights = {0x1.0d82fd5c8c9dbp-2, 0x1.74fbcae62eedcp-3, 0x1.3ccb072d84a38p-3, 0x1.9999999999999p-2},
        .budget = 0x1.7c19e2ad6c194p-8,
        .target_return = kNaN,
        .volatility = kNaN,
        .leverage = -0x0p+0,
        .bounds = {0x0p+0, 0x0p+0, 0x0p+0, -0x1.db205b58c71f8p-7},
        .residual = 0x0p+0,
    };
    expect_golden(spec, covariance_data(false), golden);
}

TEST(Determinism, RiskParityHasTheCommittedBits)
{
    Spec spec;
    spec.objective = Objective::RiskParity;
    spec.lower = 0.0;
    spec.risk_budget = {0.4, 0.3, 0.2, 0.1};
    spec.target_volatility = 0.12;
    const Golden golden{
        .status = SolveStatus::Optimal,
        .iterations = 21,
        .weights = {0x1.c417131eb7d22p-2, 0x1.1a128e75af2dp-2, 0x1.5d7650a8b83d3p-3, 0x1.cc6cd85cf3896p-4},
        .budget = -0x1.915aea403cb09p-4,
        .target_return = kNaN,
        .volatility = 0x1.a21409583f38p-1,
        .leverage = kNaN,
        .bounds = {0x0p+0, 0x0p+0, 0x0p+0, 0x0p+0},
        .residual = 0x1.839466f9e8ff6p-8,
    };
    expect_golden(spec, covariance_data(false), golden);
}

} // namespace
