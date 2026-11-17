#pragma once

#include <cstddef>
#include <cstdint>
#include <map>
#include <random>
#include <tuple>
#include <vector>

#include "core/data.hpp"
#include "random_covariance.hpp"

namespace po::bench
{

// Every instance is determined by its size and these constants, so two runs (and two builds) see the same
// numbers. A covariance of n = 1000 takes seconds to generate: each one is generated once per process and kept.
inline constexpr std::uint64_t kSeed = 1;
inline constexpr double kCondition = 1e2;
inline constexpr double kRiskFreeRate = 0.01;
inline constexpr double kMeanLow = 0.02;
inline constexpr double kMeanHigh = 0.12;

struct DataOptions
{
    bool cash = false;
    // Observations of the return series of TRACKING_ERROR; 0 means no series.
    std::size_t observations = 0;
};

inline const std::vector<double>& covariance_of(std::size_t n)
{
    static std::map<std::size_t, std::vector<double>> cache;
    auto found = cache.find(n);
    if (found == cache.end())
    {
        found = cache.emplace(n, po::testing::random_covariance(n, kSeed, kCondition)).first;
    }
    return found->second;
}

inline std::vector<double> mean_of(std::size_t n)
{
    return po::testing::random_vector(n, kSeed + 1, kMeanLow, kMeanHigh);
}

// t x n returns of scale 2% and a benchmark near their cross-sectional mean, row-major.
inline std::vector<double> returns_of(std::size_t observations, std::size_t n)
{
    std::mt19937_64 engine(kSeed + 2);
    std::vector<double> returns(observations * n);
    for (double& value : returns)
    {
        value = 0.02 * po::testing::standard_normal(engine);
    }
    return returns;
}

inline std::vector<double> benchmark_of(const std::vector<double>& returns, std::size_t observations,
                                        std::size_t n)
{
    std::mt19937_64 engine(kSeed + 3);
    std::vector<double> benchmark(observations);
    for (std::size_t t = 0; t < observations; ++t)
    {
        double sum = 0.0;
        for (std::size_t j = 0; j < n; ++j)
        {
            sum += returns[t * n + j];
        }
        benchmark[t] = sum / static_cast<double>(n) + 0.005 * po::testing::standard_normal(engine);
    }
    return benchmark;
}

inline po::Data build_data(std::size_t n, const DataOptions& options)
{
    po::Data data = po::Data::create(n, mean_of(n), covariance_of(n)).value();
    if (options.cash)
    {
        data.set_risk_free(kRiskFreeRate);
    }
    if (options.observations > 0)
    {
        std::vector<double> returns = returns_of(options.observations, n);
        std::vector<double> benchmark = benchmark_of(returns, options.observations, n);
        data.set_returns(returns, benchmark, options.observations).value();
    }
    return data;
}

// The data of a size and options, built once; the reference stays valid for the whole process.
inline const po::Data& data_of(std::size_t n, DataOptions options)
{
    static std::map<std::tuple<std::size_t, bool, std::size_t>, po::Data> cache;
    auto key = std::make_tuple(n, options.cash, options.observations);
    auto found = cache.find(key);
    if (found == cache.end())
    {
        found = cache.emplace(key, build_data(n, options)).first;
    }
    return found->second;
}

} // namespace po::bench
