// BM_CApi/<what>/<n>: the fixed cost of the C interface against the C++ calls under it. The solve cases are the
// min_variance ones of cases.hpp, whose solves are the shortest of their size, so the ratio is the worst one.

#include <benchmark/benchmark.h>

#include <cstddef>
#include <cstdint>
#include <exception>
#include <string>
#include <utility>

#include "capi/capi_case.hpp"
#include "cases.hpp"
#include "solve/solve.hpp"

namespace
{

using po::bench::DataHandle;
using po::bench::Family;
using po::bench::ResultHandle;
using po::bench::SpecHandle;

constexpr std::int64_t kSizes[] = {10, 100, 1000};
constexpr Family kFree{po::Objective::MinVariance, po::bench::Variant::Free};
constexpr Family kLongOnly{po::Objective::MinVariance, po::bench::Variant::LongOnly};

std::size_t size_of(const benchmark::State& state)
{
    return static_cast<std::size_t>(state.range(0));
}

// Runs one benchmark body; an exception of the C helpers ends the benchmark with its message.
template <typename Body>
void guarded(benchmark::State& state, Body&& body)
{
    try
    {
        body();
    }
    catch (const std::exception& error)
    {
        state.SkipWithError(error.what());
    }
}

void call_abi_version(benchmark::State& state)
{
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(po_abi_version());
    }
}

// The cheapest call that goes through the guard (clear_error, the try block, the handle check).
void call_guarded(benchmark::State& state)
{
    guarded(state,
            [&]
            {
                DataHandle data = po::bench::make_c_data(po::bench::data_of(10, po::bench::DataOptions{}));
                std::size_t count = 0;
                for (auto _ : state)
                {
                    benchmark::DoNotOptimize(po_data_asset_count(data.get(), &count));
                }
            });
}

void spec_cycle_c(benchmark::State& state)
{
    guarded(state,
            [&]
            {
                const po::Spec& spec = po::bench::case_of(kLongOnly, 10).spec;
                for (auto _ : state)
                {
                    SpecHandle handle = po::bench::make_c_spec(spec);
                    benchmark::DoNotOptimize(handle.get());
                }
            });
}

void spec_cycle_cpp(benchmark::State& state)
{
    const po::Spec& spec = po::bench::case_of(kLongOnly, 10).spec;
    for (auto _ : state)
    {
        po::Spec copy = spec;
        benchmark::DoNotOptimize(copy);
    }
}

void data_cycle_c(benchmark::State& state)
{
    guarded(state,
            [&]
            {
                const po::Data& data = po::bench::data_of(size_of(state), po::bench::DataOptions{});
                for (auto _ : state)
                {
                    DataHandle handle = po::bench::make_c_data(data);
                    benchmark::DoNotOptimize(handle.get());
                }
            });
}

void data_cycle_cpp(benchmark::State& state)
{
    const po::Data& data = po::bench::data_of(size_of(state), po::bench::DataOptions{});
    for (auto _ : state)
    {
        auto created = po::Data::create(data.asset_count(), data.mu(), data.covariance());
        benchmark::DoNotOptimize(created);
    }
}

void run_solve_c(benchmark::State& state, const Family& family)
{
    guarded(state,
            [&]
            {
                const po::bench::Case& instance = po::bench::case_of(family, size_of(state));
                SpecHandle spec = po::bench::make_c_spec(instance.spec);
                DataHandle data = po::bench::make_c_data(*instance.data);
                for (auto _ : state)
                {
                    ResultHandle result = po::bench::solve_with_handles(spec, data);
                    benchmark::DoNotOptimize(result.get());
                }
            });
}

void run_solve_cpp(benchmark::State& state, const Family& family)
{
    const po::bench::Case& instance = po::bench::case_of(family, size_of(state));
    for (auto _ : state)
    {
        auto result = po::solve(instance.spec, *instance.data);
        benchmark::DoNotOptimize(result);
    }
}

// Every reader of the C interface on one result, the copy of the weights and of the bound multipliers included.
void result_read_c(benchmark::State& state)
{
    guarded(state,
            [&]
            {
                const po::bench::Case& instance = po::bench::case_of(kLongOnly, size_of(state));
                SpecHandle spec = po::bench::make_c_spec(instance.spec);
                DataHandle data = po::bench::make_c_data(*instance.data);
                ResultHandle result = po::bench::solve_with_handles(spec, data);
                for (auto _ : state)
                {
                    po::Result read = po::bench::read_c_result(result.get());
                    benchmark::DoNotOptimize(read);
                }
            });
}

void full_cycle_c(benchmark::State& state)
{
    guarded(state,
            [&]
            {
                const po::bench::Case& instance = po::bench::case_of(kLongOnly, size_of(state));
                for (auto _ : state)
                {
                    po::Result read = po::bench::solve_through_c_api(instance.spec, *instance.data);
                    benchmark::DoNotOptimize(read);
                }
            });
}

void full_cycle_cpp(benchmark::State& state)
{
    const po::bench::Case& instance = po::bench::case_of(kLongOnly, size_of(state));
    const po::Data& source = *instance.data;
    for (auto _ : state)
    {
        po::Spec spec = instance.spec;
        auto data = po::Data::create(source.asset_count(), source.mu(), source.covariance());
        auto result = po::solve(spec, data.value());
        benchmark::DoNotOptimize(result);
    }
}

void register_sized(const std::string& name, void (*body)(benchmark::State&))
{
    benchmark::internal::Benchmark* registered = benchmark::RegisterBenchmark(name, body);
    registered->Unit(benchmark::kMicrosecond);
    for (std::int64_t n : kSizes)
    {
        registered->Arg(n);
    }
}

void register_solve(const std::string& name, const Family& family)
{
    for (auto [suffix, body] : {std::pair{"c", &run_solve_c}, std::pair{"cpp", &run_solve_cpp}})
    {
        benchmark::internal::Benchmark* registered = benchmark::RegisterBenchmark(
            name + "_" + suffix, [family, body](benchmark::State& state) { body(state, family); });
        registered->Unit(benchmark::kMicrosecond);
        for (std::int64_t n : kSizes)
        {
            registered->Arg(n);
        }
    }
}

bool register_capi_benchmarks()
{
    benchmark::RegisterBenchmark("BM_CApi/call_abi_version", call_abi_version)->Unit(benchmark::kNanosecond);
    benchmark::RegisterBenchmark("BM_CApi/call_guarded", call_guarded)->Unit(benchmark::kNanosecond);
    benchmark::RegisterBenchmark("BM_CApi/spec_cycle_c", spec_cycle_c)->Unit(benchmark::kNanosecond);
    benchmark::RegisterBenchmark("BM_CApi/spec_cycle_cpp", spec_cycle_cpp)->Unit(benchmark::kNanosecond);
    register_sized("BM_CApi/data_cycle_c", data_cycle_c);
    register_sized("BM_CApi/data_cycle_cpp", data_cycle_cpp);
    register_solve("BM_CApi/solve_free", kFree);
    register_solve("BM_CApi/solve_long_only", kLongOnly);
    register_sized("BM_CApi/result_read_c", result_read_c);
    register_sized("BM_CApi/full_cycle_c", full_cycle_c);
    register_sized("BM_CApi/full_cycle_cpp", full_cycle_cpp);
    return true;
}

[[maybe_unused]] const bool registered = register_capi_benchmarks();

} // namespace
