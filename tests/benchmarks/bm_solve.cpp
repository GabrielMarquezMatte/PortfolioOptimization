// BM_Solve/<objective>/<variant>/<n>: po::solve end to end on the cases of cases.hpp.

#include <benchmark/benchmark.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <format>
#include <string>

#include "alloc_counter.hpp"
#include "cases.hpp"
#include "solve/solve.hpp"

namespace
{

using po::bench::Case;
using po::bench::Family;

constexpr std::size_t kSlowFrom = 1000;

// Empty when the solve ended at a finite optimal point; the harness that measures a solve that failed publishes
// the time of the failure.
std::string problem_of(const std::expected<po::Result, po::Failure>& result)
{
    if (!result.has_value())
    {
        return result.error().message;
    }
    if (result->status != po::SolveStatus::Optimal)
    {
        return std::format("status {} is not OPTIMAL", static_cast<int>(result->status));
    }
    for (double weight : result->weights)
    {
        if (!std::isfinite(weight))
        {
            return "a weight is not finite";
        }
    }
    return {};
}

void run_solve(benchmark::State& state, const Family& family)
{
    const Case& instance = po::bench::case_of(family, static_cast<std::size_t>(state.range(0)));
    auto probe = po::solve(instance.spec, *instance.data);
    std::string problem = problem_of(probe);
    if (!problem.empty())
    {
        state.SkipWithError((instance.name + ": " + problem).c_str());
        return;
    }
    std::uint64_t before = po::bench::allocation_count();
    for (auto _ : state)
    {
        auto result = po::solve(instance.spec, *instance.data);
        benchmark::DoNotOptimize(result);
        if (!result.has_value() || result->status != po::SolveStatus::Optimal)
        {
            state.SkipWithError((instance.name + ": a solve did not end OPTIMAL").c_str());
            break;
        }
        benchmark::ClobberMemory();
    }
    std::uint64_t allocations = po::bench::allocation_count() - before;
    state.counters["iterations"] = probe->iterations;
    state.counters["status"] = static_cast<int>(probe->status);
    state.counters["allocations"] = static_cast<double>(allocations) / static_cast<double>(state.iterations());
    if (family.objective == po::Objective::RiskParity)
    {
        state.counters["residual"] = probe->residual;
    }
}

bool register_solve_benchmarks()
{
    for (const Family& family : po::bench::families())
    {
        std::string name = std::format("BM_Solve/{}/{}", po::bench::objective_name(family.objective),
                                       po::bench::variant_name(family.variant));
        for (std::size_t n : po::bench::sizes_of(family))
        {
            benchmark::internal::Benchmark* registered = benchmark::RegisterBenchmark(
                name, [family](benchmark::State& state) { run_solve(state, family); });
            registered->Arg(static_cast<std::int64_t>(n))->Unit(benchmark::kMicrosecond);
            if (n >= kSlowFrom)
            {
                registered->Iterations(1);
            }
        }
    }
    return true;
}

[[maybe_unused]] const bool registered = register_solve_benchmarks();

} // namespace
