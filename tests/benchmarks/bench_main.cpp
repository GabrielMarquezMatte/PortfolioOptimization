// Entry point of portopt_benchmarks: Google Benchmark's main plus `--selftest`, which checks the allocation counter.

#include <benchmark/benchmark.h>

#include <cstdio>
#include <string_view>

#include "alloc_counter.hpp"

int main(int argc, char** argv)
{
    if (argc == 2 && std::string_view(argv[1]) == "--selftest")
    {
        return po::bench::allocation_selftest() == 0 ? 0 : 1;
    }
    benchmark::Initialize(&argc, argv);
    if (benchmark::ReportUnrecognizedArguments(argc, argv))
    {
        return 1;
    }
    std::size_t ran = benchmark::RunSpecifiedBenchmarks();
    benchmark::Shutdown();
    if (ran == 0)
    {
        std::fprintf(stderr, "portopt_benchmarks: no benchmark matches the filter\n");
        return 1;
    }
    return 0;
}
