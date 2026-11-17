// Runs every case of the benchmark suite (objective x variant x size) once through po::solve and prints, per case,
// the status, the solver iterations, the time and the allocations; with --json it also dumps the solutions in
// hex-float, for `quick.py compare-solutions`. See docs/performance/measuring.md.
//
// portopt_profile [--filter TEXT] [--max-n N] [--repeat K] [--json PATH]

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <format>
#include <fstream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "alloc_counter.hpp"
#include "cases.hpp"
#include "solve/solve.hpp"

namespace
{

using po::bench::Case;

constexpr double kDefaultRiskAversion = 5.0;
constexpr std::size_t kDefaultMaxSize = 100;

struct Options
{
    std::string filter;
    std::size_t max_size = kDefaultMaxSize;
    int repeat = 3;
    std::string json;
};

struct Record
{
    std::string name;
    po::Result result;
    double objective = 0.0;
    double seconds = 0.0;
    std::uint64_t allocations = 0;
};

void apply_flag(Options& options, std::string_view flag, const std::string& value)
{
    if (flag == "--filter")
    {
        options.filter = value;
        return;
    }
    if (flag == "--max-n")
    {
        options.max_size = std::stoul(value);
        return;
    }
    if (flag == "--repeat")
    {
        options.repeat = std::stoi(value);
        return;
    }
    if (flag == "--json")
    {
        options.json = value;
        return;
    }
    throw std::invalid_argument("unknown flag: " + std::string(flag));
}

Options parse(std::span<char*> args)
{
    Options options;
    for (std::size_t index = 1; index < args.size(); index += 2)
    {
        if (index + 1 >= args.size())
        {
            throw std::invalid_argument(std::string(args[index]) + " needs a value");
        }
        apply_flag(options, args[index], args[index + 1]);
    }
    if (options.repeat < 1)
    {
        throw std::invalid_argument("--repeat must be at least 1");
    }
    return options;
}

double objective_value(const Case& instance, const po::Result& result)
{
    const po::Data& data = *instance.data;
    switch (instance.spec.objective)
    {
    case po::Objective::MinVariance:
        return std::pow(po::metrics::portfolio_volatility(data, result.weights), 2);
    case po::Objective::MaxReturn:
        return po::metrics::portfolio_return(data, result.weights);
    case po::Objective::MaxSharpe:
        return po::metrics::sharpe_ratio(data, result.weights);
    case po::Objective::Utility:
        return po::metrics::utility(data, result.weights, instance.spec.risk_aversion.value_or(kDefaultRiskAversion));
    case po::Objective::TrackingError:
        return po::metrics::tracking_error(data, result.weights);
    case po::Objective::MaxDiversification:
        return po::metrics::diversification_ratio(data, result.weights);
    case po::Objective::RiskParity:
        return result.residual;
    case po::Objective::Importance:
    case po::Objective::Kelly:
        break;
    }
    throw std::logic_error("no objective value for " + instance.name);
}

double seconds_of(const Case& instance)
{
    auto started = std::chrono::steady_clock::now();
    auto result = po::solve(instance.spec, *instance.data);
    auto finished = std::chrono::steady_clock::now();
    if (!result.has_value())
    {
        throw std::runtime_error(instance.name + ": " + result.error().message);
    }
    return std::chrono::duration<double>(finished - started).count();
}

Record run_case(const Case& instance, int repeat)
{
    Record record;
    record.name = instance.name;
    std::uint64_t before = po::bench::allocation_count();
    auto started = std::chrono::steady_clock::now();
    auto first = po::solve(instance.spec, *instance.data);
    auto finished = std::chrono::steady_clock::now();
    record.allocations = po::bench::allocation_count() - before;
    if (!first.has_value())
    {
        throw std::runtime_error(instance.name + ": " + first.error().message);
    }
    record.result = std::move(*first);
    record.objective = objective_value(instance, record.result);
    record.seconds = std::chrono::duration<double>(finished - started).count();
    for (int run = 1; run < repeat; ++run)
    {
        record.seconds = std::min(record.seconds, seconds_of(instance));
    }
    return record;
}

bool selected(const Options& options, const std::string& name, std::size_t n)
{
    return n <= options.max_size && name.find(options.filter) != std::string::npos;
}

std::vector<Record> run_all(const Options& options)
{
    std::vector<Record> records;
    for (const po::bench::Family& family : po::bench::families())
    {
        for (std::size_t n : po::bench::sizes_of(family))
        {
            if (!selected(options, po::bench::case_name(family, n), n))
            {
                continue;
            }
            records.push_back(run_case(po::bench::case_of(family, n), options.repeat));
            const Record& record = records.back();
            std::printf("%-44s status %d  iterations %6d  %10.3f ms  %8llu allocations\n", record.name.c_str(),
                        static_cast<int>(record.result.status), record.result.iterations, record.seconds * 1e3,
                        static_cast<unsigned long long>(record.allocations));
            std::fflush(stdout);
        }
    }
    return records;
}

std::string hex_list(std::span<const double> values)
{
    std::string text = "[";
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        text += std::format("{}\"{:a}\"", i == 0 ? "" : ",", values[i]);
    }
    return text + "]";
}

std::string multipliers_of(const po::Result& result)
{
    return std::format("{{\"budget\":\"{:a}\",\"target_return\":\"{:a}\",\"volatility\":\"{:a}\",\"leverage\":\"{:a}\"}}",
                       result.budget_multiplier, result.target_return_multiplier, result.volatility_multiplier,
                       result.leverage_multiplier);
}

std::string json_of(const Record& record)
{
    const po::Result& result = record.result;
    return std::format("{{\"name\":\"{}\",\"status\":{},\"iterations\":{},\"seconds\":{:.9g},\"allocations\":{},"
                       "\"objective\":\"{:a}\",\"residual\":\"{:a}\",\"multipliers\":{},\"bound_multipliers\":{},"
                       "\"weights\":{}}}",
                       record.name, static_cast<int>(result.status), result.iterations, record.seconds,
                       record.allocations, record.objective, result.residual, multipliers_of(result),
                       hex_list(result.bound_multipliers), hex_list(result.weights));
}

void write_json(const std::string& path, const std::vector<Record>& records)
{
    std::ofstream file(path);
    file << "{\"format\":\"portopt-solutions-1\",\"solutions\":[\n";
    for (std::size_t i = 0; i < records.size(); ++i)
    {
        file << json_of(records[i]) << (i + 1 < records.size() ? ",\n" : "\n");
    }
    file << "]}\n";
    if (!file)
    {
        throw std::runtime_error("cannot write " + path);
    }
}

void print_total(const std::vector<Record>& records)
{
    double seconds = 0.0;
    std::uint64_t allocations = 0;
    for (const Record& record : records)
    {
        seconds += record.seconds;
        allocations += record.allocations;
    }
    std::printf("%zu cases, %.3f s of solves, %llu allocations\n", records.size(), seconds,
                static_cast<unsigned long long>(allocations));
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        Options options = parse(std::span<char*>(argv, static_cast<std::size_t>(argc)));
        std::vector<Record> records = run_all(options);
        if (records.empty())
        {
            throw std::invalid_argument("--filter and --max-n select no case");
        }
        print_total(records);
        if (!options.json.empty())
        {
            write_json(options.json, records);
        }
    }
    catch (const std::exception& error)
    {
        std::fprintf(stderr, "portopt_profile: %s\n", error.what());
        return 1;
    }
    return 0;
}
