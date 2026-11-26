// BM_Kernel/<kernel>/<n>: one numerical kernel of the library on the covariance of size n (cond 1e2).

#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "estimation/moments.hpp"
#include "instances.hpp"
#include "linalg/cholesky.hpp"
#include "linalg/qr.hpp"
#include "qp/active_set.hpp"
#include "qp/goldfarb_idnani.hpp"
#include "qp/problem.hpp"
#include "random_covariance.hpp"
#include "riskparity/objective.hpp"

namespace
{

constexpr std::uint64_t kVectorSeed = 7;

std::size_t size_of(const benchmark::State& state)
{
    return static_cast<std::size_t>(state.range(0));
}

const po::Data& market(const benchmark::State& state)
{
    return po::bench::data_of(size_of(state), po::bench::DataOptions{});
}

void cholesky(benchmark::State& state)
{
    std::size_t n = size_of(state);
    std::span<const double> covariance = market(state).covariance();
    if (!po::linalg::cholesky(covariance, n).has_value())
    {
        state.SkipWithError("the covariance is not positive definite");
        return;
    }
    for (auto _ : state)
    {
        auto factor = po::linalg::cholesky(covariance, n);
        benchmark::DoNotOptimize(factor);
    }
}

void solve_lower(benchmark::State& state)
{
    std::size_t n = size_of(state);
    std::span<const double> factor = market(state).covariance_factor();
    std::vector<double> right_hand_side = po::testing::random_vector(n, kVectorSeed, -1.0, 1.0);
    std::vector<double> work(n);
    for (auto _ : state)
    {
        std::ranges::copy(right_hand_side, work.begin());
        po::linalg::solve_lower(factor, n, work);
        benchmark::DoNotOptimize(work.data());
        benchmark::ClobberMemory();
    }
}

// A tall matrix of 2n rows, the shape of the return series of TRACKING_ERROR.
void qr(benchmark::State& state)
{
    std::size_t n = size_of(state);
    std::size_t rows = 2 * n;
    std::vector<double> matrix = po::bench::returns_of(rows, n);
    if (!po::linalg::has_full_column_rank(po::linalg::householder_qr(matrix, rows, n)))
    {
        state.SkipWithError("the matrix is rank deficient");
        return;
    }
    for (auto _ : state)
    {
        po::linalg::Qr factorization = po::linalg::householder_qr(matrix, rows, n);
        benchmark::DoNotOptimize(factorization.upper.data());
    }
}

// A long-only QP over the covariance with a linear term, solved by po::qp::solve: the time of one solve, and
// the iterations it took.
po::qp::Problem long_only_problem(const po::Data& data)
{
    std::size_t n = data.asset_count();
    po::qp::Problem problem;
    problem.n = n;
    problem.hessian.assign(data.covariance().begin(), data.covariance().end());
    problem.hessian_factor.assign(data.covariance_factor().begin(), data.covariance_factor().end());
    problem.linear = po::bench::mean_of(n);
    for (double& value : problem.linear)
    {
        value = -value;
    }
    problem.equality_rows.assign(n, 1.0);
    problem.equality_rhs = {1.0};
    problem.inequality_rows.assign(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        problem.inequality_rows[i * n + i] = 1.0;
    }
    problem.inequality_rhs.assign(n, 0.0);
    return problem;
}

void qp_iteration(benchmark::State& state)
{
    po::qp::Problem problem = long_only_problem(market(state));
    po::qp::Solution probe = po::qp::solve(problem);
    if (probe.status != po::qp::Status::Optimal)
    {
        state.SkipWithError("the long-only QP is not optimal");
        return;
    }
    for (auto _ : state)
    {
        po::qp::Solution solution = po::qp::solve(problem);
        benchmark::DoNotOptimize(solution.x.data());
        benchmark::ClobberMemory();
    }
    state.counters["qp_iterations"] = probe.iterations;
    auto per_iteration = benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert;
    state.counters["seconds_per_qp_iteration"] =
        benchmark::Counter(probe.iterations, static_cast<benchmark::Counter::Flags>(per_iteration));
}

void compensated_dot(benchmark::State& state)
{
    std::size_t n = size_of(state);
    std::vector<double> a = po::testing::random_vector(n, kVectorSeed, -1.0, 1.0);
    std::vector<double> b = po::testing::random_vector(n, kVectorSeed + 1, -1.0, 1.0);
    for (auto _ : state)
    {
        double value = po::riskparity::compensated_dot(a, b);
        benchmark::DoNotOptimize(value);
    }
}

void covariance_times(benchmark::State& state)
{
    const po::Data& data = market(state);
    std::vector<double> x = po::testing::random_vector(size_of(state), kVectorSeed, -1.0, 1.0);
    for (auto _ : state)
    {
        std::vector<double> product = po::riskparity::covariance_times(data, x);
        benchmark::DoNotOptimize(product.data());
    }
}

void gauss_newton(benchmark::State& state)
{
    std::size_t n = size_of(state);
    std::vector<double> budget(n, 1.0 / static_cast<double>(n));
    po::riskparity::Linearization point = po::riskparity::linearize(market(state), budget, budget, false);
    for (auto _ : state)
    {
        std::vector<double> curvature = po::riskparity::gauss_newton(point.jacobian, n);
        benchmark::DoNotOptimize(curvature.data());
    }
}

void empty_active_set(benchmark::State& state)
{
    std::size_t n = size_of(state);
    std::span<const double> factor = market(state).covariance_factor();
    for (auto _ : state)
    {
        po::qp::ActiveSet set = po::qp::empty_active_set(factor, n);
        benchmark::DoNotOptimize(set.basis.data());
    }
}

// Mean and covariance of 2n observations, the panel shape of the TRACKING_ERROR series.
void sample_moments(benchmark::State& state)
{
    std::size_t n = size_of(state);
    std::size_t t = 2 * n;
    std::vector<double> returns = po::bench::returns_of(t, n);
    if (!po::estimation::sample_moments(returns, t, n).has_value())
    {
        state.SkipWithError("the return panel is rejected");
        return;
    }
    for (auto _ : state)
    {
        auto moments = po::estimation::sample_moments(returns, t, n);
        benchmark::DoNotOptimize(moments);
    }
}

void kernel_sizes(benchmark::internal::Benchmark* kernel)
{
    for (int n : {10, 100, 300, 1000})
    {
        kernel->Arg(n);
    }
    kernel->Unit(benchmark::kMicrosecond);
}

BENCHMARK(cholesky)->Name("BM_Kernel/cholesky")->Apply(kernel_sizes);
BENCHMARK(solve_lower)->Name("BM_Kernel/solve_lower")->Apply(kernel_sizes);
BENCHMARK(qr)->Name("BM_Kernel/qr")->Apply(kernel_sizes);
BENCHMARK(qp_iteration)->Name("BM_Kernel/qp_iteration")->Apply(kernel_sizes);
BENCHMARK(compensated_dot)->Name("BM_Kernel/compensated_dot")->Apply(kernel_sizes);
BENCHMARK(covariance_times)->Name("BM_Kernel/covariance_times")->Apply(kernel_sizes);
BENCHMARK(gauss_newton)->Name("BM_Kernel/gauss_newton")->Apply(kernel_sizes);
BENCHMARK(empty_active_set)->Name("BM_Kernel/empty_active_set")->Apply(kernel_sizes);
BENCHMARK(sample_moments)->Name("BM_Kernel/sample_moments")->Apply(kernel_sizes);

} // namespace
