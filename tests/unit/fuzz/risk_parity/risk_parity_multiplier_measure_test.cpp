#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "fuzz/frontier_reference.hpp"
#include "fuzz/risk_parity/risk_parity_reference.hpp"
#include "fuzz/risk_parity/risk_parity_refinement.hpp"
#include "linalg/qr.hpp"
#include "riskparity/kkt.hpp"
#include "solve/solve.hpp"

namespace
{

using po::Result;
using po::SolveStatus;
using Shape = po::testing::RiskParityShape;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// [H, -A'; A, 0] and -[stationarity; gaps] over every weight (no box), as newton_polish builds its Newton
// system: each row of A and its right-hand side scaled to the largest entry of H, so that the unknown z_k of
// equality k is dlambda_k / scales[k].
struct ScaledSystem
{
    std::size_t size;
    std::vector<double> matrix;
    std::vector<double> rhs;
    std::vector<double> scales;
};

ScaledSystem scaled_system(const po::riskparity::KktEvaluation& evaluation)
{
    std::size_t n = evaluation.stationarity.size();
    std::size_t m = evaluation.gaps.size();
    std::size_t size = n + m;
    ScaledSystem system{size, std::vector<double>(size * size, 0.0), std::vector<double>(size),
                        std::vector<double>(m)};
    for (std::size_t a = 0; a < n; ++a)
    {
        for (std::size_t b = 0; b < n; ++b)
        {
            system.matrix[a * size + b] = evaluation.lagrangian_hessian[a * n + b];
        }
        system.rhs[a] = -evaluation.stationarity[a];
    }
    double hessian_scale = po::testing::largest_magnitude(system.matrix);
    for (std::size_t k = 0; k < m; ++k)
    {
        std::span<const double> row = std::span<const double>(evaluation.rows).subspan(k * n, n);
        system.scales[k] = hessian_scale / po::testing::largest_magnitude(row);
        for (std::size_t a = 0; a < n; ++a)
        {
            system.matrix[a * size + n + k] = -system.scales[k] * row[a];
            system.matrix[(n + k) * size + a] = system.scales[k] * row[a];
        }
        system.rhs[n + k] = -system.scales[k] * evaluation.gaps[k];
    }
    return system;
}

// lambda + dlambda of one Newton-KKT solve at the point, by a Householder QR of the scaled system; the
// weights stay where they are. NaN when the system is rank deficient to working precision.
std::vector<double> newton_multipliers(const po::riskparity::KktProblem& problem,
                                       const po::riskparity::KktPoint& point)
{
    std::vector<double> multipliers = point.equality_multipliers;
    ScaledSystem system = scaled_system(po::riskparity::evaluate_kkt(problem, point));
    po::linalg::Qr qr = po::linalg::householder_qr(system.matrix, system.size, system.size);
    if (!po::linalg::has_full_column_rank(qr))
    {
        multipliers.assign(multipliers.size(), kNaN);
        return multipliers;
    }
    std::vector<double> step = po::linalg::least_squares(qr, std::move(system.rhs));
    std::size_t n = point.weights.size();
    for (std::size_t k = 0; k < multipliers.size(); ++k)
    {
        multipliers[k] += system.scales[k] * step[n + k];
    }
    return multipliers;
}

double norm_of(std::span<const double> values)
{
    double total = 0.0;
    for (double value : values)
    {
        total += value * value;
    }
    return std::sqrt(total);
}

std::vector<double> difference(std::span<const double> a, std::span<const double> b)
{
    std::vector<double> result(a.size());
    for (std::size_t k = 0; k < a.size(); ++k)
    {
        result[k] = a[k] - b[k];
    }
    return result;
}

// The KKT point refined from the result in long double, by the steps of forward_error.
po::testing::WideKkt refined(const po::Spec& spec, const po::Data& data, const Result& result)
{
    constexpr int kSteps = 6;
    std::vector<double> budget = po::testing::budget_of(spec, data.asset_count());
    po::testing::WideKkt kkt = po::testing::wide_kkt_of(spec, result);
    for (int step = 0; step < kSteps; ++step)
    {
        if (!po::testing::wide_newton_step(spec, data, budget, kkt))
        {
            kkt.multipliers.assign(kkt.multipliers.size(), static_cast<long double>(kNaN));
            return kkt;
        }
    }
    return kkt;
}

// The budget and volatility multipliers of a result without a target return.
Result with_multipliers(const Result& start, std::span<const double> multipliers)
{
    Result result = start;
    result.budget_multiplier = multipliers.front();
    result.volatility_multiplier = multipliers.back();
    return result;
}

// lambda_vol relative to |lambda_vol*|; lambda_budget relative to max(|lambda_budget*|, |lambda_vol*| s), the
// size of the volatility term of the stationarity (w'grad vol = vol).
struct MultiplierError
{
    double volatility;
    double budget;
};

MultiplierError multiplier_error(double target, const std::vector<long double>& reference,
                                 std::span<const double> multipliers)
{
    long double volatility = reference.back();
    long double budget_scale = std::fmax(std::fabs(reference.front()), std::fabs(volatility) * target);
    return MultiplierError{
        static_cast<double>(std::fabs(multipliers.back() - volatility) / std::fabs(volatility)),
        static_cast<double>(std::fabs(multipliers.front() - reference.front()) / budget_scale)};
}

// One case of the measurement: the result's least-squares multipliers lambda_LS, the Newton ones
// lambda_LS + dlambda, and what each does to the stationarity at the weights of the result.
struct Measurement
{
    double weights;
    double reference_over_floor;
    MultiplierError least_squares;
    MultiplierError newton;
    // forward_error.multipliers of each: max_k |lambda_k - lambda*_k| / max_k |lambda*_k|.
    double least_squares_combined;
    double newton_combined;
    // |dlambda| / |lambda_LS|, and |dlambda_2| / |dlambda| for a second solve from (w, lambda_LS + dlambda).
    double change;
    double contraction;
    // risk_parity_stationarity(...).over_floor: the gate of expect_risk_parity_kkt is 2.
    double least_squares_over_floor;
    double newton_over_floor;
    // residual_of(...).stationarity, in units of the noise of the library: the gate of is_kkt_point is 16.
    double least_squares_over_noise;
    double newton_over_noise;
};

// NaN when a multiplier is: risk_parity_stationarity drops a NaN multiplier's term instead.
double over_floor(const po::Spec& spec, const po::Data& data, const Result& result)
{
    if (std::isnan(result.budget_multiplier) || std::isnan(result.volatility_multiplier))
    {
        return kNaN;
    }
    return po::testing::risk_parity_stationarity(spec, data, result, 0.0).over_floor;
}

double over_noise(const po::riskparity::KktProblem& problem, const po::riskparity::KktPoint& point)
{
    po::riskparity::KktEvaluation evaluation = po::riskparity::evaluate_kkt(problem, point);
    return po::riskparity::residual_of(problem, point, evaluation).stationarity;
}

struct Point
{
    po::riskparity::KktPoint least_squares;
    po::riskparity::KktPoint newton;
    std::vector<long double> reference;
};

Point point_of(const po::testing::GeneratedRiskParity& generated, const Result& result,
               const po::riskparity::KktProblem& problem)
{
    po::riskparity::KktPoint least_squares{result.weights,
                                           {result.budget_multiplier, result.volatility_multiplier},
                                           result.bound_multipliers};
    po::riskparity::KktPoint newton = least_squares;
    newton.equality_multipliers = newton_multipliers(problem, least_squares);
    return Point{std::move(least_squares), std::move(newton),
                 refined(generated.spec, generated.data, result).multipliers};
}

Measurement measure(const po::testing::GeneratedRiskParity& generated, const Result& result,
                    const Point& point, const po::riskparity::KktProblem& problem)
{
    const po::Spec& spec = generated.spec;
    const po::Data& data = generated.data;
    std::span<const double> least_squares = point.least_squares.equality_multipliers;
    std::span<const double> newton = point.newton.equality_multipliers;
    Result corrected = with_multipliers(result, newton);
    po::testing::ForwardError forward = po::testing::forward_error(spec, data, result);
    std::vector<double> step = difference(newton, least_squares);
    std::vector<double> second = difference(newton_multipliers(problem, point.newton), newton);
    double target = *spec.target_volatility;
    return Measurement{forward.weights,
                       forward.reference_over_floor,
                       multiplier_error(target, point.reference, least_squares),
                       multiplier_error(target, point.reference, newton),
                       forward.multipliers,
                       po::testing::forward_error(spec, data, corrected).multipliers,
                       norm_of(step) / norm_of(least_squares),
                       norm_of(second) / norm_of(step),
                       over_floor(spec, data, result),
                       over_floor(spec, data, corrected),
                       over_noise(problem, point.least_squares),
                       over_noise(problem, point.newton)};
}

void print_measurement(const char* label, const Shape& shape, const Measurement& measured)
{
    std::printf("| %s | %.0e | %zu | %llu | %.1e | %.1e / %.1e | %.1e / %.1e | %.1e / %.1e | %.1e | %.1e | "
                "%.2f / %.2f | %.2f / %.2f | %.1e |\n",
                label, shape.condition, shape.n, static_cast<unsigned long long>(shape.seed),
                measured.weights, measured.least_squares.volatility, measured.least_squares.budget,
                measured.newton.volatility, measured.newton.budget, measured.least_squares_combined,
                measured.newton_combined, measured.change, measured.contraction,
                measured.least_squares_over_floor, measured.newton_over_floor,
                measured.least_squares_over_noise, measured.newton_over_noise, measured.reference_over_floor);
}

void dump_values(std::FILE* file, const char* name, std::span<const double> values)
{
    std::fprintf(file, "%s", name);
    for (double value : values)
    {
        std::fprintf(file, " %.17g", value);
    }
    std::fprintf(file, "\n");
}

// The problem, the point and the analytic Lagrangian Hessian at (w, lambda_LS), for the 80-digit reference.
void dump(const char* label, const Shape& shape, const po::testing::GeneratedRiskParity& generated,
          const po::riskparity::KktProblem& problem, const Point& point)
{
    const char* path = std::getenv("PORTOPT_MULTIPLIER_DUMP");
    if (path == nullptr)
    {
        return;
    }
    std::FILE* file = std::fopen(path, "a");
    ASSERT_NE(file, nullptr) << path;
    std::fprintf(file, "case %s %.0e %zu %llu\n", label, shape.condition, shape.n,
                 static_cast<unsigned long long>(shape.seed));
    std::fprintf(file, "target %.17g\n", *generated.spec.target_volatility);
    dump_values(file, "budget", problem.budget);
    dump_values(file, "covariance", generated.data.covariance());
    dump_values(file, "weights", point.least_squares.weights);
    dump_values(file, "least_squares", point.least_squares.equality_multipliers);
    dump_values(file, "newton", point.newton.equality_multipliers);
    std::fprintf(file, "reference %.21Lg %.21Lg\n", point.reference.front(), point.reference.back());
    po::riskparity::KktEvaluation evaluation = po::riskparity::evaluate_kkt(problem, point.least_squares);
    dump_values(file, "hessian", evaluation.lagrangian_hessian);
    std::fclose(file);
}

void measure_case(const char* label, const Shape& shape, const po::testing::GeneratedRiskParity& generated)
{
    SCOPED_TRACE(po::testing::describe(shape));
    Result result = po::solve(generated.spec, generated.data).value();
    if (result.status != SolveStatus::Optimal)
    {
        std::printf("| %s | %.0e | %zu | %llu | status %d |\n", label, shape.condition, shape.n,
                    static_cast<unsigned long long>(shape.seed), static_cast<int>(result.status));
        return;
    }
    po::riskparity::KktProblem problem{generated.spec, generated.data, generated.spec.risk_budget};
    Point point = point_of(generated, result, problem);
    print_measurement(label, shape, measure(generated, result, point, problem));
    dump(label, shape, generated, problem, point);
}

// The cases of the 80-digit probe of docs/research/2026-11-16-multiplier-accuracy.md: s = 1.2 sigma_min.
po::testing::GeneratedRiskParity probe_problem(const Shape& shape)
{
    po::testing::GeneratedRiskParity generated = po::testing::generate(shape);
    po::Result minimum = po::testing::frontier_point(generated.spec, generated.data, std::nullopt);
    EXPECT_EQ(minimum.status, SolveStatus::Optimal);
    double smallest = po::testing::portfolio_volatility(generated.data, minimum.weights);
    generated.spec.target_volatility = 1.2 * smallest;
    return generated;
}

// The multipliers of the target-volatility path at the point po::solve returns (no box, equal budget), by
// cond(Sigma): those of the result (least squares at the rounded weights) and those of one Newton-KKT solve
// there, both against the KKT point refined in long double; the measurement of docs/guide/validation.md. Run
// with --gtest_also_run_disabled_tests. PORTOPT_MULTIPLIER_DUMP=path appends every case to path for the
// 80-digit reference of the same section.
TEST(RiskParityMultiplierMeasure, DISABLED_MeasureNewtonMultipliers)
{
    if (!po::testing::kLongDoubleIsExtended)
    {
        GTEST_SKIP() << "long double is double here: the refinement has no extra precision";
    }
    std::printf("| cases | cond | n | seed | w error | LS vol / budget | Newton vol / budget | "
                "combined LS / Newton | dlambda | contraction | over floor LS / Newton | "
                "over noise LS / Newton | reference over floor |\n");
    for (double condition : {1e4, 1e6, 1e8, 1e10})
    {
        for (auto [n, seed] : {std::pair<std::size_t, std::uint64_t>{5, 7}, {5, 11}, {8, 3}})
        {
            Shape shape{n, condition, 1.0, Shape::Budget::Equal, Shape::Box::None, false, true, seed};
            measure_case("probe", shape, probe_problem(shape));
        }
        for (std::size_t n : {5u, 8u})
        {
            for (std::uint64_t seed = 402; seed <= 405; ++seed)
            {
                Shape shape{n, condition, 1.0, Shape::Budget::Equal, Shape::Box::None, false, true, seed};
                measure_case("grid", shape, po::testing::generate(shape));
            }
        }
    }
}

} // namespace
