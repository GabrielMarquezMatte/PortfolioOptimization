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
#include "random_covariance.hpp"
#include "riskparity/kkt.hpp"
#include "riskparity/newton_polish.hpp"
#include "solve/solve.hpp"

namespace
{

using po::Result;
using po::SolveStatus;
using Shape = po::testing::RiskParityShape;

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

double norm_of(std::span<const double> values)
{
    double total = 0.0;
    for (double value : values)
    {
        total += value * value;
    }
    return std::sqrt(total);
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
    // equality_multipliers is [budget, volatility] here: no target return, so no middle entry.
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
    // The safeguard of corrected_multipliers kept the correction; the Newton columns are NaN otherwise.
    bool accepted;
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

// The least-squares point of a polish of the result (no step: the weights of the result), the correction
// corrected_multipliers made there, and the refined KKT point.
struct Point
{
    po::riskparity::KktPoint least_squares;
    po::riskparity::MultiplierCorrection correction;
    std::vector<long double> reference;
};

// The corrected multipliers when the safeguard kept them, NaN otherwise.
po::riskparity::KktPoint newton_point(const Point& point)
{
    po::riskparity::KktPoint newton = point.correction.reported_point;
    if (!point.correction.accepted)
    {
        newton.equality_multipliers.assign(newton.equality_multipliers.size(), kNaN);
    }
    return newton;
}

Point point_of(const po::testing::GeneratedRiskParity& generated, const Result& result,
               const po::riskparity::KktProblem& problem)
{
    po::riskparity::KktPoint start{result.weights,
                                   {result.budget_multiplier, result.volatility_multiplier},
                                   result.bound_multipliers};
    po::riskparity::PolishResult polished = po::riskparity::newton_polish(problem, start);
    EXPECT_EQ(polished.steps, 0);
    Result least_squares = with_multipliers(result, polished.least_squares_point.equality_multipliers);
    po::riskparity::MultiplierCorrection correction =
        po::riskparity::corrected_multipliers(problem, polished.least_squares_point);
    return Point{std::move(polished.least_squares_point), std::move(correction),
                 refined(generated.spec, generated.data, least_squares).multipliers};
}

Measurement measure(const po::testing::GeneratedRiskParity& generated, const Result& result,
                    const Point& point, const po::riskparity::KktProblem& problem)
{
    const po::Spec& spec = generated.spec;
    const po::Data& data = generated.data;
    po::riskparity::KktPoint newton_kkt = newton_point(point);
    std::span<const double> least_squares = point.least_squares.equality_multipliers;
    std::span<const double> newton = newton_kkt.equality_multipliers;
    Result fitted = with_multipliers(result, least_squares);
    Result corrected = with_multipliers(result, newton);
    po::testing::ForwardError forward = po::testing::forward_error(spec, data, fitted);
    double target = *spec.target_volatility;
    return Measurement{forward.weights,
                       forward.reference_over_floor,
                       multiplier_error(target, point.reference, least_squares),
                       multiplier_error(target, point.reference, newton),
                       forward.multipliers,
                       po::testing::forward_error(spec, data, corrected).multipliers,
                       point.correction.change / norm_of(least_squares),
                       point.correction.second_change / point.correction.change,
                       point.correction.accepted,
                       over_floor(spec, data, fitted),
                       over_floor(spec, data, corrected),
                       over_noise(problem, point.least_squares),
                       over_noise(problem, newton_kkt)};
}

void print_measurement(const char* label, const Shape& shape, const Measurement& measured)
{
    std::printf("| %s | %.0e | %zu | %llu | %.1e | %.1e / %.1e | %.1e / %.1e | %.1e / %.1e | %.1e | %.1e | "
                "%d | %.2f / %.2f | %.2f / %.2f | %.1e |\n",
                label, shape.condition, shape.n, static_cast<unsigned long long>(shape.seed),
                measured.weights, measured.least_squares.volatility, measured.least_squares.budget,
                measured.newton.volatility, measured.newton.budget, measured.least_squares_combined,
                measured.newton_combined, measured.change, measured.contraction,
                static_cast<int>(measured.accepted), measured.least_squares_over_floor,
                measured.newton_over_floor,
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
    dump_values(file, "newton", newton_point(point).equality_multipliers);
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
// cond(Sigma): the least-squares ones at the rounded weights and those of one Newton-KKT solve there
// (corrected_multipliers, NaN where its safeguard refused it), both against the KKT point refined in long
// double; the measurement of docs/guide/validation.md. Run with --gtest_also_run_disabled_tests.
// PORTOPT_MULTIPLIER_DUMP=path appends every case to path for the 80-digit reference of the same section.
TEST(RiskParityMultiplierMeasure, DISABLED_MeasureNewtonMultipliers)
{
    if (!po::testing::kLongDoubleIsExtended)
    {
        GTEST_SKIP() << "long double is double here: the refinement has no extra precision";
    }
    std::printf("| cases | cond | n | seed | w error | LS vol / budget | Newton vol / budget | "
                "combined LS / Newton | dlambda | contraction | accepted | over floor LS / Newton | "
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

// The KKT point of MIN_VARIANCE with a target return on the active set of a result, in long double: the
// weights whose bound multiplier is not 0 stay where the result has them, and the free weights w_F and the
// multipliers (budget, target return) solve
//   (2 Sigma w)_F = lambda_budget 1 + lambda_return mu_F,   1'w = 1,   mu'w = r,
// a linear system, by Newton steps (iterative refinement) from the result; beta of a held weight is what its
// row leaves. last_step is the largest change of a multiplier in the last step, relative to the multipliers.
struct QpReference
{
    std::vector<long double> weights;
    std::vector<long double> multipliers;
    std::vector<long double> bounds;
    long double last_step = 0.0L;
};

std::vector<long double> qp_stationarity(const po::Data& data, const QpReference& point)
{
    std::vector<long double> residual = po::testing::wide_covariance_times(data, point.weights);
    for (std::size_t j = 0; j < residual.size(); ++j)
    {
        residual[j] = 2.0L * residual[j] - point.multipliers[0] - point.multipliers[1] * data.mu()[j];
    }
    return residual;
}

// One Newton step on the free weights and the two multipliers; false when the system is singular.
bool qp_newton_step(const po::Data& data, double target, const std::vector<std::size_t>& free,
                    QpReference& point)
{
    std::size_t size = free.size() + 2;
    std::size_t n = point.weights.size();
    std::vector<long double> matrix(size * size, 0.0L);
    std::vector<long double> rhs(size, 0.0L);
    std::vector<long double> stationarity = qp_stationarity(data, point);
    std::vector<long double> mu = po::testing::widened(data.mu());
    for (std::size_t a = 0; a < free.size(); ++a)
    {
        for (std::size_t b = 0; b < free.size(); ++b)
        {
            matrix[a * size + b] = 2.0L * data.covariance()[free[a] * n + free[b]];
        }
        matrix[a * size + free.size()] = -1.0L;
        matrix[a * size + free.size() + 1] = -mu[free[a]];
        matrix[free.size() * size + a] = 1.0L;
        matrix[(free.size() + 1) * size + a] = mu[free[a]];
        rhs[a] = -stationarity[free[a]];
    }
    long double total = -1.0L;
    for (long double weight : point.weights)
    {
        total += weight;
    }
    rhs[free.size()] = -total;
    rhs[free.size() + 1] = -(po::testing::wide_dot(mu, point.weights) - target);
    if (!po::testing::solve_wide(matrix, rhs))
    {
        return false;
    }
    for (std::size_t a = 0; a < free.size(); ++a)
    {
        point.weights[free[a]] += rhs[a];
    }
    long double scale = std::fmax(std::fabs(point.multipliers[0]), std::fabs(point.multipliers[1]));
    point.multipliers[0] += rhs[free.size()];
    point.multipliers[1] += rhs[free.size() + 1];
    point.last_step = std::fmax(std::fabs(rhs[free.size()]), std::fabs(rhs[free.size() + 1])) / scale;
    return true;
}

std::optional<QpReference> qp_reference(const po::Spec& spec, const po::Data& data, const Result& result)
{
    constexpr int kSteps = 4;
    QpReference point{po::testing::widened(result.weights),
                      {result.budget_multiplier, result.target_return_multiplier},
                      {},
                      0.0L};
    std::vector<std::size_t> free;
    for (std::size_t j = 0; j < result.weights.size(); ++j)
    {
        if (result.bound_multipliers[j] == 0.0)
        {
            free.push_back(j);
        }
    }
    for (int step = 0; step < kSteps; ++step)
    {
        if (!qp_newton_step(data, *spec.target_return, free, point))
        {
            return std::nullopt;
        }
    }
    point.bounds = qp_stationarity(data, point);
    return point;
}

// One result of the QP path against qp_reference: weights relative to max |w|, the multipliers of the
// equalities by the metric of forward_error (max_k |lambda_k - lambda*_k| / max_k |lambda*_k|), and beta of
// the held weights over max(|lambda_budget*|, |lambda_return*| max |mu|, max |beta*|), the size of a term of
// the stationarity.
void measure_qp_case(double condition, std::size_t n, std::uint64_t seed, double lower)
{
    po::Data data = po::Data::create(n, po::testing::random_vector(n, seed + 500, 0.03, 0.15),
                                     po::testing::random_covariance(n, seed, condition))
                        .value();
    po::Spec spec;
    spec.lower = lower;
    Result minimum = po::testing::frontier_point(spec, data, std::nullopt);
    double largest_mean = po::testing::largest_magnitude(data.mu());
    spec.target_return = 0.5 * (po::testing::portfolio_return(data, minimum.weights) + largest_mean);
    Result result = po::solve(spec, data).value();
    std::optional<QpReference> reference;
    if (result.status == SolveStatus::Optimal)
    {
        reference = qp_reference(spec, data, result);
    }
    if (!reference.has_value())
    {
        std::printf("| %.0e | %zu | %llu | %g | status %d or singular |\n", condition, n,
                    static_cast<unsigned long long>(seed), lower, static_cast<int>(result.status));
        return;
    }
    long double largest = 0.0L;
    long double weights = 0.0L;
    for (std::size_t j = 0; j < n; ++j)
    {
        largest = std::fmax(largest, std::fabs(static_cast<long double>(result.weights[j])));
        weights = std::fmax(weights, std::fabs(result.weights[j] - reference->weights[j]));
    }
    long double multiplier_scale =
        std::fmax(std::fabs(reference->multipliers[0]), std::fabs(reference->multipliers[1]));
    long double budget = std::fabs(result.budget_multiplier - reference->multipliers[0]);
    long double multipliers =
        std::fmax(budget, std::fabs(result.target_return_multiplier - reference->multipliers[1]));
    long double bound_scale =
        std::fmax(std::fabs(reference->multipliers[0]), std::fabs(reference->multipliers[1]) * largest_mean);
    long double bounds = 0.0L;
    int held = 0;
    for (std::size_t j = 0; j < n; ++j)
    {
        if (result.bound_multipliers[j] != 0.0)
        {
            ++held;
            bound_scale = std::fmax(bound_scale, std::fabs(reference->bounds[j]));
            bounds = std::fmax(bounds, std::fabs(result.bound_multipliers[j] - reference->bounds[j]));
        }
    }
    std::printf("| %.0e | %zu | %llu | %g | %d | %.1e | %.1e | %.1e | %.1Le |\n", condition, n,
                static_cast<unsigned long long>(seed), lower, held, static_cast<double>(weights / largest),
                static_cast<double>(multipliers / multiplier_scale),
                static_cast<double>(bounds / bound_scale), reference->last_step);
}

// The multipliers of the QP path (Goldfarb-Idnani, MIN_VARIANCE with a target return) against the KKT point
// of the result's active set solved in long double, by cond(Sigma) up to 1e8, no box and long only: the
// measurement of item 5 of docs/research/2026-11-16-multiplier-accuracy.md, recorded in
// docs/guide/validation.md. Run with --gtest_also_run_disabled_tests.
TEST(RiskParityMultiplierMeasure, DISABLED_MeasureQpMultipliers)
{
    if (!po::testing::kLongDoubleIsExtended)
    {
        GTEST_SKIP() << "long double is double here: the reference has no extra precision";
    }
    std::printf("| cond | n | seed | lower | held | w error | lambda error | beta error | last step |\n");
    for (double condition : {1e4, 1e6, 1e8})
    {
        for (double lower : {-std::numeric_limits<double>::infinity(), 0.0})
        {
            for (std::size_t n : {5u, 8u, 20u})
            {
                for (std::uint64_t seed = 1; seed <= 4; ++seed)
                {
                    measure_qp_case(condition, n, seed, lower);
                }
            }
        }
    }
}

} // namespace
