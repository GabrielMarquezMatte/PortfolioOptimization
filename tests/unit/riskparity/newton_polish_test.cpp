#include "riskparity/newton_polish.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/data.hpp"
#include "core/objective.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "fuzz/frontier_reference.hpp"
#include "fuzz/risk_parity/risk_parity_reference.hpp"
#include "fuzz/risk_parity/risk_parity_refinement.hpp"
#include "random_covariance.hpp"
#include "riskparity/kkt.hpp"
#include "riskparity/polyhedral_descent.hpp"
#include "solve/solve.hpp"

namespace
{

using po::Data;
using po::Result;
using po::Spec;
using po::riskparity::DescentOptions;
using po::riskparity::DescentProblem;
using po::riskparity::DescentResult;
using po::riskparity::DescentStop;
using po::riskparity::KktEvaluation;
using po::riskparity::KktPoint;
using po::riskparity::KktProblem;
using po::riskparity::PolishResult;
using po::testing::random_covariance;
using po::testing::random_vector;

constexpr double kBand = 1e-10;
constexpr double kStartTolerance = 1e-6;

Data sample_data(std::size_t n, std::uint64_t seed)
{
    return Data::create(n, random_vector(n, seed + 500, 0.03, 0.15), random_covariance(n, seed, 1e2)).value();
}

// Box [0, 1.2 / n] and an uneven budget, which holds weights at the cap.
Spec capped_spec(std::size_t n)
{
    Spec spec;
    spec.objective = po::Objective::RiskParity;
    spec.lower = 0.0;
    spec.upper = 1.2 / static_cast<double>(n);
    spec.risk_budget = {0.5, 0.2, 0.15, 0.1, 0.05};
    return spec;
}

// A descent stopped at a relative KKT residual of 1e-6: a point near a KKT point, with the active set of its
// last QP in the bound multipliers.
KktPoint near_kkt_point(const Spec& spec, const Data& data)
{
    std::vector<double> start = po::testing::frontier_point(spec, data, spec.target_return).weights;
    DescentProblem problem{spec, data, spec.risk_budget, std::nullopt};
    DescentResult descent =
        po::riskparity::descend(problem, start, DescentOptions{0.05, kStartTolerance, 1000});
    EXPECT_EQ(descent.stop, DescentStop::Converged);
    EXPECT_GT(descent.stationarity, 100.0 * kBand);
    return KktPoint{descent.weights, descent.equality_multipliers, descent.bound_multipliers};
}

Result as_result(const Spec& spec, const KktPoint& point)
{
    Result result;
    result.status = po::SolveStatus::Optimal;
    result.weights = point.weights;
    result.bound_multipliers = point.bound_multipliers;
    result.budget_multiplier = point.equality_multipliers.at(0);
    if (spec.target_return.has_value())
    {
        result.target_return_multiplier = point.equality_multipliers.at(1);
    }
    if (spec.target_volatility.has_value())
    {
        result.volatility_multiplier = point.equality_multipliers.back();
    }
    return result;
}

double merit_of(const KktProblem& problem, const KktPoint& point)
{
    return largest(residual_of(problem, point, evaluate_kkt(problem, point)));
}

// Every number of two evaluations, bit for bit.
void expect_same_evaluation(const KktEvaluation& first, const KktEvaluation& second)
{
    EXPECT_EQ(first.linearization.value, second.linearization.value);
    EXPECT_EQ(first.linearization.product, second.linearization.product);
    EXPECT_EQ(first.linearization.variance, second.linearization.variance);
    EXPECT_EQ(first.linearization.gradient, second.linearization.gradient);
    EXPECT_EQ(first.linearization.jacobian, second.linearization.jacobian);
    EXPECT_EQ(first.linearization.hessian, second.linearization.hessian);
    EXPECT_EQ(first.lagrangian_hessian, second.lagrangian_hessian);
    EXPECT_EQ(first.rows, second.rows);
    EXPECT_EQ(first.gaps, second.gaps);
    EXPECT_EQ(first.gap_noise, second.gap_noise);
    EXPECT_EQ(first.stationarity, second.stationarity);
    EXPECT_EQ(first.stationarity_noise, second.stationarity_noise);
}

// The multipliers of two points agree to the rounding of the least-squares solve.
void expect_same_multipliers(const KktPoint& first, const KktPoint& second)
{
    ASSERT_EQ(first.equality_multipliers.size(), second.equality_multipliers.size());
    for (std::size_t k = 0; k < first.equality_multipliers.size(); ++k)
    {
        double expected = second.equality_multipliers[k];
        EXPECT_NEAR(first.equality_multipliers[k], expected, 1e-13 * std::fabs(expected)) << "row " << k;
    }
    for (std::size_t j = 0; j < first.bound_multipliers.size(); ++j)
    {
        double scale = std::fabs(second.bound_multipliers[j]);
        EXPECT_NEAR(first.bound_multipliers[j], second.bound_multipliers[j], 1e-13 * scale) << "weight " << j;
    }
}

// From a KKT point at 1e-6, one or two Newton steps reach the band of cond(Sigma) and the rounding floor of
// the residual; a second polish takes no step from there and leaves the point where it is.
TEST(NewtonPolish, ReachesRoundingFloor)
{
    std::size_t n = 5;
    Data data = sample_data(n, 7);
    for (std::optional<double> target : {std::optional<double>(), std::optional<double>(0.09)})
    {
        SCOPED_TRACE(target.has_value());
        Spec spec = capped_spec(n);
        spec.target_return = target;
        KktProblem problem{spec, data, spec.risk_budget};
        KktPoint start = near_kkt_point(spec, data);
        double start_merit = merit_of(problem, start);
        PolishResult polished = po::riskparity::newton_polish(problem, start);
        ASSERT_LT(polished.residual, start_merit);
        EXPECT_GE(polished.steps, 1);
        EXPECT_LE(polished.steps, 2);
        EXPECT_EQ(polished.residual, merit_of(problem, polished.least_squares_point));
        expect_same_evaluation(polished.evaluation, evaluate_kkt(problem, polished.least_squares_point));
        EXPECT_LT(polished.residual, 1e-6 * start_merit);
        EXPECT_LE(polished.residual, 1.0);
        po::testing::expect_risk_parity_kkt(spec, data, as_result(spec, polished.least_squares_point), kBand);
        bool capped = false;
        for (std::size_t i = 0; i < n; ++i)
        {
            capped = capped || polished.least_squares_point.bound_multipliers[i] < 0.0;
        }
        EXPECT_TRUE(capped);
        PolishResult again = po::riskparity::newton_polish(problem, polished.least_squares_point);
        EXPECT_EQ(again.steps, 0);
        EXPECT_EQ(again.least_squares_point.weights, polished.least_squares_point.weights);
        expect_same_multipliers(again.least_squares_point, polished.least_squares_point);
        EXPECT_LE(again.residual, 1.0);
        expect_same_evaluation(again.evaluation, evaluate_kkt(problem, again.least_squares_point));
    }
}

// The multipliers of the equalities at the start are not used: they are settled by least squares at its
// weights before any step. Whatever they were, the polish ends at the same KKT point.
TEST(NewtonPolish, DoesNotDependOnTheEqualityMultipliersOfTheStart)
{
    std::size_t n = 5;
    Data data = sample_data(n, 7);
    Spec spec = capped_spec(n);
    KktProblem problem{spec, data, spec.risk_budget};
    PolishResult reference = po::riskparity::newton_polish(problem, near_kkt_point(spec, data));
    ASSERT_LE(reference.residual, 1.0);
    for (double factor : {-50.0, 0.0, 3.0, 1e6})
    {
        SCOPED_TRACE(factor);
        KktPoint start = near_kkt_point(spec, data);
        for (double& multiplier : start.equality_multipliers)
        {
            multiplier *= factor;
        }
        PolishResult polished = po::riskparity::newton_polish(problem, start);
        EXPECT_LE(polished.residual, 1.0);
        EXPECT_EQ(polished.steps, reference.steps);
        expect_same_multipliers(polished.least_squares_point, reference.least_squares_point);
        for (std::size_t i = 0; i < n; ++i)
        {
            double expected = reference.least_squares_point.weights[i];
            EXPECT_NEAR(polished.least_squares_point.weights[i], expected, 1e-14) << "weight " << i;
        }
    }
}

// The largest and the mean stationarity of the free weights (those whose bound multiplier is 0), the mean in
// units of the smallest noise among them: the component of the residual along the budget row, which the
// least-squares multipliers remove.
double budget_component(const KktProblem& problem, const KktPoint& point)
{
    po::riskparity::KktEvaluation evaluation = evaluate_kkt(problem, point);
    double total = 0.0;
    double count = 0.0;
    double smallest_noise = std::numeric_limits<double>::infinity();
    for (std::size_t j = 0; j < point.weights.size(); ++j)
    {
        if (point.bound_multipliers[j] != 0.0)
        {
            continue;
        }
        total += evaluation.stationarity[j];
        count += 1.0;
        smallest_noise = std::fmin(smallest_noise, evaluation.stationarity_noise[j]);
    }
    return std::fabs(total / count) / smallest_noise;
}

double smallest_free_noise(const KktProblem& problem, const KktPoint& point)
{
    po::riskparity::KktEvaluation evaluation = evaluate_kkt(problem, point);
    double smallest = std::numeric_limits<double>::infinity();
    for (std::size_t j = 0; j < point.weights.size(); ++j)
    {
        if (point.bound_multipliers[j] == 0.0)
        {
            smallest = std::fmin(smallest, evaluation.stationarity_noise[j]);
        }
    }
    return smallest;
}

// A start within its noise takes no Newton step and still comes back with the least-squares multipliers of
// its weights: a budget multiplier off by half the noise of the stationarity, which leaves that much of the
// residual along the budget row, is put back to what the least-squares fit at those weights leaves.
TEST(NewtonPolish, SettlesTheMultipliersWithoutAStep)
{
    std::size_t n = 5;
    Data data = sample_data(n, 7);
    Spec spec = capped_spec(n);
    KktProblem problem{spec, data, spec.risk_budget};
    PolishResult polished = po::riskparity::newton_polish(problem, near_kkt_point(spec, data));
    ASSERT_LE(polished.residual, 0.4);
    KktPoint start = polished.least_squares_point;
    double fitted = budget_component(problem, start);
    ASSERT_LT(fitted, 0.25);
    start.equality_multipliers[0] += 0.5 * smallest_free_noise(problem, start);
    ASSERT_LE(merit_of(problem, start), 1.0);
    ASSERT_GE(budget_component(problem, start), 0.25);
    PolishResult settled = po::riskparity::newton_polish(problem, start);
    EXPECT_EQ(settled.steps, 0);
    EXPECT_EQ(settled.least_squares_point.weights, start.weights);
    EXPECT_LE(budget_component(problem, settled.least_squares_point), fitted * (1.0 + 1e-6));
    EXPECT_LE(settled.residual, 1.0);
    EXPECT_EQ(settled.residual, merit_of(problem, settled.least_squares_point));
}

// The bound multipliers of the start name the active set. Here they name a lower bound that the solution does
// not touch: weight k is held 1e-4 below its optimal value, where f decreases into the box and the multiplier
// of the bound comes out negative. The polish releases the bound and returns to the solution.
TEST(NewtonPolish, ReleasesABoundWhoseMultiplierHasTheWrongSign)
{
    std::size_t n = 5;
    Data data = sample_data(n, 7);
    Spec spec = capped_spec(n);
    KktProblem problem{spec, data, spec.risk_budget};
    PolishResult solution = po::riskparity::newton_polish(problem, near_kkt_point(spec, data));
    ASSERT_LE(solution.residual, 1.0);
    std::size_t k = 0;
    while (solution.least_squares_point.bound_multipliers[k] != 0.0)
    {
        ++k;
    }
    Spec shifted = spec;
    shifted.weight_lower.assign(n, spec.lower);
    shifted.weight_upper.assign(n, spec.upper);
    shifted.weight_lower[k] = solution.least_squares_point.weights[k] - 1e-4;
    KktProblem shifted_problem{shifted, data, shifted.risk_budget};
    KktPoint start = solution.least_squares_point;
    start.bound_multipliers[k] = 1.0;
    ASSERT_GT(merit_of(shifted_problem, start), 1e3);
    PolishResult polished = po::riskparity::newton_polish(shifted_problem, start);
    EXPECT_LT(polished.residual, merit_of(shifted_problem, start));
    EXPECT_GE(polished.steps, 1);
    EXPECT_EQ(polished.least_squares_point.bound_multipliers[k], 0.0);
    EXPECT_GT(polished.least_squares_point.weights[k], shifted.weight_lower[k]);
    EXPECT_LE(polished.residual, 1.0);
    for (std::size_t i = 0; i < n; ++i)
    {
        EXPECT_NEAR(polished.least_squares_point.weights[i], solution.least_squares_point.weights[i], 1e-12)
            << "weight " << i;
    }
}

using Shape = po::testing::RiskParityShape;

KktPoint point_of(const Spec& spec, const Result& result)
{
    KktPoint point{result.weights, {result.budget_multiplier}, result.bound_multipliers};
    if (spec.target_return.has_value())
    {
        point.equality_multipliers.push_back(result.target_return_multiplier);
    }
    if (spec.target_volatility.has_value())
    {
        point.equality_multipliers.push_back(result.volatility_multiplier);
    }
    return point;
}

// The cases of the measurement of the Newton multipliers in docs/guide/validation.md at one cond(Sigma), no
// box and equal budget: three shapes at s = 1.2 sigma_min and the 8 of
// RiskParityFuzz.TargetVolatilityMultipliersByCondition.
std::vector<po::testing::GeneratedRiskParity> multiplier_cases(double condition)
{
    std::vector<po::testing::GeneratedRiskParity> cases;
    for (auto [n, seed] : {std::pair<std::size_t, std::uint64_t>{5, 7}, {5, 11}, {8, 3}})
    {
        po::testing::GeneratedRiskParity generated = po::testing::generate(
            Shape{n, condition, 1.0, Shape::Budget::Equal, Shape::Box::None, false, true, seed});
        Result minimum = po::testing::frontier_point(generated.spec, generated.data, std::nullopt);
        double smallest = po::testing::portfolio_volatility(generated.data, minimum.weights);
        generated.spec.target_volatility = 1.2 * smallest;
        cases.push_back(std::move(generated));
    }
    for (std::size_t n : {5u, 8u})
    {
        for (std::uint64_t seed = 402; seed <= 405; ++seed)
        {
            cases.push_back(po::testing::generate(
                Shape{n, condition, 1.0, Shape::Budget::Equal, Shape::Box::None, false, true, seed}));
        }
    }
    return cases;
}

// forward_error(...).multipliers of the least-squares and of the corrected multipliers of a polish.
struct MultiplierErrors
{
    double least_squares;
    double corrected;
};

// The corrected multipliers are those po::solve returned. A polish of that result takes no step and settles
// the same weights again: its least-squares multipliers are those the KKT gate checked (to their rounding,
// which depends on the multipliers the least-squares fit starts from), and it accepts the correction there.
MultiplierErrors multiplier_errors(const po::testing::GeneratedRiskParity& generated)
{
    const Spec& spec = generated.spec;
    Result solved = po::solve(spec, generated.data).value();
    EXPECT_EQ(solved.status, po::SolveStatus::Optimal);
    KktProblem problem{spec, generated.data, spec.risk_budget};
    PolishResult polished = po::riskparity::newton_polish(problem, point_of(spec, solved));
    EXPECT_EQ(polished.steps, 0);
    po::riskparity::MultiplierCorrection correction =
        po::riskparity::corrected_multipliers(problem, polished.least_squares_point, polished.evaluation);
    EXPECT_TRUE(correction.accepted);
    EXPECT_EQ(correction.reported_point.weights, solved.weights);
    Result least_squares = as_result(spec, polished.least_squares_point);
    return MultiplierErrors{po::testing::forward_error(spec, generated.data, least_squares).multipliers,
                            po::testing::forward_error(spec, generated.data, solved).multipliers};
}

struct MultiplierCeiling
{
    double condition;
    double corrected;
};

// The multipliers of the target-volatility path against those of the KKT point refined in long double, with
// the weights unchanged. The least-squares estimate at the rounded weights loses the correction H dw that the
// Newton solve at the end of the polish keeps. The ceilings of the corrected multipliers are the largest
// measured on the cases of the instances before the portable generators times 10, not derived
// (docs/guide/validation.md): 5.2e-11 at 1e6 and 7.1e-10 at 1e8. On the present cases the largest are 2.4e-12
// and 9.3e-11, 217 and 77 times below the ceilings. The least-squares ones, 4.2e-4 and 2.2e-2 (8.2e-7 and 3.7e-2
// before), stay more than 100 times above them, so that the test tells the two apart.
TEST(NewtonPolish, KeepsTheNewtonMultiplierCorrection)
{
    if (!po::testing::kLongDoubleIsExtended)
    {
        GTEST_SKIP() << "long double is double here: the refinement has no extra precision";
    }
    for (MultiplierCeiling ceiling : {MultiplierCeiling{1e6, 5.2e-10}, MultiplierCeiling{1e8, 7.2e-9}})
    {
        SCOPED_TRACE(ceiling.condition);
        double worst_least_squares = 0.0;
        double worst_corrected = 0.0;
        for (const po::testing::GeneratedRiskParity& generated : multiplier_cases(ceiling.condition))
        {
            MultiplierErrors errors = multiplier_errors(generated);
            EXPECT_LE(errors.corrected, ceiling.corrected);
            worst_least_squares = std::fmax(worst_least_squares, errors.least_squares);
            worst_corrected = std::fmax(worst_corrected, errors.corrected);
        }
        std::string exponent = std::to_string(static_cast<int>(std::log10(ceiling.condition)));
        RecordProperty("least_squares_" + exponent, ::testing::PrintToString(worst_least_squares));
        RecordProperty("corrected_" + exponent, ::testing::PrintToString(worst_corrected));
        EXPECT_GT(worst_least_squares, 100.0 * ceiling.corrected);
    }
}

// A point far from a KKT point of its problem: the solution at s = 1.1 sigma_min, multipliers included, for
// the problem at s = sigma_min (1 + 1e-4), where the volatility multiplier is much larger and the weights
// must move toward the minimum-variance portfolio. A second solve from the corrected multipliers changes them
// by more than half of the first correction, which is refused: the multipliers stay those of the point.
TEST(NewtonPolish, RejectsACorrectionThatDoesNotContract)
{
    std::size_t n = 5;
    Data data = sample_data(n, 7);
    Spec spec;
    spec.objective = po::Objective::RiskParity;
    spec.risk_budget.assign(n, 1.0 / static_cast<double>(n));
    Result minimum = po::testing::frontier_point(spec, data, std::nullopt);
    double smallest = po::testing::portfolio_volatility(data, minimum.weights);
    spec.target_volatility = 1.1 * smallest;
    Result far = po::solve(spec, data).value();
    ASSERT_EQ(far.status, po::SolveStatus::Optimal);
    spec.target_volatility = smallest * (1.0 + 1e-4);
    KktProblem problem{spec, data, spec.risk_budget};
    KktPoint point = point_of(spec, far);
    po::riskparity::MultiplierCorrection correction =
        po::riskparity::corrected_multipliers(problem, point, evaluate_kkt(problem, point));
    RecordProperty("contraction", ::testing::PrintToString(correction.second_change / correction.change));
    EXPECT_GT(correction.change, 0.0);
    EXPECT_GT(correction.second_change, 0.5 * correction.change);
    EXPECT_FALSE(correction.accepted);
    EXPECT_EQ(correction.reported_point.weights, point.weights);
    EXPECT_EQ(correction.reported_point.equality_multipliers, point.equality_multipliers);
    EXPECT_EQ(correction.reported_point.bound_multipliers, point.bound_multipliers);
}

// A held weight named at the other side of its bound: in a long-only solution with a target volatility
// (cond(Sigma) = 1e8, where the correction moves the multipliers by about 4e-8 and its second solve shrinks),
// the weight held at 0 (beta > 0) is given beta < 0, which names an upper bound. The correction gives it the
// beta of the lower bound again, of the wrong sign for the side named, and is refused: the multipliers stay
// those of the point. The same point (the least-squares point of a polish of the solution) with its own beta
// keeps the correction.
TEST(NewtonPolish, RejectsACorrectionThatChangesTheSideOfABound)
{
    po::testing::GeneratedRiskParity generated = po::testing::generate(
        Shape{5, 1e8, 1.0, Shape::Budget::Uneven, Shape::Box::LongOnly, false, true, 404});
    const Spec& spec = generated.spec;
    Result solved = po::solve(spec, generated.data).value();
    ASSERT_EQ(solved.status, po::SolveStatus::Optimal);
    KktProblem problem{spec, generated.data, spec.risk_budget};
    PolishResult polished = po::riskparity::newton_polish(problem, point_of(spec, solved));
    ASSERT_EQ(polished.steps, 0);
    ASSERT_TRUE(
        po::riskparity::corrected_multipliers(problem, polished.least_squares_point, polished.evaluation).accepted);
    KktPoint point = polished.least_squares_point;
    std::size_t k = 0;
    while (k < point.weights.size() && point.bound_multipliers[k] <= 0.0)
    {
        ++k;
    }
    ASSERT_LT(k, point.weights.size());
    point.bound_multipliers[k] = -point.bound_multipliers[k];
    po::riskparity::MultiplierCorrection correction =
        po::riskparity::corrected_multipliers(problem, point, evaluate_kkt(problem, point));
    RecordProperty("contraction", ::testing::PrintToString(correction.second_change / correction.change));
    EXPECT_LE(correction.second_change, 0.5 * correction.change);
    EXPECT_FALSE(correction.accepted);
    EXPECT_EQ(correction.reported_point.weights, point.weights);
    EXPECT_EQ(correction.reported_point.equality_multipliers, point.equality_multipliers);
    EXPECT_EQ(correction.reported_point.bound_multipliers, point.bound_multipliers);
}

// The KKT residual in units of its noise separates a KKT point from a point an error of 1e-9 of a weight
// away from it.
TEST(Kkt, ResidualSeesAnErrorOfTheWeights)
{
    std::size_t n = 5;
    Data data = sample_data(n, 7);
    Spec spec = capped_spec(n);
    spec.target_return = 0.09;
    KktProblem problem{spec, data, spec.risk_budget};
    KktPoint start = near_kkt_point(spec, data);
    PolishResult polished = po::riskparity::newton_polish(problem, start);
    ASSERT_LT(polished.residual, merit_of(problem, start));
    KktPoint moved = polished.least_squares_point;
    std::size_t first_free = n;
    std::size_t second_free = n;
    for (std::size_t i = 0; i < n; ++i)
    {
        if (moved.bound_multipliers[i] != 0.0)
        {
            continue;
        }
        second_free = first_free;
        first_free = i;
    }
    ASSERT_LT(second_free, n);
    moved.weights[first_free] += 1e-9;
    moved.weights[second_free] -= 1e-9;
    po::riskparity::KktResidual residual = residual_of(problem, moved, evaluate_kkt(problem, moved));
    EXPECT_GT(residual.stationarity, 1e3);
    EXPECT_GT(residual.equalities, 1e3);
    EXPECT_EQ(residual.bounds, 0.0);
    KktPoint outside = polished.least_squares_point;
    outside.weights[first_free] = spec.upper + 1e-9;
    EXPECT_GT(residual_of(problem, outside, evaluate_kkt(problem, outside)).bounds, 1e3);
    KktPoint wrong_sign = polished.least_squares_point;
    for (double& beta : wrong_sign.bound_multipliers)
    {
        beta = -beta;
    }
    EXPECT_GT(residual_of(problem, wrong_sign, evaluate_kkt(problem, wrong_sign)).complementarity, 1e3);
    KktPoint not_a_number = polished.least_squares_point;
    not_a_number.weights[0] = std::nan("");
    EXPECT_TRUE(std::isnan(largest(residual_of(problem, not_a_number, evaluate_kkt(problem, not_a_number)))));
}

// is_kkt_point rejects a point by the parts of its residual that the weights alone decide before it makes the
// Hessian; it decides as the full residual does, at a KKT point and at points that fail each part.
TEST(Kkt, IsKktPointDecidesAsTheFullResidual)
{
    std::size_t n = 5;
    Data data = sample_data(n, 7);
    Spec spec = capped_spec(n);
    spec.target_return = 0.09;
    KktProblem problem{spec, data, spec.risk_budget};
    KktPoint solution = po::riskparity::newton_polish(problem, near_kkt_point(spec, data)).least_squares_point;
    std::vector<KktPoint> points(5, solution);
    points[1].weights[0] += 1e-9;
    points[2].weights[0] = spec.upper + 1e-9;
    points[3].equality_multipliers[0] *= 1.0 + 1e-6;
    points[4].weights[0] = std::nan("");
    for (std::size_t p = 0; p < points.size(); ++p)
    {
        SCOPED_TRACE(p);
        KktEvaluation evaluation = evaluate_kkt(problem, points[p]);
        bool expected =
            po::riskparity::passes_kkt_gate(evaluation, largest(residual_of(problem, points[p], evaluation)));
        EXPECT_EQ(is_kkt_point(problem, points[p]), expected);
        EXPECT_EQ(expected, p == 0);
    }
}

// The evaluation of a point whose multipliers change at the same weights follows them with the bits of a new
// evaluation, the volatility multiplier included, which enters the Hessian of the Lagrangian.
TEST(Kkt, UpdateMultipliersGivesTheBitsOfANewEvaluation)
{
    std::size_t n = 5;
    Data data = sample_data(n, 7);
    Spec spec = capped_spec(n);
    spec.target_volatility = 0.2;
    KktProblem problem{spec, data, spec.risk_budget};
    std::vector<double> weights(n, 1.0 / static_cast<double>(n));
    KktPoint point{weights, {0.3, -0.7}, {0.0, 0.1, 0.0, -0.2, 0.0}};
    KktEvaluation evaluation = evaluate_kkt(problem, point);
    point.equality_multipliers = {-1.1, 2.5};
    point.bound_multipliers = {0.4, 0.0, 0.0, -0.2, 0.3};
    po::riskparity::update_multipliers(problem, point, evaluation);
    expect_same_evaluation(evaluation, evaluate_kkt(problem, point));
}

// Two weights pinned by their bounds at (1 + t, -t), with beta = grad f, so that the stationarity and the
// bounds hold exactly. At t = 1 the point is a KKT point. At t = 1e16 the first weight rounds to t and the
// weights sum to 0, not to 1: a gap of 1, which is within the rounding of weights of that size, so the
// residual alone cannot reject the point. The budget row that is no longer resolved does.
TEST(Kkt, WeightsThatDoNotResolveTheBudgetAreNotAKktPoint)
{
    Data data = Data::create(2, {}, std::vector<double>{0.04, 0.006, 0.006, 0.09}).value();
    for (double size : {1.0, 1e16})
    {
        SCOPED_TRACE(size);
        Spec spec;
        spec.objective = po::Objective::RiskParity;
        spec.risk_budget = {0.5, 0.5};
        spec.weight_lower = {1.0 + size, -size};
        spec.weight_upper = spec.weight_lower;
        KktProblem problem{spec, data, spec.risk_budget};
        KktPoint point{spec.weight_lower, {0.0}, {0.0, 0.0}};
        point.bound_multipliers = evaluate_kkt(problem, point).linearization.gradient;
        double residual = largest(residual_of(problem, point, evaluate_kkt(problem, point)));
        EXPECT_LE(residual, 1.0);
        EXPECT_EQ(is_kkt_point(problem, point), size == 1.0);
    }
    EXPECT_EQ((1.0 + 1e16) + (-1e16), 0.0);
}

} // namespace
