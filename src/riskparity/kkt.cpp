#include "riskparity/kkt.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>

#include "problems/constraints.hpp"
#include "qp/goldfarb_idnani.hpp"

namespace po::riskparity
{
namespace
{

constexpr double kEpsilon = std::numeric_limits<double>::epsilon();

// p = Sigma w, vol = sqrt(w'p) and c_i = w_i p_i / w'p at one point, w'p compensated.
struct Moments
{
    std::span<const double> product;
    double volatility;
    std::vector<double> contributions;
};

Moments moments_of(std::span<const double> product, std::span<const double> w)
{
    Moments moments{product, 0.0, std::vector<double>(w.size())};
    double variance = compensated_dot(w, moments.product);
    moments.volatility = std::sqrt(variance);
    for (std::size_t i = 0; i < w.size(); ++i)
    {
        moments.contributions[i] = w[i] * moments.product[i] / variance;
    }
    return moments;
}

// sum_j |a_j b_j|: the size of the terms of a'b, which one rounding of each moves by eps times this.
double term_size(std::span<const double> a, std::span<const double> b)
{
    double total = 0.0;
    for (std::size_t j = 0; j < a.size(); ++j)
    {
        total += std::fabs(a[j] * b[j]);
    }
    return total;
}

// A NaN candidate wins, so that a NaN residual is never hidden by a finite one.
double worst_of(double current, double candidate)
{
    if (std::isnan(candidate) || candidate > current)
    {
        return candidate;
    }
    return current;
}

// 0 for a zero residual whatever its noise, so that an exact zero over a zero noise is not NaN.
double in_noise_units(double residual, double noise)
{
    if (residual == 0.0)
    {
        return 0.0;
    }
    return std::fabs(residual) / noise;
}

void add_equality(std::span<const double> row, double gap, double noise, KktEvaluation& evaluation)
{
    evaluation.rows.insert(evaluation.rows.end(), row.begin(), row.end());
    evaluation.gaps.push_back(gap);
    evaluation.gap_noise.push_back(noise);
}

// The rounding of vol(w) = sqrt(w'p) when each w_j is rounded once is sum_j |w_j p_j| eps / vol.
void add_equalities(const KktProblem& problem, std::span<const double> w, const Moments& moments,
                    KktEvaluation& evaluation)
{
    std::vector<double> ones(w.size(), 1.0);
    add_equality(ones, compensated_dot(ones, w) - 1.0, kEpsilon * (term_size(ones, w) + 1.0), evaluation);
    if (problem.spec.target_return.has_value())
    {
        double target = *problem.spec.target_return;
        std::span<const double> mu = problem.data.mu();
        double noise = kEpsilon * (term_size(mu, w) + std::fabs(target));
        add_equality(mu, compensated_dot(mu, w) - target, noise, evaluation);
    }
    if (problem.spec.target_volatility.has_value())
    {
        double target = *problem.spec.target_volatility;
        std::vector<double> gradient(moments.product.begin(), moments.product.end());
        for (double& entry : gradient)
        {
            entry /= moments.volatility;
        }
        double noise = kEpsilon * (term_size(w, moments.product) / moments.volatility + target);
        add_equality(gradient, moments.volatility - target, noise, evaluation);
    }
}

// hess f - lambda_vol hess vol, with hess vol = Sigma / vol - p p' / vol^3. Each p is divided by vol before
// the product: p_j p_l is of the order of vol^4 and leaves the doubles once vol passes 1e77 or falls below
// 1e-77.
std::vector<double> lagrangian_hessian(const KktProblem& problem, const KktPoint& point,
                                       const Moments& moments, std::span<const double> objective_hessian)
{
    std::vector<double> hessian(objective_hessian.begin(), objective_hessian.end());
    if (!problem.spec.target_volatility.has_value())
    {
        return hessian;
    }
    std::size_t n = point.weights.size();
    std::span<const double> covariance = problem.data.covariance();
    std::span<const double> p = moments.product;
    double multiplier = point.equality_multipliers.back();
    double volatility = moments.volatility;
    for (std::size_t j = 0; j < n; ++j)
    {
        for (std::size_t l = 0; l < n; ++l)
        {
            double outer = (p[j] / volatility) * (p[l] / volatility) / volatility;
            hessian[j * n + l] -= multiplier * (covariance[j * n + l] / volatility - outer);
        }
    }
    return hessian;
}

// 2 sum_i |J_ij| (|c_i| + b_i): how far grad f_j = 2 sum_i J_ij (c_i - b_i) moves when c_i and b_i are each
// rounded once.
double gradient_noise(const KktProblem& problem, const Moments& moments, std::span<const double> jacobian,
                      std::size_t j)
{
    std::size_t n = moments.contributions.size();
    double total = 0.0;
    for (std::size_t i = 0; i < n; ++i)
    {
        double gap_size = std::fabs(moments.contributions[i]) + std::fabs(problem.budget[i]);
        total += 2.0 * std::fabs(jacobian[i * n + j]) * gap_size;
    }
    return total;
}

// The residual grad f - sum_k lambda_k grad h_k - beta and its noise: the terms of the residual, each
// rounded once, and (|H| |w|)_j, the change of the residual when each weight is rounded once, H the Hessian
// of the Lagrangian.
void add_stationarity(const KktProblem& problem, const KktPoint& point, const Moments& moments,
                      KktEvaluation& evaluation)
{
    std::size_t n = point.weights.size();
    std::span<const double> hessian = evaluation.lagrangian_hessian;
    evaluation.stationarity = evaluation.linearization.gradient;
    evaluation.stationarity_noise.assign(n, 0.0);
    for (std::size_t j = 0; j < n; ++j)
    {
        double size = term_size(hessian.subspan(j * n, n), point.weights);
        size += gradient_noise(problem, moments, evaluation.linearization.jacobian, j);
        for (std::size_t k = 0; k < evaluation.gaps.size(); ++k)
        {
            double term = point.equality_multipliers[k] * evaluation.rows[k * n + j];
            evaluation.stationarity[j] -= term;
            size += std::fabs(term);
        }
        evaluation.stationarity[j] -= point.bound_multipliers[j];
        size += std::fabs(point.bound_multipliers[j]);
        evaluation.stationarity_noise[j] = kEpsilon * size;
    }
}

// One rounding of a weight at its bound: eps max(|bound|, |w|max).
double bound_noise(double bound, double weight_scale)
{
    return kEpsilon * std::fmax(std::fabs(bound), weight_scale);
}

double weight_scale_of(std::span<const double> weights)
{
    double scale = 0.0;
    for (double weight : weights)
    {
        scale = worst_of(scale, std::fabs(weight));
    }
    return scale;
}

double bound_violation(const Spec& spec, std::span<const double> weights, std::size_t j, double weight_scale)
{
    double lower = problems::lower_of(spec, j);
    double upper = problems::upper_of(spec, j);
    if (std::isnan(weights[j]))
    {
        return weights[j];
    }
    if (weights[j] < lower)
    {
        return in_noise_units(lower - weights[j], bound_noise(lower, weight_scale));
    }
    if (weights[j] > upper)
    {
        return in_noise_units(weights[j] - upper, bound_noise(upper, weight_scale));
    }
    return 0.0;
}

// The distance of w_j to the bound that beta_j names, the lower bound when beta_j is positive and the upper
// one when it is negative; infinite when that bound is.
double distance_to_named_bound(const Spec& spec, const KktPoint& point, std::size_t j, double weight_scale)
{
    double bound = problems::upper_of(spec, j);
    if (point.bound_multipliers[j] > 0.0)
    {
        bound = problems::lower_of(spec, j);
    }
    if (std::isinf(bound))
    {
        return std::numeric_limits<double>::infinity();
    }
    return in_noise_units(point.weights[j] - bound, bound_noise(bound, weight_scale));
}

// Either w_j is at the bound that beta_j names or beta_j is noise.
double complementarity_of(const KktProblem& problem, const KktPoint& point, const KktEvaluation& evaluation,
                          std::size_t j, double weight_scale)
{
    double beta = point.bound_multipliers[j];
    if (beta == 0.0)
    {
        return 0.0;
    }
    double distance = distance_to_named_bound(problem.spec, point, j, weight_scale);
    double multiplier = in_noise_units(beta, evaluation.stationarity_noise[j]);
    if (std::isnan(distance) || std::isnan(multiplier))
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::fmin(distance, multiplier);
}

double equality_residual(const KktEvaluation& evaluation)
{
    double residual = 0.0;
    for (std::size_t k = 0; k < evaluation.gaps.size(); ++k)
    {
        residual = worst_of(residual, in_noise_units(evaluation.gaps[k], evaluation.gap_noise[k]));
    }
    return residual;
}

double bound_residual(const Spec& spec, std::span<const double> weights, double weight_scale)
{
    double residual = 0.0;
    for (std::size_t j = 0; j < weights.size(); ++j)
    {
        residual = worst_of(residual, bound_violation(spec, weights, j, weight_scale));
    }
    return residual;
}

// The rows, the gaps and their noise, which the weights alone decide. Of the linearization only Sigma w is
// set, for add_linearization.
KktEvaluation equalities_at(const KktProblem& problem, std::span<const double> w)
{
    KktEvaluation evaluation;
    evaluation.linearization.product = covariance_times(problem.data, w);
    add_equalities(problem, w, moments_of(evaluation.linearization.product, w), evaluation);
    return evaluation;
}

// The rest of evaluate_kkt after equalities_at: f with its Hessian, and what the multipliers decide.
void add_linearization(const KktProblem& problem, const KktPoint& point, KktEvaluation& evaluation)
{
    std::vector<double> product = std::move(evaluation.linearization.product);
    evaluation.linearization = linearize(problem.data, problem.budget, point.weights, std::move(product), true);
    update_multipliers(problem, point, evaluation);
}

} // namespace

KktEvaluation evaluate_kkt(const KktProblem& problem, const KktPoint& point)
{
    KktEvaluation evaluation = equalities_at(problem, point.weights);
    add_linearization(problem, point, evaluation);
    return evaluation;
}

void update_multipliers(const KktProblem& problem, const KktPoint& point, KktEvaluation& evaluation)
{
    Moments moments = moments_of(evaluation.linearization.product, point.weights);
    std::span<const double> objective_hessian = evaluation.linearization.hessian;
    evaluation.lagrangian_hessian = lagrangian_hessian(problem, point, moments, objective_hessian);
    add_stationarity(problem, point, moments, evaluation);
}

KktResidual residual_of(const KktProblem& problem, const KktPoint& point, const KktEvaluation& evaluation)
{
    double weight_scale = weight_scale_of(point.weights);
    double bounds = bound_residual(problem.spec, point.weights, weight_scale);
    KktResidual residual{0.0, equality_residual(evaluation), bounds, 0.0};
    for (std::size_t j = 0; j < point.weights.size(); ++j)
    {
        double stationarity = in_noise_units(evaluation.stationarity[j], evaluation.stationarity_noise[j]);
        residual.stationarity = worst_of(residual.stationarity, stationarity);
        double complementarity = complementarity_of(problem, point, evaluation, j, weight_scale);
        residual.complementarity = worst_of(residual.complementarity, complementarity);
    }
    return residual;
}

double largest(const KktResidual& residual)
{
    double worst = worst_of(residual.stationarity, residual.equalities);
    worst = worst_of(worst, residual.bounds);
    return worst_of(worst, residual.complementarity);
}

bool passes_kkt_gate(const KktEvaluation& evaluation, double residual, double factor)
{
    bool resolves_budget = evaluation.gap_noise.front() <= qp::kRelativeTolerance;
    return resolves_budget && residual <= factor;
}

// largest is NaN or at least each of its parts, so a part that the weights alone decide and that is not within
// the factor rejects the point before the Hessian is made.
bool is_kkt_point(const KktProblem& problem, const KktPoint& point, double factor)
{
    KktEvaluation evaluation = equalities_at(problem, point.weights);
    double bounds = bound_residual(problem.spec, point.weights, weight_scale_of(point.weights));
    double weights_part = worst_of(equality_residual(evaluation), bounds);
    if (!passes_kkt_gate(evaluation, weights_part, factor))
    {
        return false;
    }
    add_linearization(problem, point, evaluation);
    return passes_kkt_gate(evaluation, largest(residual_of(problem, point, evaluation)), factor);
}

} // namespace po::riskparity
