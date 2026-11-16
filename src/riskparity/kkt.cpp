#include "riskparity/kkt.hpp"

#include <cmath>
#include <cstddef>
#include <limits>

#include "problems/constraints.hpp"
#include "qp/goldfarb_idnani.hpp"

namespace po::riskparity
{
namespace
{

constexpr double kEpsilon = std::numeric_limits<double>::epsilon();

// p = Sigma w, vol = sqrt(w'p) and c_i = w_i p_i / w'p at one point.
struct Moments
{
    std::vector<double> product;
    double volatility;
    std::vector<double> contributions;
};

Moments moments_of(const Data& data, std::span<const double> w)
{
    Moments moments{covariance_times(data, w), 0.0, std::vector<double>(w.size())};
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
        std::vector<double> gradient = moments.product;
        for (double& entry : gradient)
        {
            entry /= moments.volatility;
        }
        double noise = kEpsilon * (term_size(w, moments.product) / moments.volatility + target);
        add_equality(gradient, moments.volatility - target, noise, evaluation);
    }
}

// hess f - lambda_vol hess vol, with hess vol = Sigma / vol - p p' / vol^3.
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
    double cube = volatility * volatility * volatility;
    for (std::size_t j = 0; j < n; ++j)
    {
        for (std::size_t l = 0; l < n; ++l)
        {
            hessian[j * n + l] -= multiplier * (covariance[j * n + l] / volatility - p[j] * p[l] / cube);
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
                          std::size_t j)
{
    double beta = point.bound_multipliers[j];
    if (beta == 0.0)
    {
        return 0.0;
    }
    double distance = distance_to_named_bound(problem.spec, point, j, weight_scale_of(point.weights));
    double multiplier = in_noise_units(beta, evaluation.stationarity_noise[j]);
    if (std::isnan(distance) || std::isnan(multiplier))
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::fmin(distance, multiplier);
}

} // namespace

KktEvaluation evaluate_kkt(const KktProblem& problem, const KktPoint& point)
{
    Moments moments = moments_of(problem.data, point.weights);
    KktEvaluation evaluation;
    evaluation.linearization = linearize(problem.data, problem.budget, point.weights, true);
    add_equalities(problem, point.weights, moments, evaluation);
    std::span<const double> objective_hessian = evaluation.linearization.hessian;
    evaluation.lagrangian_hessian = lagrangian_hessian(problem, point, moments, objective_hessian);
    add_stationarity(problem, point, moments, evaluation);
    return evaluation;
}

KktResidual residual_of(const KktProblem& problem, const KktPoint& point, const KktEvaluation& evaluation)
{
    KktResidual residual{0.0, 0.0, 0.0, 0.0};
    double weight_scale = weight_scale_of(point.weights);
    for (std::size_t j = 0; j < point.weights.size(); ++j)
    {
        double stationarity = in_noise_units(evaluation.stationarity[j], evaluation.stationarity_noise[j]);
        residual.stationarity = worst_of(residual.stationarity, stationarity);
        double violation = bound_violation(problem.spec, point.weights, j, weight_scale);
        double complementarity = complementarity_of(problem, point, evaluation, j);
        residual.bounds = worst_of(residual.bounds, violation);
        residual.complementarity = worst_of(residual.complementarity, complementarity);
    }
    for (std::size_t k = 0; k < evaluation.gaps.size(); ++k)
    {
        double gap = in_noise_units(evaluation.gaps[k], evaluation.gap_noise[k]);
        residual.equalities = worst_of(residual.equalities, gap);
    }
    return residual;
}

double largest(const KktResidual& residual)
{
    double worst = worst_of(residual.stationarity, residual.equalities);
    worst = worst_of(worst, residual.bounds);
    return worst_of(worst, residual.complementarity);
}

bool is_kkt_point(const KktProblem& problem, const KktPoint& point, double factor)
{
    KktEvaluation evaluation = evaluate_kkt(problem, point);
    bool resolves_budget = evaluation.gap_noise.front() <= qp::kRelativeTolerance;
    return resolves_budget && largest(residual_of(problem, point, evaluation)) <= factor;
}

} // namespace po::riskparity
