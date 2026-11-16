#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "fuzz/risk_parity/risk_parity_reference.hpp"

namespace po::testing
{

// The refinement needs a long double wider than double (kLongDoubleIsExtended): where it is not (MSVC), the
// refined point cannot be held beyond double and is no better than the point it refines. Its callers skip.

// Gaussian elimination with partial pivoting of a size x size system (row-major); false when a pivot is 0.
inline bool solve_wide(std::vector<long double>& matrix, std::vector<long double>& rhs)
{
    std::size_t size = rhs.size();
    for (std::size_t k = 0; k < size; ++k)
    {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < size; ++i)
        {
            if (std::fabs(matrix[i * size + k]) > std::fabs(matrix[pivot * size + k]))
            {
                pivot = i;
            }
        }
        if (matrix[pivot * size + k] == 0.0L)
        {
            return false;
        }
        for (std::size_t j = 0; j < size; ++j)
        {
            std::swap(matrix[k * size + j], matrix[pivot * size + j]);
        }
        std::swap(rhs[k], rhs[pivot]);
        for (std::size_t i = k + 1; i < size; ++i)
        {
            long double factor = matrix[i * size + k] / matrix[k * size + k];
            for (std::size_t j = k; j < size; ++j)
            {
                matrix[i * size + j] -= factor * matrix[k * size + j];
            }
            rhs[i] -= factor * rhs[k];
        }
    }
    for (std::size_t k = size; k-- > 0;)
    {
        for (std::size_t j = k + 1; j < size; ++j)
        {
            rhs[k] -= matrix[k * size + j] * rhs[j];
        }
        rhs[k] /= matrix[k * size + k];
    }
    return true;
}

// The KKT system of RISK_PARITY on the active set of a result, in long double: the weights whose bound
// multiplier is not 0 stay where the result has them, the others and the multipliers of the equalities
// (budget, target return, volatility, those of spec) are the unknowns.
struct WideKkt
{
    std::vector<long double> weights;
    std::vector<long double> multipliers;
    std::vector<std::size_t> free;
};

struct WideEqualities
{
    // One gradient per equality, n entries each.
    std::vector<std::vector<long double>> rows;
    std::vector<long double> gaps;
};

inline WideEqualities wide_equalities(const po::Spec& spec, const po::Data& data,
                                      const std::vector<long double>& w)
{
    std::size_t n = w.size();
    WideEqualities equalities;
    long double total = -1.0L;
    for (long double weight : w)
    {
        total += weight;
    }
    equalities.rows.emplace_back(n, 1.0L);
    equalities.gaps.push_back(total);
    if (spec.target_return.has_value())
    {
        std::vector<long double> mu = widened(data.mu());
        equalities.gaps.push_back(wide_dot(mu, w) - *spec.target_return);
        equalities.rows.push_back(std::move(mu));
    }
    if (spec.target_volatility.has_value())
    {
        std::vector<long double> gradient = wide_covariance_times(data, w);
        long double volatility = std::sqrt(wide_dot(w, gradient));
        for (long double& entry : gradient)
        {
            entry /= volatility;
        }
        equalities.gaps.push_back(volatility - *spec.target_volatility);
        equalities.rows.push_back(std::move(gradient));
    }
    return equalities;
}

inline long double volatility_multiplier_of(const po::Spec& spec, const WideKkt& kkt)
{
    if (!spec.target_volatility.has_value())
    {
        return 0.0L;
    }
    return kkt.multipliers.back();
}

// grad f - sum_k lambda_k grad h_k at the point, every entry.
inline std::vector<long double> wide_stationarity(const po::Data& data, std::span<const double> budget,
                                                  const WideKkt& kkt, const WideEqualities& equalities)
{
    std::vector<long double> residual = wide_gradient_at(data, budget, kkt.weights);
    for (std::size_t k = 0; k < equalities.rows.size(); ++k)
    {
        for (std::size_t j = 0; j < residual.size(); ++j)
        {
            residual[j] -= kkt.multipliers[k] * equalities.rows[k][j];
        }
    }
    return residual;
}

// One Newton step of [H, -A'; A, 0] [dw; dlambda] = -[stationarity; gaps] over the free weights, the rows of
// A scaled to the size of H. False when the system is singular.
inline bool wide_newton_step(const po::Spec& spec, const po::Data& data, std::span<const double> budget,
                             WideKkt& kkt)
{
    std::size_t n = kkt.weights.size();
    std::size_t free_count = kkt.free.size();
    WideEqualities equalities = wide_equalities(spec, data, kkt.weights);
    std::size_t size = free_count + equalities.gaps.size();
    std::vector<long double> residual = wide_stationarity(data, budget, kkt, equalities);
    std::vector<long double> hessian =
        wide_lagrangian_hessian(data, budget, volatility_multiplier_of(spec, kkt), kkt.weights);
    std::vector<long double> matrix(size * size, 0.0L);
    std::vector<long double> rhs(size);
    long double scale = 0.0L;
    for (std::size_t a = 0; a < free_count; ++a)
    {
        for (std::size_t b = 0; b < free_count; ++b)
        {
            matrix[a * size + b] = hessian[kkt.free[a] * n + kkt.free[b]];
            scale = std::fmax(scale, std::fabs(matrix[a * size + b]));
        }
        rhs[a] = -residual[kkt.free[a]];
    }
    for (std::size_t k = 0; k < equalities.gaps.size(); ++k)
    {
        for (std::size_t a = 0; a < free_count; ++a)
        {
            matrix[a * size + free_count + k] = -scale * equalities.rows[k][kkt.free[a]];
            matrix[(free_count + k) * size + a] = scale * equalities.rows[k][kkt.free[a]];
        }
        rhs[free_count + k] = -scale * equalities.gaps[k];
    }
    if (!solve_wide(matrix, rhs))
    {
        return false;
    }
    for (std::size_t a = 0; a < free_count; ++a)
    {
        kkt.weights[kkt.free[a]] += rhs[a];
    }
    for (std::size_t k = 0; k < equalities.gaps.size(); ++k)
    {
        kkt.multipliers[k] += scale * rhs[free_count + k];
    }
    return true;
}

inline WideKkt wide_kkt_of(const po::Spec& spec, const po::Result& result)
{
    WideKkt kkt{widened(result.weights), {result.budget_multiplier}, {}};
    if (spec.target_return.has_value())
    {
        kkt.multipliers.push_back(result.target_return_multiplier);
    }
    if (spec.target_volatility.has_value())
    {
        kkt.multipliers.push_back(result.volatility_multiplier);
    }
    for (std::size_t j = 0; j < result.weights.size(); ++j)
    {
        if (result.bound_multipliers[j] == 0.0)
        {
            kkt.free.push_back(j);
        }
    }
    return kkt;
}

// How far a result is from the KKT point of its own active set, w_ref, computed from it in long double.
struct ForwardError
{
    // max_i |w_i - w_ref_i| / max_i |w_i|.
    double weights;
    // |f(w) - f(w_ref)|.
    double value;
    // The largest stationarity entry of w_ref over the free weights, in units of the rounding floor of the
    // result (doubles): far below 1 when w_ref is the more accurate point. NaN when a Newton system was
    // singular.
    double reference_over_floor;
    // The largest relative gap of an equality at w_ref.
    double reference_gap;
};

inline double reference_over_floor(const po::Spec& spec, const po::Data& data, std::span<const double> budget,
                                   const WideKkt& kkt)
{
    std::size_t n = kkt.weights.size();
    WideEqualities equalities = wide_equalities(spec, data, kkt.weights);
    std::vector<long double> residual = wide_stationarity(data, budget, kkt, equalities);
    std::vector<long double> hessian =
        wide_lagrangian_hessian(data, budget, volatility_multiplier_of(spec, kkt), kkt.weights);
    double worst = 0.0;
    for (std::size_t j : kkt.free)
    {
        long double moved = 0.0L;
        for (std::size_t l = 0; l < n; ++l)
        {
            moved += std::fabs(hessian[j * n + l] * kkt.weights[l]);
        }
        for (std::size_t k = 0; k < equalities.rows.size(); ++k)
        {
            moved += std::fabs(kkt.multipliers[k] * equalities.rows[k][j]);
        }
        long double floor = std::numeric_limits<double>::epsilon() * moved;
        worst = worst_of(worst, static_cast<double>(std::fabs(residual[j]) / floor));
    }
    return worst;
}

inline double reference_gap(const po::Spec& spec, const po::Data& data, const WideKkt& kkt)
{
    WideEqualities equalities = wide_equalities(spec, data, kkt.weights);
    double worst = static_cast<double>(std::fabs(equalities.gaps[0]));
    std::size_t next = 1;
    if (spec.target_return.has_value())
    {
        double scale = worst_of(std::fabs(*spec.target_return), largest_magnitude(data.mu()));
        worst = worst_of(worst, relative_to(static_cast<double>(std::fabs(equalities.gaps[next])), scale));
        ++next;
    }
    if (spec.target_volatility.has_value())
    {
        double gap = static_cast<double>(std::fabs(equalities.gaps[next]));
        worst = worst_of(worst, gap / *spec.target_volatility);
    }
    return worst;
}

// Six Newton steps in long double from the result, on its active set, with the reference gradient and its
// Hessian by differences. The steps are not counted on to converge: reference_over_floor says whether they
// did.
inline ForwardError forward_error(const po::Spec& spec, const po::Data& data, const po::Result& result)
{
    constexpr int kSteps = 6;
    constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> budget = budget_of(spec, data.asset_count());
    WideKkt kkt = wide_kkt_of(spec, result);
    for (int step = 0; step < kSteps; ++step)
    {
        if (!wide_newton_step(spec, data, budget, kkt))
        {
            return ForwardError{kNaN, kNaN, kNaN, kNaN};
        }
    }
    long double largest = 0.0L;
    long double distance = 0.0L;
    for (std::size_t i = 0; i < kkt.weights.size(); ++i)
    {
        largest = std::fmax(largest, std::fabs(static_cast<long double>(result.weights[i])));
        distance = std::fmax(distance, std::fabs(result.weights[i] - kkt.weights[i]));
    }
    long double change =
        wide_value_at(data, budget, widened(result.weights)) - wide_value_at(data, budget, kkt.weights);
    return ForwardError{static_cast<double>(distance / largest), static_cast<double>(std::fabs(change)),
                        reference_over_floor(spec, data, budget, kkt), reference_gap(spec, data, kkt)};
}

} // namespace po::testing
