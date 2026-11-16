#pragma once

#include <algorithm>
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
    // max_k |lambda_k - lambda_ref_k| / max_k |lambda_ref_k| over the multipliers of the equalities.
    double multipliers = 0.0;
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
    std::vector<long double> returned = kkt.multipliers;
    for (int step = 0; step < kSteps; ++step)
    {
        if (!wide_newton_step(spec, data, budget, kkt))
        {
            return ForwardError{kNaN, kNaN, kNaN, kNaN, kNaN};
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
    long double multiplier_scale = 0.0L;
    long double multiplier_distance = 0.0L;
    for (std::size_t k = 0; k < returned.size(); ++k)
    {
        multiplier_scale = std::fmax(multiplier_scale, std::fabs(kkt.multipliers[k]));
        multiplier_distance = std::fmax(multiplier_distance, std::fabs(returned[k] - kkt.multipliers[k]));
    }
    return ForwardError{static_cast<double>(distance / largest), static_cast<double>(std::fabs(change)),
                        reference_over_floor(spec, data, budget, kkt), reference_gap(spec, data, kkt),
                        static_cast<double>(multiplier_distance / multiplier_scale)};
}

// The rounding of each row of wide_kkt_system at the point of kkt: for a free weight j,
// eps ((abs(H) abs(w))_j + sum_k abs(lambda_k grad h_k)_j); for an equality,
// eps (sum_i abs(grad h_k)_i abs(w_i) + abs(target)), times scale.
inline std::vector<long double> wide_kkt_floor(const po::Spec& spec, const po::Data& data,
                                               std::span<const double> budget, const WideKkt& kkt,
                                               long double scale)
{
    constexpr long double kEpsilon = std::numeric_limits<double>::epsilon();
    std::size_t n = kkt.weights.size();
    WideEqualities equalities = wide_equalities(spec, data, kkt.weights);
    std::vector<long double> hessian =
        wide_lagrangian_hessian(data, budget, volatility_multiplier_of(spec, kkt), kkt.weights);
    std::vector<long double> floor;
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
        floor.push_back(kEpsilon * moved);
    }
    std::vector<long double> targets{1.0L};
    if (spec.target_return.has_value())
    {
        targets.push_back(*spec.target_return);
    }
    if (spec.target_volatility.has_value())
    {
        targets.push_back(*spec.target_volatility);
    }
    for (std::size_t k = 0; k < equalities.rows.size(); ++k)
    {
        long double terms = std::fabs(targets[k]);
        for (std::size_t i = 0; i < n; ++i)
        {
            terms += std::fabs(equalities.rows[k][i] * kkt.weights[i]);
        }
        floor.push_back(scale * kEpsilon * terms);
    }
    return floor;
}

// The a posteriori bound of the weight error of a result whose residual is at its rounding floor:
// max over the free weights of sum_l abs(K^-1)_jl floor_l, K the KKT matrix of its active set
// (wide_kkt_system) and floor its rounding per row (wide_kkt_floor), over max_i abs(w_i). abs(K^-1) is taken
// column by column, K x = e_l, to estimate how much K amplifies the rounding; nothing is solved with it. NaN
// when K is singular.
struct WeightErrorBound
{
    // From the rounding floor of each row: what any point at that floor may be off by.
    double floor;
    // From the residual of the result itself, sum_l abs(K^-1)_jl abs(r_l).
    double residual;
};

inline WeightErrorBound weight_error_bound(const po::Spec& spec, const po::Data& data,
                                           const po::Result& result)
{
    std::vector<double> budget = budget_of(spec, data.asset_count());
    WideKkt kkt = wide_kkt_of(spec, result);
    WideSystem system = wide_kkt_system(spec, data, budget, kkt);
    std::vector<long double> floor = wide_kkt_floor(spec, data, budget, kkt, system.scale);
    std::size_t size = floor.size();
    std::vector<long double> by_floor(kkt.free.size(), 0.0L);
    std::vector<long double> by_residual(kkt.free.size(), 0.0L);
    for (std::size_t l = 0; l < size; ++l)
    {
        std::vector<long double> matrix = system.matrix;
        std::vector<long double> column(size, 0.0L);
        column[l] = 1.0L;
        if (!solve_wide(matrix, column))
        {
            constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
            return WeightErrorBound{kNaN, kNaN};
        }
        for (std::size_t a = 0; a < by_floor.size(); ++a)
        {
            by_floor[a] += std::fabs(column[a]) * floor[l];
            by_residual[a] += std::fabs(column[a] * system.rhs[l]);
        }
    }
    long double largest = 0.0L;
    for (double weight : result.weights)
    {
        largest = std::fmax(largest, std::fabs(static_cast<long double>(weight)));
    }
    auto worst = [&](const std::vector<long double>& bound)
    { return static_cast<double>(*std::ranges::max_element(bound) / largest); };
    return WeightErrorBound{worst(by_floor), worst(by_residual)};
}

} // namespace po::testing
