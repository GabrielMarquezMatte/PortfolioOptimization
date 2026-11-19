#pragma once

#include <gtest/gtest-spi.h>
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/data.hpp"
#include "core/objective.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "fuzz/frontier_reference.hpp"
#include "kkt.hpp"
#include "random_covariance.hpp"
#include "user_kkt.hpp"

namespace po::testing
{

// long double is double under MSVC. The references below are accurate there too, because every sum that
// cancels is a compensated dot product; the checks that need more than double precision (the rounding floor
// of the stationarity, the refinement of a point beyond double) test this and are skipped where it is false.
inline constexpr bool kLongDoubleIsExtended = std::numeric_limits<long double>::digits > 53;

// a'b by Dot2 (Ogita, Rump and Oishi 2005): TwoProduct by fma and TwoSum, as accurate as a sum in twice the
// precision of T rounded once, whatever T is. The same algorithm as po::riskparity::compensated_dot, written
// here so that the reference does not share the code under test.
template <typename T, typename Left, typename Right>
T compensated_dot_of(std::size_t size, const Left& a, const Right& b)
{
    T sum = 0;
    T correction = 0;
    for (std::size_t i = 0; i < size; ++i)
    {
        T left = static_cast<T>(a[i]);
        T right = static_cast<T>(b[i]);
        T product = left * right;
        T product_error = std::fma(left, right, -product);
        T total = sum + product;
        T carried = total - sum;
        correction += (sum - (total - carried)) + (product - carried) + product_error;
        sum = total;
    }
    return sum + correction;
}

inline std::vector<long double> widened(std::span<const double> values)
{
    return std::vector<long double>(values.begin(), values.end());
}

inline std::vector<long double> wide_covariance_times(const po::Data& data, const std::vector<long double>& weights)
{
    std::size_t n = data.asset_count();
    std::span<const double> covariance = data.covariance();
    std::vector<long double> product(n, 0.0L);
    for (std::size_t i = 0; i < n; ++i)
    {
        product[i] = compensated_dot_of<long double>(n, covariance.subspan(i * n, n), weights);
    }
    return product;
}

inline long double wide_dot(const std::vector<long double>& a, const std::vector<long double>& b)
{
    return compensated_dot_of<long double>(a.size(), a, b);
}

// c_i = w_i (Sigma w)_i / (w' Sigma w), in long double.
inline std::vector<long double> wide_contributions(const po::Data& data, const std::vector<long double>& weights)
{
    std::vector<long double> product = wide_covariance_times(data, weights);
    long double variance = wide_dot(weights, product);
    std::vector<long double> contributions(weights.size());
    for (std::size_t i = 0; i < weights.size(); ++i)
    {
        contributions[i] = weights[i] * product[i] / variance;
    }
    return contributions;
}

inline std::vector<double> risk_contributions(const po::Data& data, std::span<const double> weights)
{
    std::vector<long double> wide = wide_contributions(data, widened(weights));
    return std::vector<double>(wide.begin(), wide.end());
}

// f(w) = sum (c_i - b_i)^2, plain loops in long double: independent of po::metrics.
inline double risk_parity_value(const po::Data& data, std::span<const double> budget, std::span<const double> w)
{
    std::vector<long double> contributions = wide_contributions(data, widened(w));
    long double total = 0.0L;
    for (std::size_t i = 0; i < contributions.size(); ++i)
    {
        long double gap = contributions[i] - budget[i];
        total += gap * gap;
    }
    return static_cast<double>(total);
}

// grad f = 2 J'(c - b) with J = (diag(Sigma w) + diag(w) Sigma - 2 c (Sigma w)') / (w' Sigma w), as the three
// terms 2 d_j (Sigma w)_j / v, 2 (Sigma (w o d))_j / v and -4 (Sigma w)_j (c'd) / v, with d = c - b and
// v = w' Sigma w.
inline std::array<std::vector<long double>, 3> wide_gradient_terms_at(const po::Data& data,
                                                                      std::span<const double> budget,
                                                                      const std::vector<long double>& w)
{
    std::vector<long double> product = wide_covariance_times(data, w);
    long double variance = wide_dot(w, product);
    std::size_t n = w.size();
    std::vector<long double> gap(n);
    std::vector<long double> scaled(n);
    std::vector<long double> contributions(n);
    for (std::size_t i = 0; i < n; ++i)
    {
        contributions[i] = w[i] * product[i] / variance;
        gap[i] = contributions[i] - budget[i];
        scaled[i] = w[i] * gap[i];
    }
    long double weighted_gap = wide_dot(contributions, gap);
    std::vector<long double> spread = wide_covariance_times(data, scaled);
    std::array<std::vector<long double>, 3> terms{std::vector<long double>(n), std::vector<long double>(n),
                                                  std::vector<long double>(n)};
    long double factor = 2.0L / variance;
    for (std::size_t i = 0; i < n; ++i)
    {
        terms[0][i] = factor * gap[i] * product[i];
        terms[1][i] = factor * spread[i];
        terms[2][i] = -2.0L * factor * product[i] * weighted_gap;
    }
    return terms;
}

inline std::array<std::vector<long double>, 3> wide_gradient_terms(const po::Data& data,
                                                                   std::span<const double> budget,
                                                                   std::span<const double> weights)
{
    return wide_gradient_terms_at(data, budget, widened(weights));
}

inline std::vector<long double> wide_gradient_at(const po::Data& data, std::span<const double> budget,
                                                 const std::vector<long double>& w)
{
    std::array<std::vector<long double>, 3> terms = wide_gradient_terms_at(data, budget, w);
    std::vector<long double> gradient(w.size());
    for (std::size_t i = 0; i < gradient.size(); ++i)
    {
        gradient[i] = terms[0][i] + terms[1][i] + terms[2][i];
    }
    return gradient;
}

// f at weights in long double.
inline long double wide_value_at(const po::Data& data, std::span<const double> budget,
                                 const std::vector<long double>& w)
{
    std::vector<long double> contributions = wide_contributions(data, w);
    long double total = 0.0L;
    for (std::size_t i = 0; i < contributions.size(); ++i)
    {
        long double gap = contributions[i] - budget[i];
        total += gap * gap;
    }
    return total;
}

// sqrt(w' Sigma w) in long double: the plain double product loses eps cond(Sigma).
inline double wide_volatility(const po::Data& data, std::span<const double> weights)
{
    std::vector<long double> w = widened(weights);
    return static_cast<double>(std::sqrt(wide_dot(w, wide_covariance_times(data, w))));
}

inline std::vector<std::vector<double>> risk_parity_gradient_terms(const po::Data& data,
                                                                   std::span<const double> budget,
                                                                   std::span<const double> w)
{
    std::array<std::vector<long double>, 3> wide = wide_gradient_terms(data, budget, w);
    std::vector<std::vector<double>> terms;
    for (const std::vector<long double>& term : wide)
    {
        terms.emplace_back(term.begin(), term.end());
    }
    return terms;
}

inline std::vector<double> risk_parity_gradient(const po::Data& data, std::span<const double> budget,
                                                std::span<const double> w)
{
    std::array<std::vector<long double>, 3> terms = wide_gradient_terms(data, budget, w);
    std::vector<double> gradient(w.size());
    for (std::size_t i = 0; i < gradient.size(); ++i)
    {
        gradient[i] = static_cast<double>(terms[0][i] + terms[1][i] + terms[2][i]);
    }
    return gradient;
}

// With w1 = x, w2 = 1 - x and k = rho sigma1 sigma2, c1 b2 = c2 b1 is A x^2 + B x + C = 0 with
// A = b2 (sigma1^2 - k) - b1 (sigma2^2 - k), B = b2 k - b1 (k - 2 sigma2^2), C = -b1 sigma2^2. C < 0 and the
// polynomial is positive at x = 1, so exactly one root lies in (0, 1); c1 = b1 follows since c1 + c2 = 1.
inline std::vector<double> two_asset_risk_parity(double sigma1, double sigma2, double rho, double b1)
{
    double b2 = 1.0 - b1;
    double cross = rho * sigma1 * sigma2;
    double a = b2 * (sigma1 * sigma1 - cross) - b1 * (sigma2 * sigma2 - cross);
    double b = b2 * cross - b1 * (cross - 2.0 * sigma2 * sigma2);
    double c = -b1 * sigma2 * sigma2;
    double discriminant = std::fmax(b * b - 4.0 * a * c, 0.0);
    double q = -0.5 * (b + std::copysign(std::sqrt(discriminant), b));
    for (double root : {q / a, c / q})
    {
        if (root > 0.0 && root < 1.0)
        {
            return {root, 1.0 - root};
        }
    }
    return {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN()};
}

struct ScanResult
{
    double value;
    std::vector<double> weights;
};

using Triple = std::array<double, 3>;

// The points of a scan over [first, last], ends included.
constexpr std::size_t kScanPoints = 20001;

// A one-parameter family of points of the feasible set of n = 3, over the parameter interval [first, last], scanned
// at `points` evenly spaced values of the parameter.
struct ScanCurve
{
    std::function<std::vector<double>(double)> point;
    double first;
    double last;
    std::size_t points = kScanPoints;
};

inline Triple cross_product(const Triple& a, const Triple& b)
{
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

inline double triple_dot(const Triple& a, const Triple& b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline Triple covariance_row(const po::Data& data, std::size_t row)
{
    std::span<const double> covariance = data.covariance();
    return {covariance[row * 3], covariance[row * 3 + 1], covariance[row * 3 + 2]};
}

// Sigma^-1 1 of a 3 x 3 matrix, from the cross products of its rows.
inline Triple inverse_covariance_ones(const po::Data& data)
{
    Triple row0 = covariance_row(data, 0);
    Triple row1 = covariance_row(data, 1);
    Triple row2 = covariance_row(data, 2);
    Triple column0 = cross_product(row1, row2);
    Triple column1 = cross_product(row2, row0);
    Triple column2 = cross_product(row0, row1);
    double determinant = triple_dot(row0, column0);
    Triple solution{};
    for (std::size_t k = 0; k < 3; ++k)
    {
        solution[k] = (column0[k] + column1[k] + column2[k]) / determinant;
    }
    return solution;
}

inline double covariance_form(const po::Data& data, const Triple& a, const Triple& b)
{
    double total = 0.0;
    for (std::size_t i = 0; i < 3; ++i)
    {
        total += a[i] * triple_dot(covariance_row(data, i), b);
    }
    return total;
}

// A point with 1'p = 1 and mu'p = target, as a 1 + b mu.
inline Triple return_line_base(const Triple& mu, double target)
{
    double sum = mu[0] + mu[1] + mu[2];
    double squares = triple_dot(mu, mu);
    double determinant = 3.0 * squares - sum * sum;
    double a = (squares - sum * target) / determinant;
    double b = (3.0 * target - sum) / determinant;
    return {a + b * mu[0], a + b * mu[1], a + b * mu[2]};
}

inline void narrow_to_box(const po::Spec& spec, const Triple& base, const Triple& direction, ScanCurve& curve)
{
    for (std::size_t k = 0; k < 3; ++k)
    {
        if (direction[k] == 0.0)
        {
            continue;
        }
        double low = (lower_of(spec, k) - base[k]) / direction[k];
        double high = (upper_of(spec, k) - base[k]) / direction[k];
        if (direction[k] < 0.0)
        {
            std::swap(low, high);
        }
        curve.first = std::fmax(curve.first, low);
        curve.last = std::fmin(curve.last, high);
    }
}

// The line {1'w = 1, mu'w = r}, p + t d with d = 1 x mu, narrowed to the finite bounds of the box.
inline ScanCurve return_line_curve(const po::Spec& spec, const po::Data& data)
{
    Triple mu{data.mu()[0], data.mu()[1], data.mu()[2]};
    Triple base = return_line_base(mu, *spec.target_return);
    Triple direction = cross_product({1.0, 1.0, 1.0}, mu);
    double largest = std::fmax(std::fabs(direction[0]), std::fmax(std::fabs(direction[1]), std::fabs(direction[2])));
    double reach = 5.0 / largest;
    ScanCurve curve{[base, direction](double t)
                    {
                        return std::vector<double>{base[0] + t * direction[0], base[1] + t * direction[1],
                                                   base[2] + t * direction[2]};
                    },
                    -reach, reach};
    narrow_to_box(spec, base, direction, curve);
    return curve;
}

// {1'w = 1, w'Sigma w = s^2}: the minimum-variance point plus B y, B an orthonormal basis of 1'u = 0 and y such that
// u'Sigma u = radius^2 (Cholesky of B'Sigma B), by a point z = radius c of the circle of radius radius, L'y = z.
struct EllipseFrame
{
    Triple center;
    Triple first_axis;
    Triple second_axis;
    double l11;
    double l21;
    double l22;
    double radius;
};

inline EllipseFrame ellipse_frame(const po::Data& data, const Triple& center, double radius)
{
    Triple first_axis{1.0 / std::sqrt(2.0), -1.0 / std::sqrt(2.0), 0.0};
    Triple second_axis{1.0 / std::sqrt(6.0), 1.0 / std::sqrt(6.0), -2.0 / std::sqrt(6.0)};
    double m11 = covariance_form(data, first_axis, first_axis);
    double m21 = covariance_form(data, second_axis, first_axis);
    double m22 = covariance_form(data, second_axis, second_axis);
    double l11 = std::sqrt(m11);
    double l21 = m21 / l11;
    return EllipseFrame{center, first_axis, second_axis, l11, l21, std::sqrt(m22 - l21 * l21), radius};
}

// The unit circle by x in [-1, 3], rational so that it has the same bits on every platform (no sin or cos):
// c = ((1 - t^2) / (1 + t^2), 2 t / (1 + t^2)) at t = x for x <= 1, the half with c_1 >= 0 at the angle 2 atan(t),
// and -c at t = x - 2 for x > 1, the other half. Both pieces meet at x = 1 and the circle closes at x = 3; beyond
// [-1, 3] the point stays on the circle. The angle moves between 1 and 2 radians per unit of x.
inline std::array<double, 2> unit_circle_point(double x)
{
    double t = x;
    double sign = 1.0;
    if (x > 1.0)
    {
        t = x - 2.0;
        sign = -1.0;
    }
    double denominator = 1.0 + t * t;
    return {sign * (1.0 - t * t) / denominator, sign * 2.0 * t / denominator};
}

inline std::vector<double> ellipse_point(const EllipseFrame& frame, double x)
{
    std::array<double, 2> circle = unit_circle_point(x);
    double second = frame.radius * circle[1] / frame.l22;
    double first = (frame.radius * circle[0] - frame.l21 * second) / frame.l11;
    std::vector<double> weights(3);
    for (std::size_t k = 0; k < 3; ++k)
    {
        weights[k] = frame.center[k] + first * frame.first_axis[k] + second * frame.second_axis[k];
    }
    return weights;
}

// Twice the points of the other scans: a step of 4 / 40000 in x is at most 2e-4 radians, below the 2 pi / 20000 of
// a scan by angle with 20001 points.
inline ScanCurve ellipse_curve(const EllipseFrame& frame)
{
    return ScanCurve{[frame](double x) { return ellipse_point(frame, x); }, -1.0, 3.0, 2 * kScanPoints - 1};
}

// The curve of the feasible set; none when it is not a curve (neither or both of the two targets) or empty
// (s at or below the minimum volatility).
inline std::optional<ScanCurve> scan_curve_of(const po::Spec& spec, const po::Data& data)
{
    if (spec.target_return.has_value() == spec.target_volatility.has_value())
    {
        return std::nullopt;
    }
    if (spec.target_return.has_value())
    {
        return return_line_curve(spec, data);
    }
    Triple solution = inverse_covariance_ones(data);
    double total = solution[0] + solution[1] + solution[2];
    double radius_squared = *spec.target_volatility * *spec.target_volatility - 1.0 / total;
    if (radius_squared <= 0.0)
    {
        return std::nullopt;
    }
    Triple center{solution[0] / total, solution[1] / total, solution[2] / total};
    return ellipse_curve(ellipse_frame(data, center, std::sqrt(radius_squared)));
}

constexpr double kScanBoxTolerance = 1e-12;

// The value of the point at x, +inf outside the box.
inline ScanResult scan_point(const po::Spec& spec, const po::Data& data, std::span<const double> budget,
                             const ScanCurve& curve, double x)
{
    std::vector<double> weights = curve.point(x);
    for (std::size_t k = 0; k < weights.size(); ++k)
    {
        if (weights[k] < lower_of(spec, k) - kScanBoxTolerance || weights[k] > upper_of(spec, k) + kScanBoxTolerance)
        {
            return ScanResult{std::numeric_limits<double>::infinity(), {}};
        }
    }
    double value = risk_parity_value(data, budget, weights);
    return ScanResult{value, std::move(weights)};
}

inline double golden_section_minimum(const std::function<double(double)>& value, double low, double high)
{
    constexpr double kInverseGolden = 0.6180339887498949;
    constexpr int kMaxIterations = 200;
    double inner_low = high - kInverseGolden * (high - low);
    double inner_high = low + kInverseGolden * (high - low);
    double value_low = value(inner_low);
    double value_high = value(inner_high);
    for (int iteration = 0; iteration < kMaxIterations; ++iteration)
    {
        if (high - low <= 1e-14 * (1.0 + std::fabs(low) + std::fabs(high)))
        {
            break;
        }
        if (value_low < value_high)
        {
            high = inner_high;
            inner_high = inner_low;
            value_high = value_low;
            inner_low = high - kInverseGolden * (high - low);
            value_low = value(inner_low);
            continue;
        }
        low = inner_low;
        inner_low = inner_high;
        value_low = value_high;
        inner_high = low + kInverseGolden * (high - low);
        value_high = value(inner_high);
    }
    return 0.5 * (low + high);
}

inline ScanResult scan_curve_minimum(const po::Spec& spec, const po::Data& data, std::span<const double> budget,
                                     const ScanCurve& curve)
{
    ScanResult best{std::numeric_limits<double>::infinity(), {}};
    if (curve.first > curve.last)
    {
        best.value = std::numeric_limits<double>::quiet_NaN();
        return best;
    }
    double step = (curve.last - curve.first) / static_cast<double>(curve.points - 1);
    double best_x = curve.first;
    for (std::size_t k = 0; k < curve.points; ++k)
    {
        double x = curve.first + step * static_cast<double>(k);
        ScanResult candidate = scan_point(spec, data, budget, curve, x);
        if (candidate.value < best.value)
        {
            best = std::move(candidate);
            best_x = x;
        }
    }
    if (std::isinf(best.value))
    {
        best.value = std::numeric_limits<double>::quiet_NaN();
        return best;
    }
    auto value_at = [&](double x) { return scan_point(spec, data, budget, curve, x).value; };
    double refined_x = golden_section_minimum(value_at, best_x - step, best_x + step);
    ScanResult refined = scan_point(spec, data, budget, curve, refined_x);
    if (refined.value < best.value)
    {
        best = std::move(refined);
    }
    return best;
}

// Global minimum by scan, n = 3 only. Feasible set a curve: budget + target_return (a line) or budget +
// target_volatility (an ellipse around the minimum-variance point, by a point of a circle). The curve's points
// (20001, 40001 on the ellipse), then golden-section on the best bracket; points outside the box are discarded. The
// value is NaN when no scanned point is feasible or the feasible set is not a curve.
inline ScanResult scan_reference(const po::Spec& spec, const po::Data& data, std::span<const double> budget)
{
    std::optional<ScanCurve> curve = scan_curve_of(spec, data);
    if (!curve.has_value())
    {
        return ScanResult{std::numeric_limits<double>::quiet_NaN(), {}};
    }
    return scan_curve_minimum(spec, data, budget, *curve);
}

// The scan's value moved to the feasible set to first order: its point sits on the curve only to the rounding of
// the curve (a gap of up to 2e-10 of s measured) and inside the box only to kScanBoxTolerance, and f changes by
// lambda_k h_k off the curve and by beta_j h_j past an active bound, with the multipliers of the result. Without
// this a scan point just inside the ellipse, or just past a bound the result sits on, beats an exact solution.
inline double scan_value_at_target(const po::Spec& spec, const po::Data& data, const po::Result& result,
                                   const ScanResult& scan)
{
    double total = scan.weights[0] + scan.weights[1] + scan.weights[2];
    double value = scan.value - result.budget_multiplier * (total - 1.0);
    if (spec.target_return.has_value())
    {
        double mean = compensated_dot_of<double>(3, data.mu(), scan.weights);
        value -= result.target_return_multiplier * (mean - *spec.target_return);
    }
    if (spec.target_volatility.has_value())
    {
        double volatility = wide_volatility(data, scan.weights);
        value -= result.volatility_multiplier * (volatility - *spec.target_volatility);
    }
    for (std::size_t j = 0; j < 3; ++j)
    {
        double beta = result.bound_multipliers[j];
        double below = std::fmin(scan.weights[j] - lower_of(spec, j), 0.0);
        double above = std::fmax(scan.weights[j] - upper_of(spec, j), 0.0);
        if (beta > 0.0)
        {
            value -= beta * below;
        }
        if (beta < 0.0)
        {
            value -= beta * above;
        }
    }
    return value;
}

struct RiskParityShape
{
    std::size_t n;
    double condition;
    // 1 (annual) or 1e-4 (daily) on Sigma; its square root on mu.
    double scale;
    // Uneven: log-uniform b_i in [1e-6, 1], normalized.
    enum class Budget
    {
        Equal,
        Uneven
    } budget;
    // Capped: lower 0, upper 1.2 / n.
    enum class Box
    {
        None,
        LongOnly,
        Capped
    } box;
    // r = a convex combination of two means; under Capped, moved inside the achievable returns.
    bool target_return;
    // s between sigma_min and 1.5 sigma_min of the linear set.
    bool target_volatility;
    std::uint64_t seed;
};

struct GeneratedRiskParity
{
    po::Spec spec;
    po::Data data;
};

inline std::vector<double> shape_budget(const RiskParityShape& shape)
{
    std::vector<double> budget(shape.n, 1.0 / static_cast<double>(shape.n));
    if (shape.budget == RiskParityShape::Budget::Equal)
    {
        return budget;
    }
    double total = 0.0;
    std::vector<double> exponents = random_vector(shape.n, shape.seed + 700, -6.0, 0.0);
    for (std::size_t i = 0; i < shape.n; ++i)
    {
        budget[i] = exp10(exponents[i]);
        total += budget[i];
    }
    for (double& entry : budget)
    {
        entry /= total;
    }
    return budget;
}

// Means in [0.03, 0.15] at the annual scale.
inline std::vector<double> shape_means(const RiskParityShape& shape)
{
    double root = std::sqrt(shape.scale);
    return random_vector(shape.n, shape.seed + 500, 0.03 * root, 0.15 * root);
}

inline po::Spec shape_bounds(const RiskParityShape& shape)
{
    po::Spec spec;
    if (shape.box == RiskParityShape::Box::None)
    {
        return spec;
    }
    spec.lower = 0.0;
    if (shape.box == RiskParityShape::Box::Capped)
    {
        spec.upper = 1.2 / static_cast<double>(shape.n);
    }
    return spec;
}

// The return of the weights that fill the means in the given order up to the cap, total 1.
inline double filled_return(const std::vector<double>& ordered_means, double cap)
{
    double remaining = 1.0;
    double total = 0.0;
    for (double mean : ordered_means)
    {
        double weight = std::fmin(cap, remaining);
        total += weight * mean;
        remaining -= weight;
    }
    return total;
}

// The target moved into the middle 80% of the returns reachable with weights in [0, cap] and total 1.
inline double inside_capped_returns(std::vector<double> means, double cap, double target)
{
    std::ranges::sort(means);
    double lowest = filled_return(means, cap);
    std::ranges::reverse(means);
    double highest = filled_return(means, cap);
    double margin = 0.1 * (highest - lowest);
    return std::clamp(target, lowest + margin, highest - margin);
}

inline double shape_target_return(const RiskParityShape& shape, const std::vector<double>& means)
{
    std::vector<double> draws = random_vector(3, shape.seed + 800, 0.0, 1.0);
    double count = static_cast<double>(shape.n);
    std::size_t first = std::min(shape.n - 1, static_cast<std::size_t>(draws[0] * count));
    std::size_t second = (first + 1 + static_cast<std::size_t>(draws[1] * (count - 1.0))) % shape.n;
    double target = draws[2] * means[first] + (1.0 - draws[2]) * means[second];
    if (shape.box != RiskParityShape::Box::Capped)
    {
        return target;
    }
    return inside_capped_returns(means, 1.2 / count, target);
}

inline double shape_target_volatility(const po::Spec& spec, const po::Data& data, const RiskParityShape& shape)
{
    po::Result minimum = frontier_point(spec, data, spec.target_return);
    if (minimum.status != po::SolveStatus::Optimal)
    {
        throw std::runtime_error("the linear constraints of a risk parity shape are infeasible");
    }
    double fraction = random_vector(1, shape.seed + 900, 0.0, 1.0)[0];
    return portfolio_volatility(data, minimum.weights) * (1.0 + 0.5 * fraction);
}

inline GeneratedRiskParity generate(const RiskParityShape& shape)
{
    std::vector<double> covariance = random_covariance(shape.n, shape.seed, shape.condition);
    for (double& entry : covariance)
    {
        entry *= shape.scale;
    }
    std::vector<double> means = shape_means(shape);
    po::Data data = po::Data::create(shape.n, means, covariance).value();
    po::Spec spec = shape_bounds(shape);
    spec.objective = po::Objective::RiskParity;
    spec.risk_budget = shape_budget(shape);
    if (shape.target_return)
    {
        spec.target_return = shape_target_return(shape, means);
    }
    if (shape.target_volatility)
    {
        spec.target_volatility = shape_target_volatility(spec, data, shape);
    }
    return GeneratedRiskParity{std::move(spec), std::move(data)};
}

// Shape k takes combination k mod 192 of cond(Sigma) in {1e2, 1e4, 1e6, 1e8}, the two scales, the two budgets, the
// three boxes and the two targets, so any 192 consecutive shapes cover every combination; n (3 to largest_n) and
// the seed are drawn from seed.
inline std::vector<RiskParityShape> risk_parity_shapes(std::size_t count, std::uint64_t seed,
                                                       std::size_t largest_n = 24)
{
    constexpr std::array<double, 4> kConditions = {1e2, 1e4, 1e6, 1e8};
    constexpr std::array<double, 2> kScales = {1.0, 1e-4};
    std::mt19937_64 engine(seed);
    std::vector<RiskParityShape> shapes;
    shapes.reserve(count);
    for (std::size_t k = 0; k < count; ++k)
    {
        RiskParityShape shape{};
        shape.condition = kConditions[k % 4];
        shape.scale = kScales[(k / 4) % 2];
        shape.budget = static_cast<RiskParityShape::Budget>((k / 8) % 2);
        shape.box = static_cast<RiskParityShape::Box>((k / 16) % 3);
        shape.target_return = (k / 48) % 2 == 1;
        shape.target_volatility = (k / 96) % 2 == 1;
        shape.n = 3 + static_cast<std::size_t>(engine() % (largest_n - 2));
        shape.seed = engine();
        shapes.push_back(shape);
    }
    return shapes;
}

inline ::testing::Message describe(const RiskParityShape& shape)
{
    ::testing::Message message;
    message << "n " << shape.n << " cond " << shape.condition << " scale " << shape.scale << " budget "
            << static_cast<int>(shape.budget) << " box " << static_cast<int>(shape.box) << " target_return "
            << shape.target_return << " target_volatility " << shape.target_volatility << " seed " << shape.seed;
    return message;
}

inline std::vector<double> budget_of(const po::Spec& spec, std::size_t n)
{
    if (spec.risk_budget.empty())
    {
        return std::vector<double>(n, 1.0 / static_cast<double>(n));
    }
    return spec.risk_budget;
}

// The equalities at the weights: the budget, the target return and the target volatility, those of spec.
inline void expect_risk_parity_feasible(const po::Spec& spec, const po::Data& data, const po::Result& result,
                                        double tolerance)
{
    double total = 0.0;
    for (double weight : result.weights)
    {
        total += weight;
    }
    EXPECT_NEAR(total, 1.0, tolerance);
    if (spec.target_return.has_value())
    {
        double gap = std::fabs(portfolio_return(data, result.weights) - *spec.target_return);
        double scale = worst_of(std::fabs(*spec.target_return), largest_magnitude(data.mu()));
        EXPECT_LE(relative_to(gap, scale), tolerance);
    }
    if (spec.target_volatility.has_value())
    {
        double gap = std::fabs(wide_volatility(data, result.weights) - *spec.target_volatility);
        EXPECT_LE(relative_to(gap, *spec.target_volatility), tolerance);
    }
}

// Finite multipliers for the constraints of spec, NaN for the absent ones.
inline void expect_risk_parity_multiplier_pattern(const po::Spec& spec, const po::Result& result)
{
    EXPECT_TRUE(std::isfinite(result.budget_multiplier));
    EXPECT_EQ(std::isfinite(result.target_return_multiplier), spec.target_return.has_value());
    EXPECT_EQ(std::isnan(result.target_return_multiplier), !spec.target_return.has_value());
    EXPECT_EQ(std::isfinite(result.volatility_multiplier), spec.target_volatility.has_value());
    EXPECT_EQ(std::isnan(result.volatility_multiplier), !spec.target_volatility.has_value());
    EXPECT_TRUE(std::isnan(result.leverage_multiplier));
}

// The relative distance of s to the minimum volatility within which a target volatility is that minimum: the
// feasibility tolerance of the QP that computes it.
inline constexpr double kMinimumVolatilityBand = 1e-12;

// A volatility multiplier of -inf says that s is the minimum volatility over the linear rows and the bounds,
// where the only feasible point is the minimum-variance portfolio, the gradients of the constraints are
// dependent and no KKT multipliers exist in general: every other multiplier is NaN. Accepted only there: s
// within kMinimumVolatilityBand of the volatility of MIN_VARIANCE (po::solve, with the bounds and the target
// return of spec) and w that portfolio.
inline void expect_minimum_volatility_point(const po::Spec& spec, const po::Data& data,
                                            const po::Result& result, double tolerance)
{
    ASSERT_TRUE(spec.target_volatility.has_value());
    po::Result minimum = frontier_point(spec, data, spec.target_return);
    ASSERT_EQ(minimum.status, po::SolveStatus::Optimal);
    double target = *spec.target_volatility;
    double distance = std::fabs(target - wide_volatility(data, minimum.weights));
    EXPECT_LE(distance, kMinimumVolatilityBand * target) << "tolerance: s is not the minimum volatility";
    for (std::size_t k = 0; k < result.weights.size(); ++k)
    {
        EXPECT_NEAR(result.weights[k], minimum.weights[k], tolerance) << "weight " << k;
        EXPECT_TRUE(std::isnan(result.bound_multipliers[k])) << "bound multiplier " << k;
        EXPECT_GE(result.weights[k], lower_of(spec, k) - tolerance);
        EXPECT_LE(result.weights[k], upper_of(spec, k) + tolerance);
    }
    EXPECT_TRUE(std::isnan(result.budget_multiplier));
    EXPECT_TRUE(std::isnan(result.target_return_multiplier));
    EXPECT_TRUE(std::isnan(result.leverage_multiplier));
    expect_risk_parity_feasible(spec, data, result, tolerance);
}

// w'grad f = 0, since the contributions are homogeneous of degree 0, relative to sum |w_i| times the
// largest term of grad f (the scale of the stationarity check). Not relative to sum |w_i grad f_i|: at a KKT
// point with only the budget row and bounds active at w_i = 0, grad f_i = budget + beta_i and
// beta_i w_i = 0, so that sum is |budget| sum |w_i|, and budget = -sum_i beta_i w_i is itself zero up to the
// KKT residual.
inline void expect_homogeneity_identity(std::span<const double> weights, const ObjectiveGradient& gradient,
                                        double tolerance)
{
    double product = dot(weights, gradient.value);
    double gross = 0.0;
    for (double weight : weights)
    {
        gross += std::fabs(weight);
    }
    EXPECT_LE(relative_to(std::fabs(product), gross * gradient.scale), tolerance);
}

// Every term of grad f is proportional to the gap c - b, so at a minimizer with f = 0 they all vanish together and
// a residual relative to their own magnitude compares rounding with rounding. Only at that floor, a gap of at most
// the tolerance, is the scale of the problem used instead; at any larger gap the terms are the scale, so the band
// is not relaxed by 1 / gap. This is the largest |grad f_j| that a
// gap of at most 1 in every entry can produce, 2 (|p_j| + sum_i |Sigma_ji w_i| + 2 |p_j| sum_i |c_i|) / v with
// p = Sigma w and v = w' Sigma w: the scale of the problem, not of its current gap.
inline double unit_gap_gradient_scale(const po::Data& data, std::span<const double> weights)
{
    std::vector<long double> w = widened(weights);
    std::vector<long double> product = wide_covariance_times(data, w);
    long double variance = wide_dot(w, product);
    std::size_t n = w.size();
    std::span<const double> covariance = data.covariance();
    long double contribution_norm = 0.0L;
    for (std::size_t i = 0; i < n; ++i)
    {
        contribution_norm += std::fabs(w[i] * product[i] / variance);
    }
    long double largest = 0.0L;
    for (std::size_t j = 0; j < n; ++j)
    {
        long double row_sum = 0.0L;
        for (std::size_t i = 0; i < n; ++i)
        {
            row_sum += std::fabs(static_cast<long double>(covariance[j * n + i]) * w[i]);
        }
        long double bound = 2.0L * (std::fabs(product[j]) * (1.0L + 2.0L * contribution_norm) + row_sum) / variance;
        largest = std::fmax(largest, bound);
    }
    return static_cast<double>(largest);
}

inline double contribution_gap(const po::Data& data, std::span<const double> budget,
                               std::span<const double> weights)
{
    std::vector<double> contributions = risk_contributions(data, weights);
    double gap = 0.0;
    for (std::size_t i = 0; i < contributions.size(); ++i)
    {
        gap = std::fmax(gap, std::fabs(contributions[i] - budget[i]));
    }
    return gap;
}

// The Hessian of f at w (n x n, row-major) by central differences of the reference gradient. The step is
// kHessianStep of the largest weight, small against the relative distance of about 1 / cond(Sigma) over
// which the gradient turns, and large against the rounding of the gradient in long double. Measured against
// the analytic Hessian of the library within 1e-6 at cond(Sigma) up to 1e6
// (RiskParityReference.HessianByDifferencesMatchesTheObjective).
inline constexpr long double kHessianStep = 1e-12L;

inline std::vector<long double> wide_hessian(const po::Data& data, std::span<const double> budget,
                                             const std::vector<long double>& w)
{
    std::size_t n = w.size();
    long double largest = 0.0L;
    for (long double weight : w)
    {
        largest = std::fmax(largest, std::fabs(weight));
    }
    std::vector<long double> hessian(n * n);
    for (std::size_t j = 0; j < n; ++j)
    {
        std::vector<long double> above = w;
        std::vector<long double> below = w;
        above[j] += kHessianStep * largest;
        below[j] -= kHessianStep * largest;
        std::vector<long double> gradient_above = wide_gradient_at(data, budget, above);
        std::vector<long double> gradient_below = wide_gradient_at(data, budget, below);
        for (std::size_t i = 0; i < n; ++i)
        {
            hessian[i * n + j] = (gradient_above[i] - gradient_below[i]) / (above[j] - below[j]);
        }
    }
    return hessian;
}

// hess f - volatility hess vol, with hess vol = Sigma / vol - p p' / vol^3 and p = Sigma w: the derivative in
// w of the stationarity residual. Without a volatility multiplier (NaN or 0) it is hess f.
inline std::vector<long double> wide_lagrangian_hessian(const po::Data& data, std::span<const double> budget,
                                                        long double volatility_multiplier,
                                                        const std::vector<long double>& w)
{
    std::vector<long double> hessian = wide_hessian(data, budget, w);
    if (std::isnan(volatility_multiplier) || volatility_multiplier == 0.0L)
    {
        return hessian;
    }
    std::size_t n = w.size();
    std::vector<long double> product = wide_covariance_times(data, w);
    long double volatility = std::sqrt(wide_dot(w, product));
    long double cube = volatility * volatility * volatility;
    for (std::size_t j = 0; j < n; ++j)
    {
        for (std::size_t l = 0; l < n; ++l)
        {
            long double covariance = data.covariance()[j * n + l];
            long double curvature = covariance / volatility - product[j] * product[l] / cube;
            hessian[j * n + l] -= volatility_multiplier * curvature;
        }
    }
    return hessian;
}

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

// [H, -A'; A, 0] over the free weights of kkt and -[stationarity; gaps], the rows of A scaled to the size
// of H by scale (size x size, row-major).
struct WideSystem
{
    std::vector<long double> matrix;
    std::vector<long double> rhs;
    long double scale = 0.0L;
};

inline WideSystem wide_kkt_system(const po::Spec& spec, const po::Data& data, std::span<const double> budget,
                                  const WideKkt& kkt)
{
    std::size_t n = kkt.weights.size();
    std::size_t free_count = kkt.free.size();
    WideEqualities equalities = wide_equalities(spec, data, kkt.weights);
    std::size_t size = free_count + equalities.gaps.size();
    std::vector<long double> residual = wide_stationarity(data, budget, kkt, equalities);
    std::vector<long double> hessian =
        wide_lagrangian_hessian(data, budget, volatility_multiplier_of(spec, kkt), kkt.weights);
    WideSystem system{std::vector<long double>(size * size, 0.0L), std::vector<long double>(size), 0.0L};
    for (std::size_t a = 0; a < free_count; ++a)
    {
        for (std::size_t b = 0; b < free_count; ++b)
        {
            system.matrix[a * size + b] = hessian[kkt.free[a] * n + kkt.free[b]];
            system.scale = std::fmax(system.scale, std::fabs(system.matrix[a * size + b]));
        }
        system.rhs[a] = -residual[kkt.free[a]];
    }
    for (std::size_t k = 0; k < equalities.gaps.size(); ++k)
    {
        for (std::size_t a = 0; a < free_count; ++a)
        {
            system.matrix[a * size + free_count + k] = -system.scale * equalities.rows[k][kkt.free[a]];
            system.matrix[(free_count + k) * size + a] = system.scale * equalities.rows[k][kkt.free[a]];
        }
        system.rhs[free_count + k] = -system.scale * equalities.gaps[k];
    }
    return system;
}

// One Newton step on wide_kkt_system. False when the system is singular.
inline bool wide_newton_step(const po::Spec& spec, const po::Data& data, std::span<const double> budget,
                             WideKkt& kkt)
{
    std::size_t free_count = kkt.free.size();
    WideSystem system = wide_kkt_system(spec, data, budget, kkt);
    if (!solve_wide(system.matrix, system.rhs))
    {
        return false;
    }
    for (std::size_t a = 0; a < free_count; ++a)
    {
        kkt.weights[kkt.free[a]] += system.rhs[a];
    }
    for (std::size_t k = 0; k + free_count < system.rhs.size(); ++k)
    {
        kkt.multipliers[k] += system.scale * system.rhs[free_count + k];
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

// The constraint terms of the stationarity in long double: their sum, the sum of their magnitudes per entry,
// and the largest single one.
struct WideConstraintTerms
{
    std::vector<long double> total;
    std::vector<long double> magnitude;
    long double scale = 0.0L;
};

// A NaN term makes the scale NaN for good: a NaN multiplier must fail the check, not vanish from it.
inline void add_term(long double term, std::size_t entry, WideConstraintTerms& terms)
{
    terms.total[entry] += term;
    terms.magnitude[entry] += std::fabs(term);
    if (std::isnan(term) || std::fabs(term) > terms.scale)
    {
        terms.scale = std::fabs(term);
    }
}

// budget 1 + target_return mu + volatility Sigma w / vol + beta at the weights of the result, the budget and
// the target-return terms skipped when their multiplier is NaN, the volatility term when it is NaN or 0.
inline WideConstraintTerms wide_constraint_terms(const po::Data& data, const po::Result& result)
{
    std::size_t n = data.asset_count();
    std::vector<long double> w = widened(result.weights);
    std::vector<long double> product = wide_covariance_times(data, w);
    long double volatility = std::sqrt(wide_dot(w, product));
    long double volatility_multiplier = result.volatility_multiplier;
    bool has_volatility = !std::isnan(result.volatility_multiplier) && result.volatility_multiplier != 0.0;
    WideConstraintTerms terms{std::vector<long double>(n, 0.0L), std::vector<long double>(n, 0.0L), 0.0L};
    for (std::size_t j = 0; j < n; ++j)
    {
        add_term(result.bound_multipliers[j], j, terms);
        if (!std::isnan(result.budget_multiplier))
        {
            add_term(result.budget_multiplier, j, terms);
        }
        if (!std::isnan(result.target_return_multiplier))
        {
            add_term(static_cast<long double>(result.target_return_multiplier) * data.mu()[j], j, terms);
        }
        if (has_volatility)
        {
            add_term(volatility_multiplier * product[j] / volatility, j, terms);
        }
    }
    return terms;
}

// The stationarity residual r = grad f - budget 1 - target_return mu - volatility Sigma w / vol - beta of a
// result, evaluated in long double, by two measures.
struct RiskParityStationarity
{
    // max_j |r_j| over the largest term of grad f or of the constraints (at f = 0, over the unit-gap scale).
    double relative;
    // max_j |r_j| / floor_j, floor_j = eps ((|hess f - volatility hess vol| |w|)_j + the magnitudes of the
    // constraint terms of entry j): how far r_j moves when every weight and every multiplier of the result,
    // all doubles, is rounded once. A KKT point rounded to doubles sits at about 1.
    double over_floor;
};

inline double largest_over_floor(const std::vector<long double>& residual,
                                 const std::vector<long double>& hessian, const WideConstraintTerms& terms,
                                 std::span<const double> weights)
{
    std::size_t n = residual.size();
    double worst = 0.0;
    for (std::size_t j = 0; j < n; ++j)
    {
        long double moved = terms.magnitude[j];
        for (std::size_t l = 0; l < n; ++l)
        {
            moved += std::fabs(hessian[j * n + l] * weights[l]);
        }
        long double floor = std::numeric_limits<double>::epsilon() * moved;
        if (residual[j] != 0.0L)
        {
            worst = worst_of(worst, static_cast<double>(std::fabs(residual[j]) / floor));
        }
    }
    return worst;
}

inline RiskParityStationarity risk_parity_stationarity(const po::Spec& spec, const po::Data& data,
                                                       const po::Result& result, double tolerance)
{
    std::size_t n = data.asset_count();
    std::vector<double> budget = budget_of(spec, n);
    std::vector<long double> w = widened(result.weights);
    std::array<std::vector<long double>, 3> gradient = wide_gradient_terms_at(data, budget, w);
    WideConstraintTerms terms = wide_constraint_terms(data, result);
    std::vector<long double> residual(n);
    double scale = static_cast<double>(terms.scale);
    double largest = 0.0;
    for (std::size_t j = 0; j < n; ++j)
    {
        residual[j] = gradient[0][j] + gradient[1][j] + gradient[2][j] - terms.total[j];
        largest = worst_of(largest, static_cast<double>(std::fabs(residual[j])));
        for (const std::vector<long double>& term : gradient)
        {
            scale = worst_of(scale, static_cast<double>(std::fabs(term[j])));
        }
    }
    if (contribution_gap(data, budget, result.weights) <= tolerance)
    {
        scale = worst_of(scale, unit_gap_gradient_scale(data, result.weights));
    }
    std::vector<long double> hessian = wide_lagrangian_hessian(data, budget, result.volatility_multiplier, w);
    return RiskParityStationarity{relative_to(largest, scale),
                                  largest_over_floor(residual, hessian, terms, result.weights)};
}

// The multipliers of the equalities at the KKT point refined in long double from a result: the six Newton
// steps on its active set that forward_error takes. None where a Newton system is singular.
inline std::optional<std::vector<long double>> refined_multipliers(const po::Spec& spec, const po::Data& data,
                                                                   const po::Result& result)
{
    constexpr int kSteps = 6;
    std::vector<double> budget = budget_of(spec, data.asset_count());
    WideKkt kkt = wide_kkt_of(spec, result);
    for (int step = 0; step < kSteps; ++step)
    {
        if (!wide_newton_step(spec, data, budget, kkt))
        {
            return std::nullopt;
        }
    }
    return kkt.multipliers;
}

// The component of vector along the unit vector direction, removed from it.
inline long double project_out(const std::vector<long double>& direction, std::vector<long double>& vector)
{
    long double projection = 0.0L;
    for (std::size_t i = 0; i < vector.size(); ++i)
    {
        projection += direction[i] * vector[i];
    }
    for (std::size_t i = 0; i < vector.size(); ++i)
    {
        vector[i] -= projection * direction[i];
    }
    return projection;
}

inline long double wide_norm(const std::vector<long double>& vector)
{
    return std::sqrt(compensated_dot_of<long double>(vector.size(), vector, vector));
}

// min over x of |rhs - sum_k x_k columns_k|, by modified Gram-Schmidt in long double. None when a column
// outside the span of the earlier ones is at most 1e-12 of its norm, the rule of the QP solver for dependent
// rows: the coefficients are then not determined.
inline std::optional<std::vector<long double>> wide_least_squares(
    const std::vector<std::vector<long double>>& columns, std::vector<long double> rhs)
{
    constexpr long double kDependent = 1e-12L;
    std::size_t m = columns.size();
    std::vector<std::vector<long double>> basis;
    std::vector<long double> upper(m * m, 0.0L);
    for (std::size_t k = 0; k < m; ++k)
    {
        std::vector<long double> column = columns[k];
        long double norm = wide_norm(column);
        for (std::size_t i = 0; i < k; ++i)
        {
            upper[i * m + k] = project_out(basis[i], column);
        }
        long double outside = wide_norm(column);
        if (!(outside > kDependent * norm))
        {
            return std::nullopt;
        }
        for (long double& entry : column)
        {
            entry /= outside;
        }
        upper[k * m + k] = outside;
        basis.push_back(std::move(column));
    }
    std::vector<long double> solution(m);
    for (std::size_t k = 0; k < m; ++k)
    {
        solution[k] = project_out(basis[k], rhs);
    }
    for (std::size_t k = m; k-- > 0;)
    {
        for (std::size_t i = k + 1; i < m; ++i)
        {
            solution[k] -= upper[k * m + i] * solution[i];
        }
        solution[k] /= upper[k * m + k];
    }
    return solution;
}

inline std::vector<long double> over_free(const std::vector<long double>& values, const WideKkt& kkt)
{
    std::vector<long double> entries;
    for (std::size_t j : kkt.free)
    {
        entries.push_back(values[j]);
    }
    return entries;
}

inline void set_equality_multipliers(const po::Spec& spec, const std::vector<long double>& multipliers,
                                     po::Result& result)
{
    result.budget_multiplier = static_cast<double>(multipliers.front());
    if (spec.target_return.has_value())
    {
        result.target_return_multiplier = static_cast<double>(multipliers[1]);
    }
    if (spec.target_volatility.has_value())
    {
        result.volatility_multiplier = static_cast<double>(multipliers.back());
    }
}

// The result with the multipliers the KKT gate of the library checks, refit at its weights in long double:
// for the equalities, those that leave the smallest stationarity over the free weights (bound multiplier 0);
// for a held weight, what its row then leaves. None where the free weights do not determine them: fewer free
// weights than equalities, or their gradients over the free weights dependent.
inline std::optional<po::Result> least_squares_fit(const po::Spec& spec, const po::Data& data,
                                                   const po::Result& result)
{
    std::vector<double> budget = budget_of(spec, data.asset_count());
    WideKkt kkt = wide_kkt_of(spec, result);
    WideEqualities equalities = wide_equalities(spec, data, kkt.weights);
    std::vector<long double> residual = wide_stationarity(data, budget, kkt, equalities);
    std::vector<std::vector<long double>> columns;
    for (const std::vector<long double>& row : equalities.rows)
    {
        columns.push_back(over_free(row, kkt));
    }
    std::optional<std::vector<long double>> change = wide_least_squares(columns, over_free(residual, kkt));
    if (!change.has_value())
    {
        return std::nullopt;
    }
    po::Result fit = result;
    for (std::size_t k = 0; k < columns.size(); ++k)
    {
        kkt.multipliers[k] += (*change)[k];
        for (std::size_t j = 0; j < residual.size(); ++j)
        {
            residual[j] -= (*change)[k] * equalities.rows[k][j];
        }
    }
    set_equality_multipliers(spec, kkt.multipliers, fit);
    for (std::size_t j = 0; j < residual.size(); ++j)
    {
        if (result.bound_multipliers[j] != 0.0)
        {
            fit.bound_multipliers[j] = static_cast<double>(residual[j]);
        }
    }
    return fit;
}

// max_k |lambda_k - reference_k| / max_k |reference_k| over the multipliers of the equalities; absolute where
// every reference multiplier is 0 (no target and the budget alone active, f = 0).
inline double multiplier_distance(const std::vector<long double>& multipliers,
                                  const std::vector<long double>& reference)
{
    long double scale = 0.0L;
    long double distance = 0.0L;
    for (std::size_t k = 0; k < reference.size(); ++k)
    {
        scale = std::fmax(scale, std::fabs(reference[k]));
        distance = std::fmax(distance, std::fabs(multipliers[k] - reference[k]));
    }
    if (scale == 0.0L)
    {
        return static_cast<double>(distance);
    }
    return static_cast<double>(distance / scale);
}

// The multipliers of the equalities of the result are no farther from those of the KKT point refined from it
// than the least-squares fit is, but for the tolerance: a correction of the least-squares multipliers may
// only bring them closer. Not checked where the refinement meets a singular system.
inline void expect_multipliers_toward_the_reference(const po::Spec& spec, const po::Data& data,
                                                    const po::Result& result, const po::Result& fit,
                                                    double tolerance)
{
    std::optional<std::vector<long double>> reference = refined_multipliers(spec, data, result);
    if (!reference.has_value())
    {
        return;
    }
    double returned = multiplier_distance(wide_kkt_of(spec, result).multipliers, *reference);
    double fitted = multiplier_distance(wide_kkt_of(spec, fit).multipliers, *reference);
    EXPECT_LE(returned, fitted + tolerance)
        << "multipliers: " << returned << " from those of the refined KKT point, above the " << fitted
        << " of the least-squares fit plus the tolerance " << tolerance;
}

// A stationarity above the tolerance passes when every entry of the residual is within this many times its
// rounding floor: no doubles near the solution do better (docs/guide/validation.md).
inline constexpr double kRoundingFloorFactor = 2.0;

// KKT in the user's variables: equalities, bounds, sign and complementarity of beta, the NaN pattern of the
// multipliers, and
//   grad f(w) = budget * 1 + target_return * mu + volatility * Sigma w / vol + beta,
// with grad f from risk_parity_gradient, and the identity w'grad f = 0. The stationarity passes at the
// tolerance, relative to its largest term, or within kRoundingFloorFactor of its rounding floor. It is
// checked with the multipliers of least_squares_fit where the free weights determine them, as the KKT gate of
// the library does, and with those of the result elsewhere: at weights off the KKT point by more than their
// rounding, even the exact multipliers leave a stationarity of H dw, above the floor. The multipliers of the
// result are then checked against the reference instead (expect_multipliers_toward_the_reference). A
// volatility multiplier of -inf is the case s = sigma_min, which has no stationarity to check
// (expect_minimum_volatility_point). Where long double is not extended (extended false) a stationarity above
// the tolerance cannot be checked against its floor, and the test is skipped instead; nor is there a
// reference for the multipliers.
inline void expect_risk_parity_kkt(const po::Spec& spec, const po::Data& data, const po::Result& result,
                                   double tolerance, bool extended = kLongDoubleIsExtended)
{
    ASSERT_EQ(result.status, po::SolveStatus::Optimal);
    ASSERT_EQ(result.weights.size(), data.weight_count());
    ASSERT_EQ(result.bound_multipliers.size(), data.weight_count());
    if (result.volatility_multiplier == -std::numeric_limits<double>::infinity())
    {
        expect_minimum_volatility_point(spec, data, result, tolerance);
        return;
    }
    expect_risk_parity_multiplier_pattern(spec, result);
    expect_risk_parity_feasible(spec, data, result, tolerance);
    expect_user_bounds(spec, result, tolerance);
    std::vector<double> budget = budget_of(spec, data.asset_count());
    ObjectiveGradient gradient = sum_of(risk_parity_gradient_terms(data, budget, result.weights));
    if (contribution_gap(data, budget, result.weights) <= tolerance)
    {
        gradient.scale = worst_of(gradient.scale, unit_gap_gradient_scale(data, result.weights));
    }
    expect_homogeneity_identity(result.weights, gradient, tolerance);
    std::optional<po::Result> fit = least_squares_fit(spec, data, result);
    const po::Result* checked = &result;
    if (fit.has_value())
    {
        checked = &*fit;
    }
    RiskParityStationarity stationarity = risk_parity_stationarity(spec, data, *checked, tolerance);
    if (stationarity.relative > tolerance && !extended)
    {
        // The floor is eps of the magnitudes of the terms; in double their sum is rounded by about as much.
        GTEST_SKIP() << "stationarity " << stationarity.relative << " above the tolerance " << tolerance
                     << ": its rounding floor cannot be checked where long double is double";
    }
    EXPECT_TRUE(stationarity.relative <= tolerance || stationarity.over_floor <= kRoundingFloorFactor)
        << "stationarity: " << stationarity.relative << " of the largest term, above the tolerance "
        << tolerance << ", and " << stationarity.over_floor << " of its rounding floor, above "
        << kRoundingFloorFactor;
    if (fit.has_value() && extended)
    {
        expect_multipliers_toward_the_reference(spec, data, result, *fit, tolerance);
    }
}

// The message of every failure of expect_risk_parity_kkt on the result; empty when it passes. The failures
// are intercepted, not reported. A skip is not a failure: it is passed on to the running test, which is then
// reported as skipped (unless something else fails).
inline std::vector<std::string> risk_parity_kkt_failures(const po::Spec& spec, const po::Data& data,
                                                         const po::Result& result, double tolerance,
                                                         bool extended = kLongDoubleIsExtended)
{
    ::testing::TestPartResultArray failures;
    {
        ::testing::ScopedFakeTestPartResultReporter reporter(
            ::testing::ScopedFakeTestPartResultReporter::INTERCEPT_ONLY_CURRENT_THREAD, &failures);
        expect_risk_parity_kkt(spec, data, result, tolerance, extended);
    }
    std::vector<std::string> messages;
    for (int k = 0; k < failures.size(); ++k)
    {
        const ::testing::TestPartResult& part = failures.GetTestPartResult(k);
        if (part.skipped())
        {
            [&part] { GTEST_SKIP() << part.message(); }();
            continue;
        }
        messages.emplace_back(part.message());
    }
    return messages;
}

} // namespace po::testing
