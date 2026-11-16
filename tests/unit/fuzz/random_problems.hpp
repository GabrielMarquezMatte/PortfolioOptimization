#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "core/data.hpp"
#include "core/spec.hpp"
#include "random_covariance.hpp"

namespace po::testing
{

struct ProblemShape
{
    std::size_t n;
    // cond(Sigma)
    double condition;
    // 1 (annual) or 1e-4 (daily) on Sigma; its square root on mu.
    double scale;
    enum class Means
    {
        Spread,
        Equal,
        NearlyEqual,
        BelowRiskFree
    } means;
    // Tight: upper = 2 / n; Fixed: some lower == upper.
    enum class Box
    {
        None,
        LongOnly,
        Tight,
        Fixed
    } box;
    bool cash;
    std::uint64_t seed;
};

// spec holds the bounds only; the objective and its setters are the caller's.
struct GeneratedProblem
{
    po::Spec spec;
    po::Data data;
};

// The rate of the risk-free asset at the annual scale: below every mean, or above every mean.
constexpr double kShapeRiskFree = 0.02;
constexpr double kShapeRiskFreeAboveMeans = 0.2;

// NearlyEqual: mu_i = m (1 + spread u_i), u_i uniform in [-1, 1], spread log-uniform in [1e-8, 1e-5].
inline std::vector<double> nearly_equal_mu(double level, std::size_t n, std::uint64_t seed)
{
    double spread = std::pow(10.0, random_vector(1, seed, -8.0, -5.0)[0]);
    std::vector<double> mu = random_vector(n, seed + 1, -1.0, 1.0);
    for (double& value : mu)
    {
        value = level * (1.0 + spread * value);
    }
    return mu;
}

// Means in [0.03, 0.15] at the annual scale.
inline std::vector<double> spread_mu(const ProblemShape& shape)
{
    double root = std::sqrt(shape.scale);
    return random_vector(shape.n, shape.seed + 500, 0.03 * root, 0.15 * root);
}

inline std::vector<double> negated(std::vector<double> values)
{
    for (double& value : values)
    {
        value = -value;
    }
    return values;
}

// Without cash rf = 0, so BelowRiskFree negates the means.
inline std::vector<double> shape_mu(const ProblemShape& shape)
{
    if (shape.means == ProblemShape::Means::Equal)
    {
        return std::vector<double>(shape.n, spread_mu(shape).front());
    }
    if (shape.means == ProblemShape::Means::NearlyEqual)
    {
        return nearly_equal_mu(spread_mu(shape).front(), shape.n, shape.seed + 600);
    }
    if (shape.means == ProblemShape::Means::BelowRiskFree && !shape.cash)
    {
        return negated(spread_mu(shape));
    }
    return spread_mu(shape);
}

inline po::Data shape_data(const ProblemShape& shape)
{
    std::vector<double> covariance = random_covariance(shape.n, shape.seed, shape.condition);
    for (double& entry : covariance)
    {
        entry *= shape.scale;
    }
    po::Data data = po::Data::create(shape.n, shape_mu(shape), covariance).value();
    if (!shape.cash)
    {
        return data;
    }
    double rate = kShapeRiskFree;
    if (shape.means == ProblemShape::Means::BelowRiskFree)
    {
        rate = kShapeRiskFreeAboveMeans;
    }
    data.set_risk_free(rate * std::sqrt(shape.scale));
    return data;
}

// Long-only with a third of the risky weights (at least one) fixed at 0.5 / n each.
inline po::Spec fixed_box(const ProblemShape& shape, std::size_t weight_count)
{
    po::Spec spec;
    spec.asset_lower.assign(weight_count, 0.0);
    spec.asset_upper.assign(weight_count, 1.0);
    std::size_t fixed = std::max<std::size_t>(1, shape.n / 3);
    double value = 0.5 / static_cast<double>(shape.n);
    for (std::size_t i = 0; i < fixed; ++i)
    {
        spec.asset_lower[i] = value;
        spec.asset_upper[i] = value;
    }
    return spec;
}

inline po::Spec bounded_box(double upper)
{
    po::Spec spec;
    spec.lower = 0.0;
    spec.upper = upper;
    return spec;
}

inline po::Spec shape_box(const ProblemShape& shape, std::size_t weight_count)
{
    if (shape.box == ProblemShape::Box::Fixed)
    {
        return fixed_box(shape, weight_count);
    }
    if (shape.box == ProblemShape::Box::None)
    {
        return po::Spec{};
    }
    if (shape.box == ProblemShape::Box::Tight)
    {
        return bounded_box(2.0 / static_cast<double>(shape.n));
    }
    return bounded_box(1.0);
}

inline GeneratedProblem generate(const ProblemShape& shape)
{
    po::Data data = shape_data(shape);
    return GeneratedProblem{shape_box(shape, data.weight_count()), data};
}

// Shape k takes combination k mod 256 of cond(Sigma) in {1e2, 1e4, 1e6, 1e8}, the two scales, the four mean
// patterns, the four boxes and cash, so any 256 consecutive shapes cover every combination; n (3 to 24) and the
// seed are drawn from seed.
inline std::vector<ProblemShape> fuzz_shapes(std::size_t count, std::uint64_t seed)
{
    constexpr std::array<double, 4> kConditions = {1e2, 1e4, 1e6, 1e8};
    constexpr std::array<double, 2> kScales = {1.0, 1e-4};
    std::mt19937_64 engine(seed);
    std::vector<ProblemShape> shapes;
    shapes.reserve(count);
    for (std::size_t k = 0; k < count; ++k)
    {
        ProblemShape shape{};
        shape.condition = kConditions[k % 4];
        shape.scale = kScales[(k / 4) % 2];
        shape.means = static_cast<ProblemShape::Means>((k / 8) % 4);
        shape.box = static_cast<ProblemShape::Box>((k / 32) % 4);
        shape.cash = (k / 128) % 2 == 1;
        shape.n = 3 + static_cast<std::size_t>(engine() % 22);
        shape.seed = engine();
        shapes.push_back(shape);
    }
    return shapes;
}

} // namespace po::testing
