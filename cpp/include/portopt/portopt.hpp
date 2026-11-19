// PortOpt for C++23: a header-only wrapper over the C interface (portopt/portopt.h). It translates and nothing else:
// a solve through it is the C solve bit for bit. Errors are std::expected values carrying the status and the message
// of po_last_error(); handles are freed by RAII. The names follow the Python package.
#ifndef PORTOPT_CPP_PORTOPT_HPP
#define PORTOPT_CPP_PORTOPT_HPP

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <memory>
#include <optional>
#include <portopt/portopt.h>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace po
{

enum class Status : std::int32_t
{
    ok = PO_OK,
    invalid_argument = PO_INVALID_ARGUMENT,
    invalid_handle = PO_INVALID_HANDLE,
    unsupported = PO_UNSUPPORTED,
    out_of_memory = PO_OUT_OF_MEMORY,
    internal_error = PO_INTERNAL_ERROR,
};

// A failed call: its status and the message of po_last_error(), copied when the call returned.
struct Error
{
    Status status;
    std::string message;
};

template <class T>
using Expected = std::expected<T, Error>;

enum class Objective : std::int32_t
{
    min_variance = PO_OBJECTIVE_MIN_VARIANCE,
    max_return = PO_OBJECTIVE_MAX_RETURN,
    max_sharpe = PO_OBJECTIVE_MAX_SHARPE,
    utility = PO_OBJECTIVE_UTILITY,
    importance = PO_OBJECTIVE_IMPORTANCE,
    kelly = PO_OBJECTIVE_KELLY,
    tracking_error = PO_OBJECTIVE_TRACKING_ERROR,
    max_diversification = PO_OBJECTIVE_MAX_DIVERSIFICATION,
    risk_parity = PO_OBJECTIVE_RISK_PARITY,
};

enum class SolveStatus : std::int32_t
{
    optimal = PO_SOLVE_OPTIMAL,
    infeasible = PO_SOLVE_INFEASIBLE,
    unbounded = PO_SOLVE_UNBOUNDED,
    max_iter = PO_SOLVE_MAX_ITER,
};

// One pair for every weight, cash included (po_spec_set_bounds).
struct Bounds
{
    double lower;
    double upper;
};

// One pair per weight, cash last (po_spec_set_asset_bounds).
struct AssetBounds
{
    std::vector<double> lower;
    std::vector<double> upper;
};

// What to solve, as values: each field that is set is passed to its po_spec_set_* setter, whose header comment
// gives its domain and the objectives it applies to. An empty risk_budget is the default 1 / n.
struct Spec
{
    Objective objective = Objective::min_variance;
    std::variant<std::monostate, Bounds, AssetBounds> bounds{};
    std::optional<double> target_return{};
    std::optional<double> risk_aversion{};
    std::optional<double> importance{};
    std::optional<double> max_volatility{};
    std::optional<double> leverage{};
    std::optional<double> target_volatility{};
    std::vector<double> risk_budget{};
};

// The result of a solve, read out of the native result; the meaning of each field is in portopt.h (po_result_*).
struct Result
{
    SolveStatus status;
    std::vector<double> weights;
    std::int32_t iterations;
    double budget_multiplier;
    double target_return_multiplier;
    double volatility_multiplier;
    double leverage_multiplier;
    std::vector<double> bound_multipliers;
    double residual;
};

// The version of the library, "major.minor.patch".
inline std::string_view version()
{
    return po_version();
}

namespace detail
{

// The message of po_last_error() is valid only until the next call on this thread, so it is copied here.
inline Expected<void> check(po_status status)
{
    if (status == PO_OK)
    {
        return {};
    }
    const char* message = po_last_error();
    return std::unexpected(Error{static_cast<Status>(status), message != nullptr ? message : ""});
}

inline std::unexpected<Error> invalid(std::string message)
{
    return std::unexpected(Error{Status::invalid_argument, std::move(message)});
}

inline Expected<void> check_abi_version(std::int32_t native, std::int32_t expected = PO_ABI_VERSION)
{
    if (native == expected)
    {
        return {};
    }
    return std::unexpected(Error{Status::unsupported, "the PortOpt library has ABI version " +
                                                          std::to_string(native) + ", but these headers need version " +
                                                          std::to_string(expected)});
}

// A span the C side reads `expected` entries from: an error unless it has exactly that many.
inline Expected<void> check_size(std::string_view name, std::size_t size, std::size_t expected, std::string what)
{
    if (size == expected)
    {
        return {};
    }
    return invalid(std::string(name) + " has " + std::to_string(size) + " entries, expected " + what);
}

inline const double* data_or_null(std::span<const double> values)
{
    if (values.empty())
    {
        return nullptr;
    }
    return values.data();
}

template <class Handle, void (*Free)(Handle*)>
struct Deleter
{
    void operator()(Handle* handle) const noexcept
    {
        Free(handle);
    }
};

using SpecHandle = std::unique_ptr<po_spec, Deleter<po_spec, po_spec_free>>;
using DataHandle = std::unique_ptr<po_data, Deleter<po_data, po_data_free>>;
using ResultHandle = std::unique_ptr<po_result, Deleter<po_result, po_result_free>>;

inline Expected<void> set_bounds(po_spec* spec, const std::variant<std::monostate, Bounds, AssetBounds>& bounds)
{
    if (const auto* scalar = std::get_if<Bounds>(&bounds))
    {
        return check(po_spec_set_bounds(spec, scalar->lower, scalar->upper));
    }
    const auto* per_weight = std::get_if<AssetBounds>(&bounds);
    if (per_weight == nullptr)
    {
        return {};
    }
    if (per_weight->lower.size() != per_weight->upper.size())
    {
        return invalid("the lower bounds have " + std::to_string(per_weight->lower.size()) +
                       " entries and the upper bounds " + std::to_string(per_weight->upper.size()));
    }
    return check(po_spec_set_asset_bounds(spec, per_weight->lower.data(), per_weight->upper.data(),
                                          per_weight->lower.size()));
}

struct ScalarSetter
{
    const std::optional<double>* value;
    po_status (*set)(po_spec*, double);
};

// The native spec of a Spec, freed when the handle goes out of scope.
inline Expected<SpecHandle> make_spec(const Spec& spec)
{
    po_spec* raw = nullptr;
    if (auto created = check(po_spec_create(static_cast<std::int32_t>(spec.objective), &raw)); !created)
    {
        return std::unexpected(std::move(created.error()));
    }
    SpecHandle handle(raw);
    if (auto bounds = set_bounds(handle.get(), spec.bounds); !bounds)
    {
        return std::unexpected(std::move(bounds.error()));
    }
    const std::array<ScalarSetter, 6> setters{{
        {&spec.target_return, po_spec_set_target_return},
        {&spec.risk_aversion, po_spec_set_risk_aversion},
        {&spec.importance, po_spec_set_importance},
        {&spec.max_volatility, po_spec_set_max_volatility},
        {&spec.leverage, po_spec_set_leverage},
        {&spec.target_volatility, po_spec_set_target_volatility},
    }};
    for (const ScalarSetter& setter : setters)
    {
        if (!setter.value->has_value())
        {
            continue;
        }
        if (auto set = check(setter.set(handle.get(), **setter.value)); !set)
        {
            return std::unexpected(std::move(set.error()));
        }
    }
    if (spec.risk_budget.empty())
    {
        return handle;
    }
    if (auto budget = check(po_spec_set_risk_budget(handle.get(), spec.risk_budget.data(), spec.risk_budget.size()));
        !budget)
    {
        return std::unexpected(std::move(budget.error()));
    }
    return handle;
}

inline Expected<Result> read_result(const po_result* handle)
{
    Result result{};
    std::size_t count = 0;
    std::int32_t status = 0;
    const auto multiplier = [handle](std::int32_t constraint, double& out) {
        return check(po_result_multiplier(handle, constraint, &out));
    };
    auto read = check(po_result_weight_count(handle, &count))
                    .and_then([&] {
                        result.weights.resize(count);
                        result.bound_multipliers.resize(count);
                        return check(po_result_status(handle, &status));
                    })
                    .and_then([&] { return check(po_result_weights(handle, result.weights.data(), count)); })
                    .and_then([&] { return check(po_result_iterations(handle, &result.iterations)); })
                    .and_then([&] { return multiplier(PO_CONSTRAINT_BUDGET, result.budget_multiplier); })
                    .and_then([&] { return multiplier(PO_CONSTRAINT_TARGET_RETURN, result.target_return_multiplier); })
                    .and_then([&] { return multiplier(PO_CONSTRAINT_VOLATILITY, result.volatility_multiplier); })
                    .and_then([&] { return multiplier(PO_CONSTRAINT_LEVERAGE, result.leverage_multiplier); })
                    .and_then([&] {
                        return check(po_result_bound_multipliers(handle, result.bound_multipliers.data(), count));
                    })
                    .and_then([&] { return check(po_result_residual(handle, &result.residual)); });
    if (!read)
    {
        return std::unexpected(std::move(read.error()));
    }
    result.status = static_cast<SolveStatus>(status);
    return result;
}

} // namespace detail

// Market data: expected returns, covariance, an optional risk-free asset and return series. Move-only; a moved-from
// Data holds no handle, and every call on it is Status::invalid_handle.
class Data
{
public:
    // mu has `assets` entries and cov `assets x assets` (row-major); an empty span is a missing input
    // (po_data_create). The ABI version of the library is checked here.
    static Expected<Data> create(std::size_t assets, std::span<const double> mu, std::span<const double> cov)
    {
        if (auto abi = detail::check_abi_version(po_abi_version()); !abi)
        {
            return std::unexpected(std::move(abi.error()));
        }
        if (!mu.empty())
        {
            if (auto size = detail::check_size("mu", mu.size(), assets, std::to_string(assets)); !size)
            {
                return std::unexpected(std::move(size.error()));
            }
        }
        if (!cov.empty())
        {
            const std::string what =
                std::to_string(assets) + " x " + std::to_string(assets) + " = " + std::to_string(assets * assets);
            if (auto size = detail::check_size("cov", cov.size(), assets * assets, what); !size)
            {
                return std::unexpected(std::move(size.error()));
            }
        }
        po_data* raw = nullptr;
        if (auto created = detail::check(po_data_create(assets, detail::data_or_null(mu), detail::data_or_null(cov),
                                                        &raw));
            !created)
        {
            return std::unexpected(std::move(created.error()));
        }
        return Data(raw);
    }

    Expected<void> set_risk_free(double rate)
    {
        return detail::check(po_data_set_risk_free(handle_.get(), rate));
    }

    // returns is periods x assets (row-major) and benchmark has one entry per period (po_data_set_returns).
    Expected<void> set_returns(std::span<const double> returns, std::span<const double> benchmark)
    {
        Expected<std::size_t> assets = asset_count();
        if (!assets)
        {
            return std::unexpected(std::move(assets.error()));
        }
        const std::size_t periods = benchmark.size();
        const std::string what = std::to_string(periods) + " x " + std::to_string(*assets) + " = " +
                                 std::to_string(periods * *assets);
        if (auto size = detail::check_size("returns", returns.size(), periods * *assets, what); !size)
        {
            return std::unexpected(std::move(size.error()));
        }
        return detail::check(po_data_set_returns(handle_.get(), detail::data_or_null(returns),
                                                 detail::data_or_null(benchmark), periods));
    }

    Expected<std::size_t> asset_count() const
    {
        return count(po_data_asset_count);
    }

    // The asset count, plus one with a risk-free asset.
    Expected<std::size_t> weight_count() const
    {
        return count(po_data_weight_count);
    }

    // The metrics of a weight vector of weight_count() entries, cash last (po_portfolio_return and the rest).
    Expected<double> portfolio_return(std::span<const double> weights) const
    {
        return metric(po_portfolio_return, weights);
    }

    Expected<double> volatility(std::span<const double> weights) const
    {
        return metric(po_portfolio_volatility, weights);
    }

    Expected<double> utility(std::span<const double> weights, double risk_aversion) const
    {
        double value = 0.0;
        if (auto read = detail::check(po_utility(handle_.get(), weights.data(), weights.size(), risk_aversion, &value));
            !read)
        {
            return std::unexpected(std::move(read.error()));
        }
        return value;
    }

    Expected<double> tracking_error(std::span<const double> weights) const
    {
        return metric(po_tracking_error, weights);
    }

    Expected<double> sharpe_ratio(std::span<const double> weights) const
    {
        return metric(po_sharpe_ratio, weights);
    }

    Expected<double> diversification_ratio(std::span<const double> weights) const
    {
        return metric(po_diversification_ratio, weights);
    }

    Expected<std::vector<double>> risk_contributions(std::span<const double> weights) const
    {
        std::vector<double> contributions(weights.size());
        if (auto read = detail::check(
                po_risk_contributions(handle_.get(), weights.data(), weights.size(), contributions.data()));
            !read)
        {
            return std::unexpected(std::move(read.error()));
        }
        return contributions;
    }

    const po_data* get() const noexcept
    {
        return handle_.get();
    }

private:
    using Metric = po_status (*)(const po_data*, const double*, std::size_t, double*);

    explicit Data(po_data* handle) noexcept : handle_(handle)
    {
    }

    Expected<std::size_t> count(po_status (*read)(const po_data*, std::size_t*)) const
    {
        std::size_t value = 0;
        if (auto counted = detail::check(read(handle_.get(), &value)); !counted)
        {
            return std::unexpected(std::move(counted.error()));
        }
        return value;
    }

    // The C side checks that the weight count matches, so a short span is an error there, not a read past it.
    Expected<double> metric(Metric function, std::span<const double> weights) const
    {
        double value = 0.0;
        if (auto read = detail::check(function(handle_.get(), weights.data(), weights.size(), &value)); !read)
        {
            return std::unexpected(std::move(read.error()));
        }
        return value;
    }

    detail::DataHandle handle_;
};

// Solves the spec on the data. An infeasible or unbounded problem is a Result with that status, not an error.
inline Expected<Result> solve(const Spec& spec, const Data& data)
{
    Expected<detail::SpecHandle> native = detail::make_spec(spec);
    if (!native)
    {
        return std::unexpected(std::move(native.error()));
    }
    po_result* raw = nullptr;
    if (auto solved = detail::check(po_solve(native->get(), data.get(), &raw)); !solved)
    {
        return std::unexpected(std::move(solved.error()));
    }
    const detail::ResultHandle result(raw);
    return detail::read_result(result.get());
}

} // namespace po

#endif
