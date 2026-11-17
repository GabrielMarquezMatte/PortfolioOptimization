#pragma once

#include <portopt/portopt.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>

#include "core/data.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"

namespace po::bench
{

struct SpecDeleter
{
    void operator()(po_spec* spec) const noexcept
    {
        po_spec_free(spec);
    }
};

struct DataDeleter
{
    void operator()(po_data* data) const noexcept
    {
        po_data_free(data);
    }
};

struct ResultDeleter
{
    void operator()(po_result* result) const noexcept
    {
        po_result_free(result);
    }
};

using SpecHandle = std::unique_ptr<po_spec, SpecDeleter>;
using DataHandle = std::unique_ptr<po_data, DataDeleter>;
using ResultHandle = std::unique_ptr<po_result, ResultDeleter>;

inline void require_ok(po_status status, const char* call)
{
    if (status != PO_OK)
    {
        throw std::runtime_error(std::string(call) + ": " + po_last_error());
    }
}

inline void set_scalar(po_status (*setter)(po_spec*, double), po_spec* spec, const std::optional<double>& value,
                       const char* name)
{
    if (value.has_value())
    {
        require_ok(setter(spec, *value), name);
    }
}

// The C handle of a spec of the benchmark cases (scalar bounds and scalar options; no asset bounds or budget).
inline SpecHandle make_c_spec(const po::Spec& spec)
{
    po_spec* raw = nullptr;
    require_ok(po_spec_create(static_cast<std::int32_t>(spec.objective), &raw), "po_spec_create");
    SpecHandle handle(raw);
    require_ok(po_spec_set_bounds(raw, spec.lower, spec.upper), "po_spec_set_bounds");
    set_scalar(po_spec_set_target_return, raw, spec.target_return, "po_spec_set_target_return");
    set_scalar(po_spec_set_risk_aversion, raw, spec.risk_aversion, "po_spec_set_risk_aversion");
    set_scalar(po_spec_set_importance, raw, spec.importance, "po_spec_set_importance");
    set_scalar(po_spec_set_max_volatility, raw, spec.max_volatility, "po_spec_set_max_volatility");
    set_scalar(po_spec_set_leverage, raw, spec.leverage, "po_spec_set_leverage");
    set_scalar(po_spec_set_target_volatility, raw, spec.target_volatility, "po_spec_set_target_volatility");
    return handle;
}

inline DataHandle make_c_data(const po::Data& data)
{
    po_data* raw = nullptr;
    require_ok(po_data_create(data.asset_count(), data.mu().data(), data.covariance().data(), &raw),
               "po_data_create");
    DataHandle handle(raw);
    if (data.risk_free().has_value())
    {
        require_ok(po_data_set_risk_free(raw, *data.risk_free()), "po_data_set_risk_free");
    }
    if (data.has_returns())
    {
        require_ok(po_data_set_returns(raw, data.returns().data(), data.benchmark().data(),
                                       data.observation_count()),
                   "po_data_set_returns");
    }
    return handle;
}

inline double read_multiplier(const po_result* result, std::int32_t constraint)
{
    double value = 0.0;
    require_ok(po_result_multiplier(result, constraint, &value), "po_result_multiplier");
    return value;
}

// Every field of the result, through the readers of the C interface.
inline po::Result read_c_result(const po_result* result)
{
    po::Result out;
    std::int32_t status = 0;
    require_ok(po_result_status(result, &status), "po_result_status");
    out.status = static_cast<po::SolveStatus>(status);
    std::size_t count = 0;
    require_ok(po_result_weight_count(result, &count), "po_result_weight_count");
    out.weights.resize(count);
    require_ok(po_result_weights(result, out.weights.data(), count), "po_result_weights");
    std::int32_t iterations = 0;
    require_ok(po_result_iterations(result, &iterations), "po_result_iterations");
    out.iterations = iterations;
    out.budget_multiplier = read_multiplier(result, PO_CONSTRAINT_BUDGET);
    out.target_return_multiplier = read_multiplier(result, PO_CONSTRAINT_TARGET_RETURN);
    out.volatility_multiplier = read_multiplier(result, PO_CONSTRAINT_VOLATILITY);
    out.leverage_multiplier = read_multiplier(result, PO_CONSTRAINT_LEVERAGE);
    out.bound_multipliers.resize(count);
    require_ok(po_result_bound_multipliers(result, out.bound_multipliers.data(), count),
               "po_result_bound_multipliers");
    require_ok(po_result_residual(result, &out.residual), "po_result_residual");
    return out;
}

inline ResultHandle solve_with_handles(const SpecHandle& spec, const DataHandle& data)
{
    po_result* raw = nullptr;
    require_ok(po_solve(spec.get(), data.get(), &raw), "po_solve");
    return ResultHandle(raw);
}

// The whole cycle of a C caller: handles from the C++ spec and data, po_solve, every reader, every free.
inline po::Result solve_through_c_api(const po::Spec& spec, const po::Data& data)
{
    SpecHandle c_spec = make_c_spec(spec);
    DataHandle c_data = make_c_data(data);
    ResultHandle c_result = solve_with_handles(c_spec, c_data);
    return read_c_result(c_result.get());
}

} // namespace po::bench
