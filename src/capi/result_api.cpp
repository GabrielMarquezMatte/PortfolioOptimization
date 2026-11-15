#include <algorithm>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "capi/error.hpp"
#include "capi/handles.hpp"

using po::capi::guard;
using po::capi::null_argument;
using po::capi::null_handle;
using po::capi::set_error;

static_assert(PO_SOLVE_OPTIMAL == static_cast<int32_t>(po::SolveStatus::Optimal));
static_assert(PO_SOLVE_INFEASIBLE == static_cast<int32_t>(po::SolveStatus::Infeasible));
static_assert(PO_SOLVE_UNBOUNDED == static_cast<int32_t>(po::SolveStatus::Unbounded));
static_assert(PO_SOLVE_MAX_ITER == static_cast<int32_t>(po::SolveStatus::MaxIterations));

namespace
{

po_status copy_vector(const std::vector<double>& values, double* out, std::size_t n)
{
    if (out == nullptr)
    {
        return null_argument("out");
    }
    if (n != values.size())
    {
        return set_error(PO_INVALID_ARGUMENT, "n must be " + std::to_string(values.size()) + ", got " +
                                                  std::to_string(n));
    }
    std::ranges::copy(values, out);
    return PO_OK;
}

po_status multiplier_of(const po::Result& result, int32_t constraint, double& value)
{
    switch (constraint)
    {
    case PO_CONSTRAINT_BUDGET:
        value = result.budget_multiplier;
        return PO_OK;
    case PO_CONSTRAINT_TARGET_RETURN:
        value = result.target_return_multiplier;
        return PO_OK;
    case PO_CONSTRAINT_VOLATILITY:
        value = result.volatility_multiplier;
        return PO_OK;
    case PO_CONSTRAINT_LEVERAGE:
        value = result.leverage_multiplier;
        return PO_OK;
    default:
        return set_error(PO_INVALID_ARGUMENT, "constraint is not one of PO_CONSTRAINT_*");
    }
}

template <typename Read>
po_status read_scalar(const po_result* r, Read&& read)
{
    return guard(
        [&]
        {
            if (r == nullptr)
            {
                return null_handle("r");
            }
            return read(r->result);
        });
}

} // namespace

extern "C" po_status po_result_status(const po_result* r, int32_t* out)
{
    return read_scalar(r,
                       [&](const po::Result& result)
                       {
                           if (out == nullptr)
                           {
                               return null_argument("out");
                           }
                           *out = static_cast<int32_t>(result.status);
                           return PO_OK;
                       });
}

extern "C" po_status po_result_weight_count(const po_result* r, size_t* out)
{
    return read_scalar(r,
                       [&](const po::Result& result)
                       {
                           if (out == nullptr)
                           {
                               return null_argument("out");
                           }
                           *out = result.weights.size();
                           return PO_OK;
                       });
}

extern "C" po_status po_result_weights(const po_result* r, double* out, size_t n)
{
    return read_scalar(r, [&](const po::Result& result) { return copy_vector(result.weights, out, n); });
}

extern "C" po_status po_result_iterations(const po_result* r, int32_t* out)
{
    return read_scalar(r,
                       [&](const po::Result& result)
                       {
                           if (out == nullptr)
                           {
                               return null_argument("out");
                           }
                           *out = result.iterations;
                           return PO_OK;
                       });
}

extern "C" po_status po_result_multiplier(const po_result* r, int32_t constraint, double* out)
{
    return read_scalar(r,
                       [&](const po::Result& result)
                       {
                           if (out == nullptr)
                           {
                               return null_argument("out");
                           }
                           return multiplier_of(result, constraint, *out);
                       });
}

extern "C" po_status po_result_bound_multipliers(const po_result* r, double* out, size_t n)
{
    return read_scalar(r, [&](const po::Result& result) { return copy_vector(result.bound_multipliers, out, n); });
}

extern "C" po_status po_result_residual(const po_result* r, double* out)
{
    return read_scalar(r,
                       [&](const po::Result& result)
                       {
                           if (out == nullptr)
                           {
                               return null_argument("out");
                           }
                           *out = result.residual;
                           return PO_OK;
                       });
}

extern "C" void po_result_free(po_result* r)
{
    std::unique_ptr<po_result> owner(r);
}
