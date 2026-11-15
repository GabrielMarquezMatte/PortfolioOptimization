#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "capi/error.hpp"
#include "capi/handles.hpp"

using po::capi::guard;
using po::capi::null_argument;
using po::capi::null_handle;
using po::capi::set_error;

static_assert(PO_OBJECTIVE_MIN_VARIANCE == static_cast<int32_t>(po::Objective::MinVariance));
static_assert(PO_OBJECTIVE_MAX_RETURN == static_cast<int32_t>(po::Objective::MaxReturn));
static_assert(PO_OBJECTIVE_MAX_SHARPE == static_cast<int32_t>(po::Objective::MaxSharpe));
static_assert(PO_OBJECTIVE_UTILITY == static_cast<int32_t>(po::Objective::Utility));
static_assert(PO_OBJECTIVE_IMPORTANCE == static_cast<int32_t>(po::Objective::Importance));
static_assert(PO_OBJECTIVE_KELLY == static_cast<int32_t>(po::Objective::Kelly));
static_assert(PO_OBJECTIVE_TRACKING_ERROR == static_cast<int32_t>(po::Objective::TrackingError));
static_assert(PO_OBJECTIVE_MAX_DIVERSIFICATION == static_cast<int32_t>(po::Objective::MaxDiversification));
static_assert(PO_OBJECTIVE_RISK_PARITY == static_cast<int32_t>(po::Objective::RiskParity));

namespace
{

constexpr double kInfinity = std::numeric_limits<double>::infinity();

bool valid_bounds(double lower, double upper)
{
    if (std::isnan(lower) || std::isnan(upper))
    {
        return false;
    }
    return lower <= upper && lower != kInfinity && upper != -kInfinity;
}

po_status invalid_bounds()
{
    return set_error(PO_INVALID_ARGUMENT,
                     "bounds must not be NaN and must satisfy lower <= upper, lower < +inf and upper > -inf");
}

bool all_valid(const double* lower, const double* upper, std::size_t n)
{
    for (std::size_t i = 0; i < n; ++i)
    {
        if (!valid_bounds(lower[i], upper[i]))
        {
            return false;
        }
    }
    return true;
}

} // namespace

extern "C" po_status po_spec_create(int32_t objective, po_spec** out)
{
    if (out != nullptr)
    {
        *out = nullptr;
    }
    return guard(
        [&]
        {
            if (out == nullptr)
            {
                return null_argument("out");
            }
            if (objective < PO_OBJECTIVE_MIN_VARIANCE || objective > PO_OBJECTIVE_RISK_PARITY)
            {
                return set_error(PO_INVALID_ARGUMENT, "objective is not one of PO_OBJECTIVE_*");
            }
            auto handle = std::make_unique<po_spec>();
            handle->spec.objective = static_cast<po::Objective>(objective);
            *out = handle.release();
            return PO_OK;
        });
}

extern "C" po_status po_spec_set_bounds(po_spec* s, double lower, double upper)
{
    return guard(
        [&]
        {
            if (s == nullptr)
            {
                return null_handle("s");
            }
            if (!valid_bounds(lower, upper))
            {
                return invalid_bounds();
            }
            s->spec.lower = lower;
            s->spec.upper = upper;
            s->spec.asset_lower.clear();
            s->spec.asset_upper.clear();
            return PO_OK;
        });
}

extern "C" po_status po_spec_set_asset_bounds(po_spec* s, const double* lower, const double* upper, size_t n)
{
    return guard(
        [&]
        {
            if (s == nullptr)
            {
                return null_handle("s");
            }
            if (lower == nullptr || upper == nullptr)
            {
                return null_argument("lower and upper");
            }
            if (n == 0)
            {
                return set_error(PO_INVALID_ARGUMENT, "n must be positive");
            }
            if (!all_valid(lower, upper, n))
            {
                return invalid_bounds();
            }
            std::vector<double> lower_copy(lower, lower + n);
            std::vector<double> upper_copy(upper, upper + n);
            s->spec.asset_lower = std::move(lower_copy);
            s->spec.asset_upper = std::move(upper_copy);
            s->spec.lower = -kInfinity;
            s->spec.upper = kInfinity;
            return PO_OK;
        });
}

extern "C" po_status po_spec_set_target_return(po_spec* s, double r)
{
    return guard(
        [&]
        {
            if (s == nullptr)
            {
                return null_handle("s");
            }
            if (!std::isfinite(r))
            {
                return set_error(PO_INVALID_ARGUMENT, "r must be finite");
            }
            s->spec.target_return = r;
            return PO_OK;
        });
}

extern "C" po_status po_spec_set_risk_aversion(po_spec* s, double delta)
{
    return guard(
        [&]
        {
            if (s == nullptr)
            {
                return null_handle("s");
            }
            if (!std::isfinite(delta) || delta <= 0.0)
            {
                return set_error(PO_INVALID_ARGUMENT, "delta must be positive and finite");
            }
            s->spec.risk_aversion = delta;
            return PO_OK;
        });
}

extern "C" po_status po_spec_set_importance(po_spec* s, double lambda)
{
    return guard(
        [&]
        {
            if (s == nullptr)
            {
                return null_handle("s");
            }
            if (!(lambda >= 0.0 && lambda < 1.0))
            {
                return set_error(PO_INVALID_ARGUMENT, "lambda must satisfy 0 <= lambda < 1");
            }
            s->spec.importance = lambda;
            return PO_OK;
        });
}

extern "C" po_status po_spec_set_max_volatility(po_spec* s, double sigma)
{
    return guard(
        [&]
        {
            if (s == nullptr)
            {
                return null_handle("s");
            }
            if (!std::isfinite(sigma) || sigma <= 0.0)
            {
                return set_error(PO_INVALID_ARGUMENT, "sigma must be positive and finite");
            }
            s->spec.max_volatility = sigma;
            return PO_OK;
        });
}

extern "C" void po_spec_free(po_spec* s)
{
    std::unique_ptr<po_spec> owner(s);
}
