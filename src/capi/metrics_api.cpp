#include <cmath>
#include <cstddef>
#include <span>
#include <string>

#include "capi/error.hpp"
#include "capi/handles.hpp"
#include "metrics/portfolio.hpp"

using po::capi::guard;
using po::capi::null_argument;
using po::capi::null_handle;
using po::capi::set_error;

namespace
{

po_status check_weights(const po::Data& data, const double* w, std::size_t n)
{
    if (w == nullptr)
    {
        return null_argument("w");
    }
    if (n != data.weight_count())
    {
        return set_error(PO_INVALID_ARGUMENT, "n must be " + std::to_string(data.weight_count()) + ", got " +
                                                  std::to_string(n));
    }
    return PO_OK;
}

} // namespace

extern "C" po_status po_portfolio_return(const po_data* d, const double* w, size_t n, double* out)
{
    return guard(
        [&]
        {
            if (d == nullptr)
            {
                return null_handle("d");
            }
            if (out == nullptr)
            {
                return null_argument("out");
            }
            po_status status = check_weights(d->data, w, n);
            if (status != PO_OK)
            {
                return status;
            }
            if (!d->data.has_mu())
            {
                return set_error(PO_INVALID_ARGUMENT, "the return needs mu");
            }
            *out = po::metrics::portfolio_return(d->data, std::span<const double>(w, n));
            return PO_OK;
        });
}

extern "C" po_status po_portfolio_volatility(const po_data* d, const double* w, size_t n, double* out)
{
    return guard(
        [&]
        {
            if (d == nullptr)
            {
                return null_handle("d");
            }
            if (out == nullptr)
            {
                return null_argument("out");
            }
            po_status status = check_weights(d->data, w, n);
            if (status != PO_OK)
            {
                return status;
            }
            if (!d->data.has_covariance())
            {
                return set_error(PO_INVALID_ARGUMENT, "the volatility needs the covariance");
            }
            *out = po::metrics::portfolio_volatility(d->data, std::span<const double>(w, n));
            return PO_OK;
        });
}

extern "C" po_status po_utility(const po_data* d, const double* w, size_t n, double delta, double* out)
{
    return guard(
        [&]
        {
            if (d == nullptr)
            {
                return null_handle("d");
            }
            if (out == nullptr)
            {
                return null_argument("out");
            }
            po_status status = check_weights(d->data, w, n);
            if (status != PO_OK)
            {
                return status;
            }
            if (!d->data.has_mu())
            {
                return set_error(PO_INVALID_ARGUMENT, "the utility needs mu");
            }
            if (!d->data.has_covariance())
            {
                return set_error(PO_INVALID_ARGUMENT, "the utility needs the covariance");
            }
            if (!std::isfinite(delta) || delta <= 0.0)
            {
                return set_error(PO_INVALID_ARGUMENT, "delta must be positive and finite");
            }
            *out = po::metrics::utility(d->data, std::span<const double>(w, n), delta);
            return PO_OK;
        });
}

extern "C" po_status po_tracking_error(const po_data* d, const double* w, size_t n, double* out)
{
    return guard(
        [&]
        {
            if (d == nullptr)
            {
                return null_handle("d");
            }
            if (out == nullptr)
            {
                return null_argument("out");
            }
            po_status status = check_weights(d->data, w, n);
            if (status != PO_OK)
            {
                return status;
            }
            if (!d->data.has_returns())
            {
                return set_error(PO_INVALID_ARGUMENT, "the tracking error needs returns");
            }
            *out = po::metrics::tracking_error(d->data, std::span<const double>(w, n));
            return PO_OK;
        });
}
