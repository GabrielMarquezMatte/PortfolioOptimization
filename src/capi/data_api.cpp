#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <span>
#include <utility>

#include "capi/error.hpp"
#include "capi/handles.hpp"

using po::capi::guard;
using po::capi::null_argument;
using po::capi::null_handle;
using po::capi::set_error;

namespace
{

std::span<const double> view_of(const double* values, std::size_t count)
{
    if (values == nullptr)
    {
        return {};
    }
    return std::span<const double>(values, count);
}

bool square_fits(std::size_t n)
{
    return n == 0 || n <= std::numeric_limits<std::size_t>::max() / n;
}

} // namespace

extern "C" po_status po_data_create(size_t n, const double* mu, const double* cov, po_data** out)
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
            if (!square_fits(n))
            {
                return set_error(PO_INVALID_ARGUMENT, "n is too large");
            }
            auto created = po::Data::create(n, view_of(mu, n), view_of(cov, n * n));
            if (!created.has_value())
            {
                return po::capi::from_failure(created.error());
            }
            *out = std::make_unique<po_data>(std::move(*created)).release();
            return PO_OK;
        });
}

extern "C" po_status po_data_set_risk_free(po_data* d, double rate)
{
    return guard(
        [&]
        {
            if (d == nullptr)
            {
                return null_handle("d");
            }
            if (!std::isfinite(rate))
            {
                return set_error(PO_INVALID_ARGUMENT, "rate must be finite");
            }
            d->data.set_risk_free(rate);
            return PO_OK;
        });
}

extern "C" po_status po_data_set_returns(po_data* d, const double* returns, const double* benchmark, size_t t)
{
    return guard(
        [&]
        {
            if (d == nullptr)
            {
                return null_handle("d");
            }
            if (returns == nullptr || benchmark == nullptr)
            {
                return null_argument("returns and benchmark");
            }
            std::size_t n = d->data.asset_count();
            if (t > std::numeric_limits<std::size_t>::max() / n)
            {
                return set_error(PO_INVALID_ARGUMENT, "returns size is too large");
            }
            auto stored = d->data.set_returns(view_of(returns, t * n), view_of(benchmark, t), t);
            if (!stored.has_value())
            {
                return po::capi::from_failure(stored.error());
            }
            return PO_OK;
        });
}

extern "C" po_status po_data_asset_count(const po_data* d, size_t* out)
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
            *out = d->data.asset_count();
            return PO_OK;
        });
}

extern "C" po_status po_data_weight_count(const po_data* d, size_t* out)
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
            *out = d->data.weight_count();
            return PO_OK;
        });
}

extern "C" void po_data_free(po_data* d)
{
    std::unique_ptr<po_data> owner(d);
}
