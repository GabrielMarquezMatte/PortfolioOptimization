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

bool product_fits(std::size_t a, std::size_t b)
{
    return a == 0 || b <= std::numeric_limits<std::size_t>::max() / a;
}

} // namespace

extern "C" po_status po_views_create(size_t k, size_t n, const double* p, const double* q, po_views** out)
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
            if (!product_fits(k, n))
            {
                return set_error(PO_INVALID_ARGUMENT, "views size is too large");
            }
            if (k >= 1 && p == nullptr)
            {
                return null_argument("p");
            }
            if (k >= 1 && q == nullptr)
            {
                return null_argument("q");
            }
            auto created = po::estimation::Views::create(k, n, view_of(p, k * n), view_of(q, k));
            if (!created.has_value())
            {
                return po::capi::from_failure(created.error());
            }
            *out = std::make_unique<po_views>(std::move(*created)).release();
            return PO_OK;
        });
}

extern "C" po_status po_views_set_omega(po_views* v, const double* omega, size_t k)
{
    return guard(
        [&]
        {
            if (v == nullptr)
            {
                return null_handle("v");
            }
            if (omega == nullptr && k != 0)
            {
                return null_argument("omega");
            }
            if (!product_fits(k, k))
            {
                return set_error(PO_INVALID_ARGUMENT, "omega size is too large");
            }
            auto stored = v->views.set_omega(view_of(omega, k * k));
            if (!stored.has_value())
            {
                return po::capi::from_failure(stored.error());
            }
            return PO_OK;
        });
}

extern "C" void po_views_free(po_views* v)
{
    std::unique_ptr<po_views> owner(v);
}
