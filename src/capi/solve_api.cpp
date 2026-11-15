#include <memory>
#include <utility>

#include "capi/error.hpp"
#include "capi/handles.hpp"
#include "problems/solve.hpp"

using po::capi::guard;
using po::capi::null_argument;
using po::capi::null_handle;

extern "C" po_status po_solve(const po_spec* s, const po_data* d, po_result** out)
{
    if (out != nullptr)
    {
        *out = nullptr;
    }
    return guard(
        [&]
        {
            if (s == nullptr)
            {
                return null_handle("s");
            }
            if (d == nullptr)
            {
                return null_handle("d");
            }
            if (out == nullptr)
            {
                return null_argument("out");
            }
            auto solved = po::solve(s->spec, d->data);
            if (!solved.has_value())
            {
                return po::capi::from_failure(solved.error());
            }
            *out = std::make_unique<po_result>(std::move(*solved)).release();
            return PO_OK;
        });
}
