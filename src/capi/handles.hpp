#pragma once

#include <portopt/portopt.h>
#include <string>

#include "capi/error.hpp"
#include "core/data.hpp"
#include "core/failure.hpp"
#include "core/result.hpp"
#include "core/spec.hpp"
#include "estimation/views.hpp"

struct po_spec
{
    po::Spec spec;
};

struct po_data
{
    po::Data data;
};

struct po_result
{
    po::Result result;
};

struct po_views
{
    po::estimation::Views views;
};

namespace po::capi
{

inline po_status null_handle(const char* name)
{
    return set_error(PO_INVALID_HANDLE, std::string(name) + " must not be null");
}

inline po_status null_argument(const char* name)
{
    return set_error(PO_INVALID_ARGUMENT, std::string(name) + " must not be null");
}

inline po_status from_failure(const Failure& failure)
{
    if (failure.kind == FailureKind::Unsupported)
    {
        return set_error(PO_UNSUPPORTED, failure.message);
    }
    return set_error(PO_INVALID_ARGUMENT, failure.message);
}

} // namespace po::capi
