#include <portopt/portopt.h>

#include "capi/error.hpp"

extern "C" int32_t po_abi_version(void)
{
    return PO_ABI_VERSION;
}

extern "C" const char* po_version(void)
{
    return PO_VERSION_STRING;
}

extern "C" const char* po_last_error(void)
{
    return po::capi::last_error();
}
