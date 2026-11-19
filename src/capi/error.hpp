#pragma once

#include <exception>
#include <new>
#include <portopt/portopt.h>
#include <string_view>

namespace po::capi
{

po_status set_error(po_status status, std::string_view message) noexcept;
void clear_error() noexcept;
const char* last_error();

template <typename F>
po_status guard(F&& callable)
{
    clear_error();
    try
    {
        return callable();
    }
    catch (const std::bad_alloc&)
    {
        return set_error(PO_OUT_OF_MEMORY, "out of memory");
    }
    catch (const std::exception& error)
    {
        return set_error(PO_INTERNAL_ERROR, error.what());
    }
    catch (...)
    {
        return set_error(PO_INTERNAL_ERROR, "unknown error");
    }
}

}
