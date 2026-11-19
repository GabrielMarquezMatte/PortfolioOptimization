#include "capi/error.hpp"

#include <string>

namespace po::capi
{

namespace
{

thread_local std::string last_error_message;

}

po_status set_error(po_status status, std::string_view message) noexcept
{
    try
    {
        last_error_message.assign(message);
    }
    catch (...)
    {
        // A std::string never shrinks below its small buffer, which holds this literal: it cannot allocate.
        last_error_message.assign("out of memory");
    }
    return status;
}

void clear_error() noexcept
{
    last_error_message.clear();
}

const char* last_error()
{
    return last_error_message.c_str();
}

}
