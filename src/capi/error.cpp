#include "capi/error.hpp"

#include <string>

namespace po::capi
{

namespace
{

thread_local std::string last_error_message;

}

po_status set_error(po_status status, std::string_view message)
{
    last_error_message.assign(message);
    return status;
}

void clear_error()
{
    last_error_message.clear();
}

const char* last_error()
{
    return last_error_message.c_str();
}

}
