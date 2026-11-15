#pragma once

#include <string>

namespace po
{

enum class FailureKind
{
    InvalidArgument,
    Unsupported
};

struct Failure
{
    FailureKind kind;
    std::string message;
};

} // namespace po
