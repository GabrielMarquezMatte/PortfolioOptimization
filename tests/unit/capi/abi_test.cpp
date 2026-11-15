#include <gtest/gtest.h>
#include <new>
#include <portopt/portopt.h>
#include <stdexcept>

#include "capi/error.hpp"

TEST(Abi, VersionMatchesHeader)
{
    EXPECT_EQ(po_abi_version(), PO_ABI_VERSION);
}

TEST(Abi, VersionStringIsProjectVersion)
{
    EXPECT_STREQ(po_version(), "0.1.0");
}

TEST(Abi, LastErrorStartsEmpty)
{
    EXPECT_STREQ(po_last_error(), "");
}

TEST(Guard, BadAllocBecomesOutOfMemory)
{
    EXPECT_EQ(po::capi::guard([]() -> po_status { throw std::bad_alloc(); }), PO_OUT_OF_MEMORY);
    EXPECT_STREQ(po_last_error(), "out of memory");
}

TEST(Guard, OtherExceptionBecomesInternalError)
{
    EXPECT_EQ(po::capi::guard([]() -> po_status { throw std::runtime_error("boom"); }), PO_INTERNAL_ERROR);
    EXPECT_STREQ(po_last_error(), "boom");
}

TEST(Guard, SuccessClearsLastError)
{
    po::capi::set_error(PO_INVALID_ARGUMENT, "x");
    EXPECT_EQ(po::capi::guard([] { return PO_OK; }), PO_OK);
    EXPECT_STREQ(po_last_error(), "");
}
