#pragma once

#include <cstdint>

namespace po::bench
{

// Calls of the global operator new since the start of the process. Defined only in the benchmark executables,
// which replace the global allocation functions to count them; the library and the tests are not affected.
std::uint64_t allocation_count();

// Checks the counter against allocations whose number is known; prints each failure and returns the number of
// failures.
int allocation_selftest();

} // namespace po::bench
