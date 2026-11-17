# PortOpt for C++

`include/portopt/portopt.hpp` is a header-only C++23 wrapper over the C interface (`include/portopt/portopt.h`):
`std::expected` for errors (the status and the message of `po_last_error()`), `std::span` for vectors and matrices,
and a move-only `po::Data` that frees its handle. The names are those of the Python package. It translates and
nothing else, so a solve through it is the C solve bit for bit.

```cpp
#include <portopt/portopt.hpp>

const auto data = po::Data::create(3, mu, cov);                  // mu: 3 entries, cov: 3 x 3 row-major
const auto result = po::solve(po::Spec{.objective = po::Objective::min_variance, .bounds = po::Bounds{0.0, 1.0}},
                              *data);
if (!result)
{
    std::fprintf(stderr, "%s\n", result.error().message.c_str());
}
```

Link it through CMake as the target `PortOpt::cpp`, which brings the shared library `PortOpt::portopt` with it,
either from a build of this repository (`add_subdirectory`) or from an install (`find_package(PortOpt)`).

The tests are `cpp/tests/wrapper_test.cpp` (built with the library's tests) and `cpp/scripts/check_abi_coverage.py`,
which fails if the wrapper leaves a function of the C header unused.
