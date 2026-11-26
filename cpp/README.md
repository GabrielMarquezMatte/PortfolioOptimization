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

Estimation, from a return panel to a solve with one relative view (asset 0 beats asset 1 by 1%); each call returns an
`Expected`, dereferenced here for brevity, and `po::Views` owns its handle like `po::Data`:

```cpp
auto moments = po::sample_moments(returns, 3);                     // returns: periods x 3, daily, row-major
for (double& value : moments->cov)
{
    value *= 252;                                                  // annualize after the estimation
}
const auto sigma = po::Data::create(3, {}, moments->cov);          // covariance only: no mu, no risk-free rate
const auto pi = po::bl_equilibrium(*sigma, w_mkt, 2.5);            // in excess of the rate; w_mkt sums to 1
const auto prior = po::Data::create(3, *pi, moments->cov);        // its mu is Pi
auto views = po::Views::create(3, p, q);                           // p = {1, -1, 0}, q = {0.01}
if (!views->set_omega(*po::bl_omega_he_litterman(*sigma, *views, 0.05)))
{
    return 1;
}
auto posterior = po::bl_posterior(*prior, *views, 0.05);           // posterior->mu, posterior->cov
for (double& value : posterior->mu)
{
    value += risk_free;                                            // the posterior is in excess of the rate
}
auto data = po::Data::create(3, posterior->mu, posterior->cov);
if (!data->set_risk_free(risk_free))
{
    return 1;
}
const auto result = po::solve(po::Spec{.objective = po::Objective::max_sharpe, .bounds = po::Bounds{0.0, 1.0}}, *data);
```

A `po::Data` with a risk-free rate is refused by every estimation function, so set it only on the data you solve.
`po::bayes_stein(data, periods)` shrinks the `mu` of a `Data` built from `sample_moments` (compute it on the daily data
and annualize afterwards: its weight depends on the frequency).

Link it through CMake as the target `PortOpt::cpp`, which brings the shared library `PortOpt::portopt` with it,
either from a build of this repository (`add_subdirectory`) or from an install (`find_package(PortOpt)`).

The tests are `cpp/tests/wrapper_test.cpp` (built with the library's tests) and `cpp/scripts/check_abi_coverage.py`,
which fails if the wrapper leaves a function of the C header unused.
