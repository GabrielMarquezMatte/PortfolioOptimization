# PortOpt

PortOpt is a portfolio optimization library. The numerical core is C++23; its only stable contract is a C API
(`include/portopt/portopt.h`), so other languages can bind to it, and a Python package, `portopt`, is the first
consumer (`pip install pyportopt`), with a header-only C++ wrapper as the second. Optimization problems are
quadratic programs, solved by a dense Goldfarb-Idnani dual active-set method written from the paper; the library has
no dependency besides GoogleTest (for its tests).

## Status

Version 0.1.0 (pre-release). Nine objectives are solved: minimum variance (`MIN_VARIANCE`), maximum return under a
volatility cap (`MAX_RETURN`), the maximum Sharpe ratio (`MAX_SHARPE`), mean-variance utility (`UTILITY`), the same
with an importance weight between risk and return (`IMPORTANCE`), the Kelly criterion (`KELLY`), tracking error
against a benchmark (`TRACKING_ERROR`), the maximum diversification ratio (`MAX_DIVERSIFICATION`) and risk parity with
a risk budget (`RISK_PARITY`, a local method: its result is a KKT point, not a proven global minimum). All take a box on
the weights (one pair for all of them or one pair per asset); all but `RISK_PARITY` take an optional gross exposure
limit (the leverage, `sum |w_i| <= L` over the risky weights) and an optional risk-free asset (cash, always the last
weight);
`MIN_VARIANCE` and `TRACKING_ERROR` also take a target return, and `RISK_PARITY` takes neither cash nor a leverage but a
target return and a target volatility. The result carries the weights and the multipliers of
the budget, target-return, volatility and leverage constraints and of the bounds, under the convention in
[docs/guide/math.md](docs/guide/math.md). A solve that ends without a point it can vouch for has status `NOT_CONVERGED`,
with NaN weights, multipliers and residual, for every objective. A `MAX_SHARPE` or
`MAX_DIVERSIFICATION` problem whose supremum is not attained, such as `MAX_SHARPE` with unbounded short selling and a
risk-free rate above the minimum-variance return, has status `UNBOUNDED` and NaN weights. Portfolio return, volatility,
utility, tracking error and the Sharpe and diversification ratios of a weight vector are available too. Until 1.0 the C
ABI may change between releases. See [CHANGELOG.md](CHANGELOG.md).

## Install

Python 3.11 or newer, with `numpy`:

```
pip install pyportopt
```

The distribution is `pyportopt` (the name `portopt` is taken on PyPI); the module is `portopt`. Wheels are built for
Windows x64, macOS arm64 and x64, and Linux (glibc and musl) x64 and arm64, each with the native library inside.

.NET 10 (C#), with the same seven native libraries inside, trimming- and Native AOT-compatible:

```
dotnet add package PortOpt
```

C and C++: every [GitHub Release](https://github.com/GabrielMarquezMatte/PortfolioOptimization/releases) has the
library of each platform, `portopt-native-<platform>.so`, `.dylib` or `.dll`, to be saved as `libportopt.so`,
`libportopt.dylib` or `libportopt.dll` (the name the loader looks for), and for Windows the import libraries
`portopt-native-win-x64.dll.a` (MinGW) and `portopt-native-win-x64.lib` (MSVC). The headers are
`include/portopt/portopt.h` (C) and `cpp/include/portopt/portopt.hpp` (C++23) of the tag. From a source tree,
`cmake --install` installs the library, both headers and a CMake package:

```
cmake --preset native-release
cmake --build --preset native-release
cmake --install build/native-release --prefix <prefix>
```

and a CMake project then uses it with `find_package(PortOpt 0.1 REQUIRED)` and the target `PortOpt::cpp` (or
`PortOpt::portopt` for the C API alone). Until 1.0 the package accepts only the same minor version, since the ABI may
change between minor versions. `add_subdirectory` of a checkout gives the same targets.

## Python example

```python
from dataclasses import replace

import numpy as np
import portopt as po

mu = np.array([0.05, 0.08, 0.11])
cov = np.array([[0.04, 0.006, 0.002], [0.006, 0.09, 0.012], [0.002, 0.012, 0.16]])
data = po.Data(mu, cov)

# Long-only minimum variance
result = po.solve(po.Spec(po.Objective.MIN_VARIANCE, bounds=(0.0, 1.0)), data)
print(result.status.name, result.weights)      # OPTIMAL [0.6234 0.2356 0.141 ] (rounded)

# The efficient frontier: one solve per target return
base = po.Spec(po.Objective.MIN_VARIANCE, bounds=(0.0, 1.0))
frontier = []
for target in np.linspace(0.06, 0.10, 5):
    weights = po.solve(replace(base, target_return=target), data).weights
    frontier.append((data.volatility(weights), data.portfolio_return(weights)))

# Maximum return with the volatility capped at 20%
spec = po.Spec(po.Objective.MAX_RETURN, bounds=(0.0, 1.0), max_volatility=0.20)
result = po.solve(spec, data)
print(result.weights)                          # [0.2838 0.3578 0.3584] (rounded)
print(data.portfolio_return(result.weights))   # 0.0822 (rounded)
print(result.volatility_multiplier)            # -0.2501 (rounded): the change of -return per unit of cap

# Maximum Sharpe ratio with short selling, and the gross exposure capped at 2
mu = np.array([0.10, 0.04, 0.08])
cov = np.array([[0.04, 0.057, 0.0], [0.057, 0.09, 0.0], [0.0, 0.0, 0.16]])
data = po.Data(mu, cov, risk_free=0.02)
free = po.solve(po.Spec(po.Objective.MAX_SHARPE), data)
capped = po.solve(po.Spec(po.Objective.MAX_SHARPE, leverage=2.0), data)
print(np.abs(free.weights[:3]).sum())          # 4.0926 (rounded): the tangency portfolio, uncapped
print(capped.weights)                          # [ 1.3631 -0.5  0.1369  0. ] (rounded): cash weight 0
print(data.sharpe_ratio(capped.weights))       # 0.7211 (rounded), against 1.0906 uncapped
print(capped.leverage_multiplier)              # -0.3105 (rounded): the change of -Sharpe per unit of L

# Risk parity with a risk budget (50% / 30% / 20% of the risk) and the volatility fixed at 17%
# (the minimum volatility of this covariance is 16.3%, so a target below it is INFEASIBLE)
cov = np.array([[0.04, 0.006, 0.002], [0.006, 0.09, 0.012], [0.002, 0.012, 0.16]])
data = po.Data(cov=cov)
budget = [0.5, 0.3, 0.2]
spec = po.Spec(po.Objective.RISK_PARITY, risk_budget=budget, target_volatility=0.17)
result = po.solve(spec, data)
print(result.status.name, result.weights)      # OPTIMAL [0.4754 0.3046 0.22  ] (rounded)
print(data.volatility(result.weights))         # 0.17
print(data.risk_contributions(result.weights)) # [0.3501 0.3468 0.3031] (rounded): the best the target volatility allows
print(result.residual)                         # 0.0353 (rounded), the sum of (c_i - b_i)^2; 0 without the target
exact = po.solve(po.Spec(po.Objective.RISK_PARITY, risk_budget=budget), data)
print(data.risk_contributions(exact.weights), exact.residual)   # [0.5 0.3 0.2], ~0
```

`Spec` is a frozen dataclass: `dataclasses.replace` builds a new, validated spec (`from dataclasses import replace`).
Its keyword fields are `bounds`, `target_return`, `risk_aversion` (`UTILITY`), `importance` (`IMPORTANCE`),
`max_volatility` (`MAX_RETURN`), `leverage` (every objective but `RISK_PARITY`), `target_volatility` and `risk_budget`
(`RISK_PARITY`; the budget is stored as a tuple, so a `Spec` still compares and hashes by value); a value the C setter rejects
raises `ValueError` when the spec is built, and a
setter or input that does not apply to the objective raises it from `po.solve`. `TRACKING_ERROR` takes the
return series instead of `mu` and `cov`: `po.Data(returns=returns, benchmark=benchmark)`, with `returns`
a periods-by-assets array. `data.utility(weights, risk_aversion)`, `data.tracking_error(weights)`,
`data.sharpe_ratio(weights)`, `data.diversification_ratio(weights)` and `data.risk_contributions(weights)` evaluate the
metrics, and the result has `leverage_multiplier` next to the other multipliers (NaN without a leverage, 0 when it is
inactive) and `residual`, the `RISK_PARITY` sum of squares `sum_i (c_i - b_i)^2` (0 for the other objectives).

An infeasible or unbounded problem is a result with `status == po.SolveStatus.INFEASIBLE` or
`po.SolveStatus.UNBOUNDED` and NaN weights, not an exception. With
`po.Data(mu, cov, risk_free=0.02)` the weights have one more entry, the cash, last. Invalid input raises
`ValueError`. `RISK_PARITY` with cash (`risk_free`) or a `leverage` is invalid input too. `portopt.__version__` is the
installed version.

## C# example

```csharp
using PortOpt;

double[] mu = [0.08, 0.12, 0.10];
double[,] covariance = { { 0.040, 0.006, 0.010 }, { 0.006, 0.090, 0.012 }, { 0.010, 0.012, 0.060 } };
using var data = new MarketData(mu, covariance, riskFree: 0.02);

// Long-only maximum Sharpe ratio
var result = Optimizer.Solve(new Spec { Objective = Objective.MaxSharpe, Bounds = new Bounds(0, 1) }, data);
Console.WriteLine($"{result.Status}: {string.Join(", ", result.Weights)}, Sharpe {data.SharpeRatio(result.Weights)}");

// Risk parity with a risk budget, no cash
using var risky = new MarketData(mu, covariance);
var parity = Optimizer.Solve(new Spec { Objective = Objective.RiskParity, RiskBudget = [0.5, 0.3, 0.2] }, risky);
Console.WriteLine(string.Join(", ", risky.RiskContributions(parity.Weights)));
```

`Spec` is an immutable record of the settings, `MarketData` owns the native data (dispose it), `Result` is plain
values, and a native failure is a `PortOptException` with its `Status` and the library's message; an infeasible problem
is a result with `SolveStatus.Infeasible`. More in [csharp/src/PortOpt/README.md](csharp/src/PortOpt/README.md).

## C++ example

```cpp
#include <cstdio>
#include <portopt/portopt.hpp>
#include <vector>

int main()
{
    const std::vector<double> mu{0.08, 0.12, 0.10};
    const std::vector<double> cov{0.040, 0.006, 0.010, 0.006, 0.090, 0.012, 0.010, 0.012, 0.060};
    po::Expected<po::Data> data = po::Data::create(3, mu, cov);   // cov is 3 x 3, row-major
    if (!data || !data->set_risk_free(0.02))
    {
        return 1;
    }
    const po::Expected<po::Result> result =
        po::solve(po::Spec{.objective = po::Objective::max_sharpe, .bounds = po::Bounds{0.0, 1.0}}, *data);
    if (!result)
    {
        std::fprintf(stderr, "%s\n", result.error().message.c_str()); // the status and the message of the C API
        return 1;
    }
    for (const double weight : result->weights) // 3 assets and the cash, last
    {
        std::printf("%.4f\n", weight); // 0.3749 0.3031 0.3221 0.0000 (rounded)
    }
}
```

`po::Spec` takes the same fields as the Python `Spec`, `po::Data` the same metrics as the Python `Data`, and every
call returns `std::expected`: an infeasible problem is a `Result` with `status == po::SolveStatus::infeasible`, not
an error. See [cpp/README.md](cpp/README.md).

## Build

CMake 3.28 or newer and a C++23 compiler (GCC 15 or newer, MSVC, or Apple Clang 15 or newer). Presets are in
`CMakePresets.json`:

```
cmake --preset gcc-release
cmake --build --preset gcc-release
ctest --preset gcc-release
```

`native-release` uses the default compiler and sets no `-march` (it is what CI uses; for a `-march=native` build see
[docs/performance/measuring.md](docs/performance/measuring.md#the-native-build)). The library is `libportopt.dll` (MinGW),
`Release/portopt.dll` (MSVC), `libportopt.so` or `libportopt.dylib`, under `build/<preset>/`.

## Performance

Measured on one desktop (AMD Ryzen 7 5700X, GCC 16.2, generic x86-64 build), the solve time against the same code before
the optimizations of version 0.1, same session, results bit for bit identical: the geometric mean over the 35 benchmark
cases is 0.71 at 10 assets, 0.29 at 100 and 0.23 at 1000. A 1000-asset `RISK_PARITY` with a box takes 21.5 s (86.8 s
before), a `MIN_VARIANCE` with a box 0.19 s (1.03 s), a 100-asset `MAX_RETURN` with a leverage 7.7 ms (18.8 ms). The C
interface adds under 0.3 µs to a call; the Python package adds about 30 µs per `Spec`, `Data` and `solve`, which matters
only for many solves of about 10 assets. Method, tables and what was left out:
[docs/performance](docs/performance/benchmarks.md). Compare only numbers from the same machine.

## Python package

The package is managed with [uv](https://docs.astral.sh/uv/) and needs Python 3.11 or newer and `numpy`. It loads the
library you built, found through `PORTOPT_LIB`:

```
PORTOPT_LIB=$PWD/build/gcc-release/libportopt.dll uv run --project python pytest python/tests -q
```

Without `PORTOPT_LIB` it looks in `python/src/portopt/_lib`.

## C# binding

`csharp/` holds the .NET 10 binding (`csharp/PortOpt.slnx`): the package, its xUnit tests and a BenchmarkDotNet project.
The tests load the library named by `PORTOPT_LIB`:

```
PORTOPT_LIB=$PWD/build/gcc-release/libportopt.dll dotnet test --solution csharp/PortOpt.slnx -c Release
```

## Documents

- [ARCHITECTURE.md](ARCHITECTURE.md): the layers, the files and the anatomy of a solve
- [docs/guide/math.md](docs/guide/math.md): notation, the QP of each objective, the leverage, the parametric path of
  `MAX_RETURN`, the homogenization of `MAX_SHARPE` and `MAX_DIVERSIFICATION`, `RISK_PARITY` and the multiplier convention
- [docs/guide/validation.md](docs/guide/validation.md): how results are checked, the tolerances and the fuzz battery
- [docs/guide/versioning.md](docs/guide/versioning.md): what a minor version may change, and the ABI version
- [docs/performance/benchmarks.md](docs/performance/benchmarks.md), [measuring.md](docs/performance/measuring.md) and
  [python-vs-arch.md](docs/performance/python-vs-arch.md): the benchmark results, how to measure and the cost of the C
  and Python layers
- [CONTRIBUTING.md](CONTRIBUTING.md) and [STYLEGUIDE.md](STYLEGUIDE.md)

## License

MIT, see [LICENSE](LICENSE).
