# PortOpt

PortOpt is a portfolio optimization library. The numerical core is C++23; its only stable contract is a C API
(`include/portopt/portopt.h`), so other languages can bind to it, and a Python package, `portopt`, is the first
consumer. Optimization problems are quadratic programs, solved by a dense Goldfarb-Idnani dual active-set method
written from the paper; the library has no dependency besides GoogleTest (for its tests).

## Status

Version 0.1.0 (pre-release). One objective is solved so far: minimum variance (`MIN_VARIANCE`), with a box on the
weights (one pair for all of them or one pair per asset), a target return and an optional risk-free asset (cash,
always the last weight). The result carries the weights and the multipliers of the budget and target-return
constraints and of the bounds, under the convention in [docs/guide/math.md](docs/guide/math.md). The other objectives
of the enum return `PO_UNSUPPORTED`. Portfolio return and volatility of a weight vector are available too.
Until 1.0 the C ABI may change between releases. See [CHANGELOG.md](CHANGELOG.md).

## Python example

```python
import numpy as np
import portopt as po

mu = np.array([0.05, 0.08, 0.11])
cov = np.array([[0.04, 0.006, 0.002], [0.006, 0.09, 0.012], [0.002, 0.012, 0.16]])
data = po.Data(mu, cov)

# Long-only minimum variance
result = po.solve(po.Spec(po.Objective.MIN_VARIANCE, bounds=(0.0, 1.0)), data)
print(result.status.name, result.weights)      # OPTIMAL [0.6234 0.2356 0.141 ] (rounded)

# The efficient frontier: one solve per target return
frontier = []
for target in np.linspace(0.06, 0.10, 5):
    spec = po.Spec(po.Objective.MIN_VARIANCE, bounds=(0.0, 1.0), target_return=target)
    weights = po.solve(spec, data).weights
    frontier.append((data.volatility(weights), data.portfolio_return(weights)))
```

An infeasible problem is a result with `status == po.SolveStatus.INFEASIBLE` and NaN weights, not an exception. With
`po.Data(mu, cov, risk_free=0.02)` the weights have one more entry, the cash, last. Invalid input raises
`ValueError`, an objective that is not implemented `NotImplementedError`.

## Build

CMake 3.28 or newer and a C++23 compiler (GCC 15 or newer, MSVC, or Apple Clang 15 or newer). Presets are in
`CMakePresets.json`:

```
cmake --preset gcc-release
cmake --build --preset gcc-release
ctest --preset gcc-release
```

`native-release` uses the default compiler and is what CI uses. The library is `libportopt.dll` (MinGW),
`Release/portopt.dll` (MSVC), `libportopt.so` or `libportopt.dylib`, under `build/<preset>/`.

## Python package

The package is managed with [uv](https://docs.astral.sh/uv/) and needs Python 3.11 or newer and `numpy`. It loads the
library you built, found through `PORTOPT_LIB`:

```
PORTOPT_LIB=$PWD/build/gcc-release/libportopt.dll uv run --project python pytest python/tests -q
```

Without `PORTOPT_LIB` it looks in `python/src/portopt/_lib`.

## Documents

- [ARCHITECTURE.md](ARCHITECTURE.md): the layers, the files and the anatomy of a solve
- [docs/guide/math.md](docs/guide/math.md): notation, the QP of `MIN_VARIANCE` and the multiplier convention
- [docs/guide/validation.md](docs/guide/validation.md): how results are checked, and the tolerances
- [CONTRIBUTING.md](CONTRIBUTING.md) and [STYLEGUIDE.md](STYLEGUIDE.md)

## License

MIT, see [LICENSE](LICENSE).
