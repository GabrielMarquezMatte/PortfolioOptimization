# PortOpt

PortOpt is a portfolio optimization library. The numerical core is C++23; its only stable contract is a C API
(`include/portopt/portopt.h`), so other languages can bind to it, and a Python package, `portopt`, is the first
consumer. Optimization problems are quadratic programs, solved by a dense Goldfarb-Idnani dual active-set method
written from the paper; the library has no dependency besides GoogleTest (for its tests).

## Status

Version 0.1.0 (pre-release). Six objectives are solved: minimum variance (`MIN_VARIANCE`), maximum return under a
volatility cap (`MAX_RETURN`), mean-variance utility (`UTILITY`), the same with an importance weight between risk
and return (`IMPORTANCE`), the Kelly criterion (`KELLY`) and tracking error against a benchmark
(`TRACKING_ERROR`). All take a box on the weights (one pair for all of them or one pair per asset) and an optional
risk-free asset (cash, always the last weight); `MIN_VARIANCE` and `TRACKING_ERROR` also take a target return. The
result carries the weights and the multipliers of the budget, target-return and volatility constraints and of the
bounds, under the convention in [docs/guide/math.md](docs/guide/math.md). The other objectives of the enum return
`PO_UNSUPPORTED`. Portfolio return, volatility, utility and tracking error of a weight vector are available too.
Until 1.0 the C ABI may change between releases. See [CHANGELOG.md](CHANGELOG.md).

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
```

`Spec` is a frozen dataclass: `dataclasses.replace` builds a new, validated spec (`from dataclasses import replace`).
Its keyword fields are `bounds`, `target_return`, `risk_aversion` (`UTILITY`), `importance` (`IMPORTANCE`) and
`max_volatility` (`MAX_RETURN`); a value the C setter rejects raises `ValueError` when the spec is built, and a
setter or input that does not apply to the objective raises it from `po.solve`. `TRACKING_ERROR` takes the
return series instead of `mu` and `cov`: `po.Data(None, None, returns=returns, benchmark=benchmark)`, with `returns`
a periods-by-assets array. `data.utility(weights, risk_aversion)` and `data.tracking_error(weights)` evaluate the
two metrics.

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
- [docs/guide/math.md](docs/guide/math.md): notation, the QP of each objective, the parametric path of `MAX_RETURN`
  and the multiplier convention
- [docs/guide/validation.md](docs/guide/validation.md): how results are checked, and the tolerances
- [CONTRIBUTING.md](CONTRIBUTING.md) and [STYLEGUIDE.md](STYLEGUIDE.md)

## License

MIT, see [LICENSE](LICENSE).
