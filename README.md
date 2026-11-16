# PortOpt

PortOpt is a portfolio optimization library. The numerical core is C++23; its only stable contract is a C API
(`include/portopt/portopt.h`), so other languages can bind to it, and a Python package, `portopt`, is the first
consumer. Optimization problems are quadratic programs, solved by a dense Goldfarb-Idnani dual active-set method
written from the paper; the library has no dependency besides GoogleTest (for its tests).

## Status

Version 0.1.0 (pre-release). Nine objectives are solved: minimum variance (`MIN_VARIANCE`), maximum return under a
volatility cap (`MAX_RETURN`), the maximum Sharpe ratio (`MAX_SHARPE`), mean-variance utility (`UTILITY`), the same
with an importance weight between risk and return (`IMPORTANCE`), the Kelly criterion (`KELLY`), tracking error
against a benchmark (`TRACKING_ERROR`), the maximum diversification ratio (`MAX_DIVERSIFICATION`) and risk parity with
a risk budget (`RISK_PARITY`, a local method: its result is a KKT point, not a proven global minimum). All take a box on
the weights (one pair for all of them or one pair per asset), an optional gross exposure limit (the leverage,
`sum |w_i| <= L` over the risky weights) and an optional risk-free asset (cash, always the last weight);
`MIN_VARIANCE` and `TRACKING_ERROR` also take a target return, and `RISK_PARITY` takes neither cash nor a leverage but a
target return and a target volatility. The result carries the weights and the multipliers of
the budget, target-return, volatility and leverage constraints and of the bounds, under the convention in
[docs/guide/math.md](docs/guide/math.md). A solve that ends without a point it can vouch for has status `MAX_ITER`,
with NaN weights, multipliers and residual, for every objective. A `MAX_SHARPE` or
`MAX_DIVERSIFICATION` problem whose supremum is not attained, such as `MAX_SHARPE` with unbounded short selling and a
risk-free rate above the minimum-variance return, has status `UNBOUNDED` and NaN weights. Portfolio return, volatility,
utility, tracking error and the Sharpe and diversification ratios of a weight vector are available too. Until 1.0 the C
ABI may change between releases. See [CHANGELOG.md](CHANGELOG.md).

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
# (the minimum volatility of this covariance is 16.4%, so a target below it is INFEASIBLE)
cov = np.array([[0.04, 0.006, 0.002], [0.006, 0.09, 0.012], [0.002, 0.012, 0.16]])
data = po.Data(None, cov)
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
return series instead of `mu` and `cov`: `po.Data(None, None, returns=returns, benchmark=benchmark)`, with `returns`
a periods-by-assets array. `data.utility(weights, risk_aversion)`, `data.tracking_error(weights)`,
`data.sharpe_ratio(weights)`, `data.diversification_ratio(weights)` and `data.risk_contributions(weights)` evaluate the
metrics, and the result has `leverage_multiplier` next to the other multipliers (NaN without a leverage, 0 when it is
inactive) and `residual`, the `RISK_PARITY` sum of squares `sum_i (c_i - b_i)^2` (0 for the other objectives).

An infeasible or unbounded problem is a result with `status == po.SolveStatus.INFEASIBLE` or
`po.SolveStatus.UNBOUNDED` and NaN weights, not an exception. With
`po.Data(mu, cov, risk_free=0.02)` the weights have one more entry, the cash, last. Invalid input raises
`ValueError`. `RISK_PARITY` with cash (`risk_free`) or a `leverage` is invalid input too.

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
- [docs/guide/math.md](docs/guide/math.md): notation, the QP of each objective, the leverage, the parametric path of
  `MAX_RETURN`, the homogenization of `MAX_SHARPE` and `MAX_DIVERSIFICATION`, `RISK_PARITY` and the multiplier convention
- [docs/guide/validation.md](docs/guide/validation.md): how results are checked, the tolerances and the fuzz battery
- [CONTRIBUTING.md](CONTRIBUTING.md) and [STYLEGUIDE.md](STYLEGUIDE.md)

## License

MIT, see [LICENSE](LICENSE).
