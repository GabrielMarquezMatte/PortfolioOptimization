# PortOpt

Portfolio optimization for .NET on a native C++23 engine: minimum variance, maximum return under a volatility cap,
maximum Sharpe, utility, importance, Kelly, tracking error, maximum diversification and risk parity, with scalar or
per-asset weight bounds, a gross exposure limit (leverage) and an optional risk-free asset. The package carries the
native library for linux-x64, linux-arm64, linux-musl-x64, linux-musl-arm64, osx-arm64, osx-x64 and win-x64, and is
trimming- and Native AOT-compatible.

```csharp
using PortOpt;

double[] mu = [0.08, 0.12, 0.10];
double[,] covariance =
{
    { 0.040, 0.006, 0.010 },
    { 0.006, 0.090, 0.012 },
    { 0.010, 0.012, 0.060 },
};
using var data = new MarketData(mu, covariance, riskFree: 0.02);
var result = Optimizer.Solve(new Spec { Objective = Objective.MaxSharpe, Bounds = new Bounds(0, 1) }, data);
Console.WriteLine($"{result.Status}: {string.Join(", ", result.Weights)}");
Console.WriteLine($"Sharpe {data.SharpeRatio(result.Weights)}");
```

- `Spec` is an immutable record: the objective and the settings it takes (`TargetReturn`, `RiskAversion`,
  `Importance`, `MaxVolatility`, `Leverage`, `TargetVolatility`, `RiskBudget`, `Bounds` or `AssetBounds`).
- `MarketData` holds a native handle: dispose it (`using`). With a risk-free rate the cash weight comes last.
- `Result` is plain values: status, weights, iterations, the multipliers and the residual.
- A native failure is a `PortOptException` with its `Status` and the library's message; an infeasible problem is a
  result with `SolveStatus.Infeasible`, not an exception.
- `PORTOPT_LIB` loads a local build of the native library instead of the packaged one.

The C header `portopt.h` is the reference for every objective, setter and multiplier:
https://github.com/GabrielMarquezMatte/PortfolioOptimization
