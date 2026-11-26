# PortOpt

Portfolio optimization for .NET on a native C++23 engine: minimum variance, maximum return under a volatility cap,
maximum Sharpe, utility, importance, Kelly, tracking error, maximum diversification and risk parity, with scalar or
per-weight bounds, a gross exposure limit (leverage) and an optional risk-free asset. The package carries the
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
using var data = new MarketData(new MarketInputs { Mu = mu, Covariance = covariance, RiskFree = 0.02 });
var result = Optimizer.Solve(new Spec { Objective = Objective.MaxSharpe, Bounds = new Bounds(0, 1) }, data);
Console.WriteLine($"{result.Status}: {string.Join(", ", result.Weights)}");
Console.WriteLine($"Sharpe {data.SharpeRatio(result.Weights)}");
```

Estimation: the same flow from a return panel to a solve, with one relative view (asset 0 beats asset 1 by 1%).

```csharp
double[,] returns = LoadDailyReturns();                                    // periods by assets, daily
var covariance = Estimation.SampleMoments(returns).Covariance;
for (var i = 0; i < 3; i++)
{
    for (var j = 0; j < 3; j++)
    {
        covariance[i, j] *= 252;                                           // annualize after the estimation
    }
}
var pi = Estimation.Equilibrium(covariance, [0.5, 0.3, 0.2], delta: 2.5);  // in excess of the rate
double[,] p = { { 1.0, -1.0, 0.0 } };
var omega = Estimation.OmegaHeLitterman(covariance, p, tau: 0.05);
var posterior = Estimation.Posterior(covariance, pi, 0.05, new Views(p, [0.01], omega));
using var data = new MarketData(posterior.ToMarketInputs(riskFree: 0.02)); // adds the rate back, sets the cash
var result = Optimizer.Solve(new Spec { Objective = Objective.MaxSharpe, Bounds = new Bounds(0, 1) }, data);
```

`Pi`, `Q` and the posterior mean are in excess of the risk-free rate; `ToMarketInputs` adds it back (without a rate the
mean is unchanged and there is no cash). `Estimation.BayesStein(muHat, covariance, periods)` shrinks a mean (compute it
on the daily data and annualize afterwards: its weight depends on the frequency), and `OmegaIdzorek` is the other
`Omega`. Each function also has an overload on a `MarketData` (with no risk-free rate) in place of the covariance, and
of the mean where it reads one, so that a covariance used by several calls is validated and factored once.

- `Spec` is an immutable record: the objective and the settings it takes (`TargetReturn`, `RiskAversion`,
  `Importance`, `MaxVolatility`, `Leverage`, `TargetVolatility`, `RiskBudget`, `Bounds` or `WeightBounds`).
- `MarketData` holds a native handle: dispose it (`using`). With a risk-free rate the cash weight comes last.
- `Result` is plain values: status, weights, iterations, the multipliers and the residual.
- A native failure is a `PortOptException` with its `Status` and the library's message; an infeasible problem is a
  result with `SolveStatus.Infeasible`, not an exception.
- `PORTOPT_LIB` loads a local build of the native library instead of the packaged one.

The C header `portopt.h` is the reference for every objective, setter and multiplier:
https://github.com/GabrielMarquezMatte/PortfolioOptimization
