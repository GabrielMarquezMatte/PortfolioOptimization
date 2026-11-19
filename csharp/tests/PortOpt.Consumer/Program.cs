using System.Globalization;
using PortOpt;

// Solves a long-only maximum Sharpe portfolio and exits 0 only when it is optimal: the package found its native library
// by itself.
double[] mu = [0.08, 0.12, 0.10];
double[,] covariance =
{
    { 0.040, 0.006, 0.010 },
    { 0.006, 0.090, 0.012 },
    { 0.010, 0.012, 0.060 },
};
using var data = new MarketData(mu, covariance, riskFree: 0.02);
var result = Optimizer.Solve(new Spec { Objective = Objective.MaxSharpe, Bounds = new Bounds(0, 1) }, data);
Console.Out.WriteLine(string.Create(CultureInfo.InvariantCulture,
    $"PortOpt {NativeAbi.LibraryVersion} (ABI {NativeAbi.Version}): {result.Status}, Sharpe {data.SharpeRatio(result.Weights)}"));
return result.Status == SolveStatus.Optimal ? 0 : 1;
