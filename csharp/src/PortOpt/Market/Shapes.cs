namespace PortOpt;

// The checks of the input shapes the native side cannot see (it takes pointers and one count), as the Python package
// does before it calls po_data_create.
internal static class Shapes
{
    public static int AssetCount(double[]? mu, double[,]? covariance, double[,]? returns, double[]? benchmark)
    {
        CheckSeries(returns, benchmark);
        var sizes = new HashSet<int>();
        if (mu is not null)
        {
            sizes.Add(mu.Length);
        }
        if (covariance is not null)
        {
            if (covariance.GetLength(0) != covariance.GetLength(1))
            {
                throw new ArgumentException(
                    $"covariance must be square, got {covariance.GetLength(0)} x {covariance.GetLength(1)}", nameof(covariance));
            }
            sizes.Add(covariance.GetLength(0));
        }
        if (returns is not null)
        {
            sizes.Add(returns.GetLength(1));
        }
        if (sizes.Count == 0)
        {
            throw new ArgumentException("MarketData needs mu, covariance or returns", nameof(mu));
        }
        if (sizes.Count > 1)
        {
            throw new ArgumentException($"mu, covariance and returns disagree on the asset count: {string.Join(", ", sizes)}",
                                        nameof(mu));
        }
        return sizes.Single();
    }

    private static void CheckSeries(double[,]? returns, double[]? benchmark)
    {
        if (returns is null && benchmark is null)
        {
            return;
        }
        if (returns is null || benchmark is null)
        {
            throw new ArgumentException("returns and benchmark must be given together", nameof(benchmark));
        }
        if (benchmark.Length != returns.GetLength(0))
        {
            throw new ArgumentException(
                $"benchmark has {benchmark.Length} entries, expected {returns.GetLength(0)}", nameof(benchmark));
        }
    }
}
