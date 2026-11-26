namespace PortOpt;

// The checks of the input shapes the native side cannot see (it takes pointers and one count), as the Python package
// does before it calls po_data_create.
internal static class Shapes
{
    public static int AssetCount(MarketInputs inputs)
    {
        CheckSeries(inputs);
        var sizes = new HashSet<int>();
        if (inputs.Mu is { } mu)
        {
            sizes.Add(mu.Length);
        }
        if (inputs.Covariance is { } covariance)
        {
            if (covariance.GetLength(0) != covariance.GetLength(1))
            {
                throw new ArgumentException(
                    $"Covariance must be square, got {covariance.GetLength(0)} x {covariance.GetLength(1)}", nameof(inputs));
            }
            sizes.Add(covariance.GetLength(0));
        }
        if (inputs.Returns is { } returns)
        {
            sizes.Add(returns.GetLength(1));
        }
        if (sizes.Count == 0)
        {
            throw new ArgumentException("MarketData needs Mu, Covariance or Returns", nameof(inputs));
        }
        if (sizes.Count > 1)
        {
            throw new ArgumentException($"Mu, Covariance and Returns disagree on the asset count: {string.Join(", ", sizes)}",
                                        nameof(inputs));
        }
        return sizes.Single();
    }

    public static void CheckPick(double[,] p, int assets, string name)
    {
        if (p.GetLength(1) != assets)
        {
            throw new ArgumentException($"P has {p.GetLength(1)} columns, expected {assets} assets", name);
        }
    }

    public static void CheckViews(Views views, int assets)
    {
        ArgumentNullException.ThrowIfNull(views.P, nameof(views));
        ArgumentNullException.ThrowIfNull(views.Q, nameof(views));
        CheckPick(views.P, assets, nameof(views));
        var k = views.P.GetLength(0);
        if (views.Q.Length != k)
        {
            throw new ArgumentException($"Q has {views.Q.Length} entries, expected {k}", nameof(views));
        }
        if (views.Omega is null && k > 0)
        {
            throw new ArgumentException("Omega is required with at least one view", nameof(views));
        }
        if (views.Omega is { } omega && (omega.GetLength(0) != k || omega.GetLength(1) != k))
        {
            throw new ArgumentException(
                $"Omega is {omega.GetLength(0)} x {omega.GetLength(1)}, expected {k} x {k}", nameof(views));
        }
    }

    private static void CheckSeries(MarketInputs inputs)
    {
        var (returns, benchmark) = (inputs.Returns, inputs.Benchmark);
        if (returns is null && benchmark is null)
        {
            return;
        }
        if (returns is null || benchmark is null)
        {
            throw new ArgumentException("Returns and Benchmark must be given together", nameof(inputs));
        }
        if (benchmark.Length != returns.GetLength(0))
        {
            throw new ArgumentException(
                $"Benchmark has {benchmark.Length} entries, expected {returns.GetLength(0)}", nameof(inputs));
        }
    }
}
