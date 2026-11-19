using PortOpt.Interop;

namespace PortOpt;

/// <summary>
/// Market data on a native handle: expected returns, covariance, an optional risk-free asset (cash, the last weight)
/// and return series against a benchmark. Every input is copied at construction; the data does not change afterwards.
/// Dispose it when done (<c>using</c>); a call after <see cref="Dispose"/> throws <see cref="ObjectDisposedException"/>.
/// Solves and metrics may run in parallel on one <see cref="MarketData"/>.
/// </summary>
public sealed class MarketData : IDisposable
{
    private readonly DataHandle _handle;

    /// <summary>
    /// Data for n assets, n taken from whichever input is given; the objective decides which it needs
    /// (portopt.h). The covariance must be symmetric positive definite and every value finite.
    /// </summary>
    /// <param name="mu">The expected return of each asset (n entries).</param>
    /// <param name="covariance">The covariance, n x n.</param>
    /// <param name="riskFree">The rate of a risk-free asset, which becomes the last weight.</param>
    /// <param name="returns">Return series for <see cref="Objective.TrackingError"/>, t periods x n assets.</param>
    /// <param name="benchmark">The benchmark return of each period (t entries); given with <paramref name="returns"/>.</param>
    /// <exception cref="ArgumentException">No input, sizes that do not agree, a covariance that is not square, or only
    /// one of <paramref name="returns"/> and <paramref name="benchmark"/>.</exception>
    /// <exception cref="PortOptException">The library rejects a value (non-finite, an asymmetric covariance or one
    /// that is not positive definite).</exception>
    public MarketData(double[]? mu = null, double[,]? covariance = null, double? riskFree = null, double[,]? returns = null,
                      double[]? benchmark = null)
    {
        var assets = Shapes.AssetCount(mu, covariance, returns, benchmark);
        Library.EnsureLoaded();
        _handle = Create(assets, mu, covariance);
        try
        {
            if (riskFree is { } rate)
            {
                Calls.Check(NativeMethods.po_data_set_risk_free(_handle, rate));
            }
            SetReturns(returns, benchmark);
            Calls.Check(NativeMethods.po_data_asset_count(_handle, out var counted));
            Calls.Check(NativeMethods.po_data_weight_count(_handle, out var weights));
            AssetCount = checked((int)counted);
            WeightCount = checked((int)weights);
        }
        catch
        {
            _handle.Dispose();
            throw;
        }
    }

    /// <summary>The number of assets n, cash not counted.</summary>
    public int AssetCount { get; }

    /// <summary>The number of weights: <see cref="AssetCount"/>, plus one with a risk-free asset.</summary>
    public int WeightCount { get; }

    /// <summary>mu' w, plus the rate times the cash weight.</summary>
    /// <param name="weights">One weight per <see cref="WeightCount"/>.</param>
    /// <exception cref="PortOptException">A wrong weight count, or no expected returns.</exception>
    public unsafe double PortfolioReturn(ReadOnlySpan<double> weights)
    {
        return Metric(weights, NativeMethods.po_portfolio_return);
    }

    /// <summary>sqrt(w' Sigma w) over the risky weights.</summary>
    /// <param name="weights">One weight per <see cref="WeightCount"/>.</param>
    /// <exception cref="PortOptException">A wrong weight count, or no covariance.</exception>
    public unsafe double Volatility(ReadOnlySpan<double> weights)
    {
        return Metric(weights, NativeMethods.po_portfolio_volatility);
    }

    /// <summary>The portfolio return less (delta / 2) w' Sigma w.</summary>
    /// <param name="weights">One weight per <see cref="WeightCount"/>.</param>
    /// <param name="riskAversion">delta, positive and finite.</param>
    /// <exception cref="PortOptException">A wrong weight count, a missing input or delta outside its domain.</exception>
    public unsafe double Utility(ReadOnlySpan<double> weights, double riskAversion)
    {
        ObjectDisposedException.ThrowIf(_handle.IsClosed, this);
        fixed (double* w = weights)
        {
            Calls.Check(NativeMethods.po_utility(_handle, w, (nuint)weights.Length, riskAversion, out var value));
            return value;
        }
    }

    /// <summary>The root mean square of the gap between the portfolio and the benchmark over the periods.</summary>
    /// <param name="weights">One weight per <see cref="WeightCount"/>.</param>
    /// <exception cref="PortOptException">A wrong weight count, or no return series.</exception>
    public unsafe double TrackingError(ReadOnlySpan<double> weights)
    {
        return Metric(weights, NativeMethods.po_tracking_error);
    }

    /// <summary>(return - rate) / volatility; NaN when the volatility is 0.</summary>
    /// <param name="weights">One weight per <see cref="WeightCount"/>.</param>
    /// <exception cref="PortOptException">A wrong weight count, or a missing input.</exception>
    public unsafe double SharpeRatio(ReadOnlySpan<double> weights)
    {
        return Metric(weights, NativeMethods.po_sharpe_ratio);
    }

    /// <summary>sigma' w / volatility, sigma the asset volatilities; NaN when the volatility is 0.</summary>
    /// <param name="weights">One weight per <see cref="WeightCount"/>.</param>
    /// <exception cref="PortOptException">A wrong weight count, or no covariance.</exception>
    public unsafe double DiversificationRatio(ReadOnlySpan<double> weights)
    {
        return Metric(weights, NativeMethods.po_diversification_ratio);
    }

    /// <summary>c_i = w_i (Sigma w)_i / w' Sigma w, 0 for the cash entry; all NaN when the volatility is 0.</summary>
    /// <param name="weights">One weight per <see cref="WeightCount"/>.</param>
    /// <returns>One contribution per weight.</returns>
    /// <exception cref="PortOptException">A wrong weight count, or no covariance.</exception>
    public unsafe double[] RiskContributions(ReadOnlySpan<double> weights)
    {
        ObjectDisposedException.ThrowIf(_handle.IsClosed, this);
        var contributions = new double[weights.Length];
        fixed (double* w = weights)
        fixed (double* output = contributions)
        {
            Calls.Check(NativeMethods.po_risk_contributions(_handle, w, (nuint)weights.Length, output));
        }
        return contributions;
    }

    /// <summary>Frees the native data.</summary>
    public void Dispose()
    {
        _handle.Dispose();
    }

    internal DataHandle Handle
    {
        get
        {
            ObjectDisposedException.ThrowIf(_handle.IsClosed, this);
            return _handle;
        }
    }

    private unsafe delegate int Measure(DataHandle data, double* weights, nuint count, out double value);

    private unsafe double Metric(ReadOnlySpan<double> weights, Measure measure)
    {
        ObjectDisposedException.ThrowIf(_handle.IsClosed, this);
        fixed (double* w = weights)
        {
            Calls.Check(measure(_handle, w, (nuint)weights.Length, out var value));
            return value;
        }
    }

    private static unsafe DataHandle Create(int assets, double[]? mu, double[,]? covariance)
    {
        fixed (double* expected = mu)
        fixed (double* matrix = covariance)
        {
            Calls.Check(NativeMethods.po_data_create((nuint)assets, expected, matrix, out var handle));
            return handle;
        }
    }

    private unsafe void SetReturns(double[,]? returns, double[]? benchmark)
    {
        if (returns is null)
        {
            return;
        }
        fixed (double* series = returns)
        fixed (double* reference = benchmark)
        {
            var periods = (nuint)returns.GetLength(0);
            Calls.Check(NativeMethods.po_data_set_returns(_handle, series, reference, periods));
        }
    }
}
