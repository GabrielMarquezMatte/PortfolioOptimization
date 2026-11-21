namespace PortOpt;

/// <summary>
/// The inputs of a <see cref="MarketData"/>, n assets taken from whichever is given; the objective decides which it
/// needs (portopt.h). <see cref="MarketData"/> copies them, so changing an array afterwards changes nothing. A later
/// version adds an input as a new property.
/// </summary>
public sealed class MarketInputs
{
    /// <summary>The expected return of each asset (n entries).</summary>
    public double[]? Mu { get; init; }

    /// <summary>The covariance, n x n, symmetric positive definite.</summary>
    public double[,]? Covariance { get; init; }

    /// <summary>The rate of a risk-free asset, which becomes the last weight.</summary>
    public double? RiskFree { get; init; }

    /// <summary>Return series for <see cref="Objective.TrackingError"/>, t periods x n assets.</summary>
    public double[,]? Returns { get; init; }

    /// <summary>The benchmark return of each period (t entries); given with <see cref="Returns"/>.</summary>
    public double[]? Benchmark { get; init; }
}
