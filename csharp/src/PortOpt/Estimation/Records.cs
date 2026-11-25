namespace PortOpt;

/// <summary>The sample mean and covariance of a return panel.</summary>
/// <param name="Mu">The mean of each column (n entries).</param>
/// <param name="Covariance">The covariance with divisor t - 1, n x n, exactly symmetric.</param>
public sealed record Moments(double[] Mu, double[,] Covariance);

/// <summary>The Bayes-Stein (Jorion 1986) shrinkage of a mean.</summary>
/// <param name="Mu">The shrunk mean (n entries).</param>
/// <param name="GrandMean">The grand mean mu_0 the mean is shrunk toward.</param>
/// <param name="Weight">The weight w of the grand mean.</param>
public sealed record BayesSteinEstimate(double[] Mu, double GrandMean, double Weight);

/// <summary>The k views of a Black-Litterman model, in the unit of the equilibrium return (excess of the rate).</summary>
/// <param name="P">The pick matrix, k x n, of full row rank.</param>
/// <param name="Q">The view returns (k entries).</param>
/// <param name="Omega">The uncertainty of the views, k x k; it can be null only with no view (k = 0).</param>
public sealed record Views(double[,] P, double[] Q, double[,]? Omega);

/// <summary>The Black-Litterman posterior, in excess of the risk-free rate.</summary>
/// <param name="Mu">The posterior mean (n entries).</param>
/// <param name="Covariance">The posterior covariance, n x n.</param>
public sealed record BlackLittermanPosterior(double[] Mu, double[,] Covariance)
{
    /// <summary>
    /// The inputs of an optimization. The posterior mean is in excess of the rate, so with a rate it is added to every
    /// entry and the rate becomes the risk-free asset (cash); without one the mean stays as it is and there is no cash.
    /// The covariance is the posterior one. The arrays are copied.
    /// </summary>
    /// <param name="riskFree">The risk-free rate, if the optimization has cash.</param>
    /// <returns>The inputs of a <see cref="MarketData"/>.</returns>
    public MarketInputs ToMarketInputs(double? riskFree = null)
    {
        var rate = riskFree ?? 0.0;
        return new MarketInputs
        {
            Mu = [.. Mu.Select(value => value + rate)],
            Covariance = (double[,])Covariance.Clone(),
            RiskFree = riskFree,
        };
    }
}
