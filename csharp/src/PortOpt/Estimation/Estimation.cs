using PortOpt.Interop;

namespace PortOpt;

/// <summary>
/// The estimation utilities: sample moments, the Bayes-Stein shrinkage and the Black-Litterman model. Pure functions on
/// arrays (row-major <c>double[,]</c> matrices); each call builds its native data and frees it before it returns. Every
/// function that takes a covariance also has an overload on a <see cref="MarketData"/>, which reuses the data instead
/// of copying, validating and factoring the covariance on every call (the factor costs n^3 / 3); that data must have no
/// risk-free rate, and its return series are not used.
/// </summary>
public static class Estimation
{
    /// <summary>The mean and the covariance (divisor t - 1) of a t x n return panel; needs t of at least 2.</summary>
    /// <param name="returns">The returns, one row per period.</param>
    /// <returns>The moments; the covariance is exactly symmetric and not checked for positive definiteness.</returns>
    /// <exception cref="PortOptException">Fewer than two periods, no asset or a non-finite value.</exception>
    public static unsafe Moments SampleMoments(double[,] returns)
    {
        ArgumentNullException.ThrowIfNull(returns);
        Library.EnsureLoaded();
        var assets = returns.GetLength(1);
        var (mu, covariance) = (new double[assets], new double[assets, assets]);
        fixed (double* panel = returns)
        fixed (double* outMu = mu)
        fixed (double* outCov = covariance)
        {
            var periods = (nuint)returns.GetLength(0);
            Calls.Check(NativeMethods.po_sample_moments(periods, (nuint)assets, panel, outMu, outCov));
        }
        return new Moments(mu, covariance);
    }

    /// <summary>
    /// The Bayes-Stein shrinkage of the mean (Jorion 1986). The weight depends on the frequency of the data: estimate at
    /// the original frequency and annualize the resulting mean afterwards.
    /// </summary>
    /// <param name="muHat">The sample mean (n entries).</param>
    /// <param name="covariance">The sample covariance with divisor t - 1, n x n.</param>
    /// <param name="periods">The number of periods t, more than n + 2.</param>
    /// <returns>The shrunk mean, the grand mean and the weight.</returns>
    /// <exception cref="ArgumentException">The sizes do not agree.</exception>
    /// <exception cref="PortOptException">t is not above n + 2, or the covariance is rejected.</exception>
    public static BayesSteinEstimate BayesStein(double[] muHat, double[,] covariance, int periods)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(periods);
        using var data = new MarketData(new MarketInputs { Mu = muHat, Covariance = covariance });
        return BayesStein(data, periods);
    }

    /// <summary>
    /// The Bayes-Stein shrinkage of the mean of <paramref name="data"/>, whose covariance is the sample covariance with
    /// divisor t - 1: <see cref="BayesStein(double[], double[,], int)"/> on the mean and the covariance of the data.
    /// </summary>
    /// <param name="data">The sample mean and covariance, with no risk-free rate; it is not disposed.</param>
    /// <param name="periods">The number of periods t, more than n + 2.</param>
    /// <returns>The shrunk mean, the grand mean and the weight.</returns>
    /// <exception cref="ObjectDisposedException">The data is disposed.</exception>
    /// <exception cref="PortOptException">t is not above n + 2, or the data has no mean or has a risk-free
    /// rate.</exception>
    public static unsafe BayesSteinEstimate BayesStein(MarketData data, int periods)
    {
        ArgumentNullException.ThrowIfNull(data);
        ArgumentOutOfRangeException.ThrowIfNegative(periods);
        var mu = new double[data.AssetCount];
        fixed (double* output = mu)
        {
            Calls.Check(NativeMethods.po_bayes_stein(data.Handle, (nuint)periods, output, out var grandMean, out var weight));
            return new BayesSteinEstimate(mu, grandMean, weight);
        }
    }

    /// <summary>The equilibrium returns Pi = delta Sigma w, in excess of the risk-free rate.</summary>
    /// <param name="covariance">Sigma, n x n.</param>
    /// <param name="marketWeights">The market weights (n entries, summing to 1).</param>
    /// <param name="delta">The risk aversion, positive.</param>
    /// <returns>Pi (n entries).</returns>
    /// <exception cref="ArgumentException">The sizes do not agree.</exception>
    /// <exception cref="PortOptException">Weights that do not sum to 1, or a delta outside its domain.</exception>
    public static double[] Equilibrium(double[,] covariance, double[] marketWeights, double delta)
    {
        ArgumentNullException.ThrowIfNull(marketWeights);
        using var data = new MarketData(new MarketInputs { Covariance = covariance });
        return Equilibrium(data, marketWeights, delta);
    }

    /// <summary>The equilibrium returns Pi = delta Sigma w, Sigma the covariance of <paramref name="data"/>.</summary>
    /// <param name="data">Sigma, with no risk-free rate; its mean, if any, is not used, and it is not disposed.</param>
    /// <param name="marketWeights">The market weights (n entries, summing to 1).</param>
    /// <param name="delta">The risk aversion, positive.</param>
    /// <returns>Pi (n entries).</returns>
    /// <exception cref="ArgumentException">The sizes do not agree.</exception>
    /// <exception cref="ObjectDisposedException">The data is disposed.</exception>
    /// <exception cref="PortOptException">Weights that do not sum to 1, a delta outside its domain, or data with no
    /// covariance or with a risk-free rate.</exception>
    public static unsafe double[] Equilibrium(MarketData data, double[] marketWeights, double delta)
    {
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(marketWeights);
        if (marketWeights.Length != data.AssetCount)
        {
            throw new ArgumentException(
                $"marketWeights has {marketWeights.Length} entries, expected {data.AssetCount}", nameof(marketWeights));
        }
        var pi = new double[data.AssetCount];
        fixed (double* weights = marketWeights)
        fixed (double* output = pi)
        {
            Calls.Check(NativeMethods.po_bl_equilibrium(data.Handle, weights, (nuint)pi.Length, delta, output));
        }
        return pi;
    }

    /// <summary>The He and Litterman Omega, diagonal with omega_kk = tau p_k Sigma p_k'.</summary>
    /// <param name="covariance">Sigma, n x n.</param>
    /// <param name="p">The pick matrix, k x n.</param>
    /// <param name="tau">The scale of the prior covariance, positive.</param>
    /// <returns>Omega, k x k.</returns>
    /// <exception cref="ArgumentException">The sizes do not agree.</exception>
    /// <exception cref="PortOptException">Redundant views, or a tau outside its domain.</exception>
    public static double[,] OmegaHeLitterman(double[,] covariance, double[,] p, double tau)
    {
        using var data = new MarketData(new MarketInputs { Covariance = covariance });
        return OmegaHeLitterman(data, p, tau);
    }

    /// <summary>The He and Litterman Omega, with Sigma the covariance of <paramref name="data"/>.</summary>
    /// <param name="data">Sigma, with no risk-free rate; its mean, if any, is not used, and it is not disposed.</param>
    /// <param name="p">The pick matrix, k x n.</param>
    /// <param name="tau">The scale of the prior covariance, positive.</param>
    /// <returns>Omega, k x k.</returns>
    /// <exception cref="ArgumentException">The sizes do not agree.</exception>
    /// <exception cref="ObjectDisposedException">The data is disposed.</exception>
    /// <exception cref="PortOptException">Redundant views, a tau outside its domain, or data with no covariance or with
    /// a risk-free rate.</exception>
    public static unsafe double[,] OmegaHeLitterman(MarketData data, double[,] p, double tau)
    {
        ArgumentNullException.ThrowIfNull(data);
        using var views = CreateViews(p, data.AssetCount, null);
        var omega = new double[p.GetLength(0), p.GetLength(0)];
        fixed (double* output = omega)
        {
            Calls.Check(NativeMethods.po_bl_omega_he_litterman(data.Handle, views, tau, output));
        }
        return omega;
    }

    /// <summary>The Idzorek Omega, diagonal with omega_kk = ((1 - c_k) / c_k) tau p_k Sigma p_k'; c_k = 1 gives 0.</summary>
    /// <param name="covariance">Sigma, n x n.</param>
    /// <param name="p">The pick matrix, k x n.</param>
    /// <param name="tau">The scale of the prior covariance, positive.</param>
    /// <param name="confidence">The confidence of each view (k entries), above 0 and at most 1.</param>
    /// <returns>Omega, k x k.</returns>
    /// <exception cref="ArgumentException">The sizes do not agree.</exception>
    /// <exception cref="PortOptException">Redundant views, or a value outside its domain.</exception>
    public static double[,] OmegaIdzorek(double[,] covariance, double[,] p, double tau, double[] confidence)
    {
        ArgumentNullException.ThrowIfNull(confidence);
        using var data = new MarketData(new MarketInputs { Covariance = covariance });
        return OmegaIdzorek(data, p, tau, confidence);
    }

    /// <summary>The Idzorek Omega, with Sigma the covariance of <paramref name="data"/>.</summary>
    /// <param name="data">Sigma, with no risk-free rate; its mean, if any, is not used, and it is not disposed.</param>
    /// <param name="p">The pick matrix, k x n.</param>
    /// <param name="tau">The scale of the prior covariance, positive.</param>
    /// <param name="confidence">The confidence of each view (k entries), above 0 and at most 1.</param>
    /// <returns>Omega, k x k.</returns>
    /// <exception cref="ArgumentException">The sizes do not agree.</exception>
    /// <exception cref="ObjectDisposedException">The data is disposed.</exception>
    /// <exception cref="PortOptException">Redundant views, a value outside its domain, or data with no covariance or
    /// with a risk-free rate.</exception>
    public static unsafe double[,] OmegaIdzorek(MarketData data, double[,] p, double tau, double[] confidence)
    {
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(confidence);
        using var views = CreateViews(p, data.AssetCount, null);
        var k = p.GetLength(0);
        if (confidence.Length != k)
        {
            throw new ArgumentException($"confidence has {confidence.Length} entries, expected {k}", nameof(confidence));
        }
        var omega = new double[k, k];
        fixed (double* levels = confidence)
        fixed (double* output = omega)
        {
            Calls.Check(NativeMethods.po_bl_omega_idzorek(data.Handle, views, tau, levels, output));
        }
        return omega;
    }

    /// <summary>
    /// The Black-Litterman posterior. It is in excess of the risk-free rate: use
    /// <see cref="BlackLittermanPosterior.ToMarketInputs"/> to add the rate before an optimization.
    /// </summary>
    /// <param name="covariance">Sigma, n x n.</param>
    /// <param name="pi">The equilibrium returns Pi (n entries), in excess of the rate.</param>
    /// <param name="tau">The scale of the prior covariance, positive.</param>
    /// <param name="views">The views; with none, the posterior is Pi and (1 + tau) Sigma.</param>
    /// <returns>The posterior mean and covariance.</returns>
    /// <exception cref="ArgumentException">The sizes do not agree.</exception>
    /// <exception cref="PortOptException">Redundant views, an Omega that makes the solve fail, or a tau outside its domain.</exception>
    public static BlackLittermanPosterior Posterior(double[,] covariance, double[] pi, double tau, Views views)
    {
        ArgumentNullException.ThrowIfNull(views);
        using var data = new MarketData(new MarketInputs { Mu = pi, Covariance = covariance });
        return Posterior(data, tau, views);
    }

    /// <summary>
    /// The Black-Litterman posterior with Pi the mean of <paramref name="data"/> and Sigma its covariance:
    /// <see cref="Posterior(double[,], double[], double, Views)"/> on them.
    /// </summary>
    /// <param name="data">Pi (in excess of the rate) and Sigma, with no risk-free rate; it is not disposed.</param>
    /// <param name="tau">The scale of the prior covariance, positive.</param>
    /// <param name="views">The views; with none, the posterior is Pi and (1 + tau) Sigma.</param>
    /// <returns>The posterior mean and covariance.</returns>
    /// <exception cref="ArgumentException">The sizes do not agree.</exception>
    /// <exception cref="ObjectDisposedException">The data is disposed.</exception>
    /// <exception cref="PortOptException">Redundant views, an Omega that makes the solve fail, a tau outside its domain,
    /// or data with no mean, no covariance or a risk-free rate.</exception>
    public static unsafe BlackLittermanPosterior Posterior(MarketData data, double tau, Views views)
    {
        ArgumentNullException.ThrowIfNull(data);
        ArgumentNullException.ThrowIfNull(views);
        Shapes.CheckViews(views, data.AssetCount);
        using var handle = CreateViews(views.P, data.AssetCount, views.Q);
        if (views.Omega is { Length: > 0 } omega)
        {
            SetOmega(handle, omega);
        }
        var (mu, posterior) = (new double[data.AssetCount], new double[data.AssetCount, data.AssetCount]);
        fixed (double* outMu = mu)
        fixed (double* outCov = posterior)
        {
            Calls.Check(NativeMethods.po_bl_posterior(data.Handle, handle, tau, outMu, outCov));
        }
        return new BlackLittermanPosterior(mu, posterior);
    }

    // The Omega helpers run before there are view returns, so their q is zeros.
    private static unsafe ViewsHandle CreateViews(double[,] p, int assets, double[]? q)
    {
        ArgumentNullException.ThrowIfNull(p);
        Shapes.CheckPick(p, assets, nameof(p));
        var k = p.GetLength(0);
        fixed (double* pick = p)
        fixed (double* returns = q ?? new double[k])
        {
            Calls.Check(NativeMethods.po_views_create((nuint)k, (nuint)assets, pick, returns, out var handle));
            return handle;
        }
    }

    private static unsafe void SetOmega(ViewsHandle views, double[,] omega)
    {
        fixed (double* matrix = omega)
        {
            Calls.Check(NativeMethods.po_views_set_omega(views, matrix, (nuint)omega.GetLength(0)));
        }
    }
}
