namespace PortOpt.Tests;

// The cases of python/tests/test_estimation.py; the references are the values numpy computes there (np.mean, np.cov
// with ddof=1, np.linalg.solve), copied as literals.
public sealed class EstimationTests
{
    private const double Tolerance = 1e-14;
    private const double Tau = 0.05;

    private static readonly double[,] Covariance = new[,]
    {
        { 0.04, 0.012, 0.006 },
        { 0.012, 0.09, 0.018 },
        { 0.006, 0.018, 0.0225 },
    };

    private static readonly double[] MarketWeights = [0.5, 0.3, 0.2];
    private static readonly double[] Pi = [0.062000000000000006, 0.09149999999999998, 0.032249999999999994];
    private static readonly double[,] RelativeP = { { 1.0, -1.0, 0.0 } };
    private static readonly double[,] TwoViewsP = { { 1.0, -1.0, 0.0 }, { 0.0, 1.0, 1.0 } };

    [Fact]
    public void SampleMoments_match_numpy()
    {
        double[,] returns =
        {
            { 0.01, 0.02, -0.01 }, { 0.03, -0.01, 0.02 }, { -0.02, 0.015, 0.01 }, { 0.005, 0.0, -0.005 }, { 0.02, 0.01, 0.03 },
        };
        var moments = Estimation.SampleMoments(returns);
        Close([0.009, 0.007000000000000001, 0.009], moments.Mu);
        Close(
            [0.00035499999999999996, -0.00012874999999999996, 0.00011749999999999998,
             -0.00012874999999999996, 0.000145, -6.625000000000001e-05,
             0.00011749999999999998, -6.625000000000001e-05, 0.00028],
            moments.Covariance);
    }

    [Fact]
    public void BayesStein_matches_numpy()
    {
        var covariance = (double[,])Covariance.Clone();
        for (var i = 0; i < 9; i++)
        {
            covariance[i / 3, i % 3] /= 100;
        }
        var estimate = Estimation.BayesStein([0.001, 0.002, 0.0005], covariance, 50);
        Close([0.0006922637739080789, 0.0007200970596631461, 0.0006783471310305453], estimate.Mu);
        Assert.InRange(Math.Abs(estimate.GrandMean - 0.0006834532374100721), 0, Tolerance);
        Assert.InRange(Math.Abs(estimate.Weight - 0.9721667142449328), 0, Tolerance);
    }

    [Fact]
    public void Equilibrium_is_delta_sigma_w()
    {
        Close(Pi, Estimation.Equilibrium(Covariance, MarketWeights, 2.5));
    }

    [Fact]
    public void Omegas_are_diagonal()
    {
        Close([0.0053], Estimation.OmegaHeLitterman(Covariance, RelativeP, Tau));
        Close([0.0053, 0.0, 0.0, 0.0], Estimation.OmegaIdzorek(Covariance, TwoViewsP, Tau, [0.5, 1.0]));
    }

    [Fact]
    public void Omegas_without_views_are_empty()
    {
        var none = new double[0, 3];
        Assert.Equal(new double[0, 0], Estimation.OmegaHeLitterman(Covariance, none, Tau));
        Assert.Equal(new double[0, 0], Estimation.OmegaIdzorek(Covariance, none, Tau, []));
    }

    [Fact]
    public void BayesStein_rejects_a_negative_period_count()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Estimation.BayesStein([0.001, 0.002, 0.0005], Covariance, -1));
    }

    [Fact]
    public void Posterior_matches_numpy_for_a_relative_view()
    {
        var omega = Estimation.OmegaHeLitterman(Covariance, RelativeP, Tau);
        var posterior = Estimation.Posterior(Covariance, Pi, Tau, new Views(RelativeP, [0.02], omega));
        Close([0.0685377358490566, 0.0732877358490566, 0.029448113207547165], posterior.Mu);
        Close(
            [0.04181509433962265, 0.013115094339622642, 0.006379245283018868,
             0.013115094339622642, 0.09306509433962264, 0.01867924528301887,
             0.006379245283018868, 0.01867924528301887, 0.023591037735849057],
            posterior.Covariance);
    }

    [Fact]
    public void Posterior_reproduces_a_full_confidence_view()
    {
        double[] q = [0.02, 0.05];
        var omega = Estimation.OmegaIdzorek(Covariance, TwoViewsP, Tau, [1.0, 1.0]);
        var mu = Estimation.Posterior(Covariance, Pi, Tau, new Views(TwoViewsP, q, omega)).Mu;
        Close(q, [mu[0] - mu[1], mu[1] + mu[2]]);
    }

    [Fact]
    public void No_views_return_the_prior()
    {
        var posterior = Estimation.Posterior(Covariance, Pi, Tau, new Views(new double[0, 3], [], null));
        Assert.Equal(Pi, posterior.Mu);
        Close(Covariance.Cast<double>().Select(value => value * 1.05).ToArray(), posterior.Covariance);
    }

    [Fact]
    public void ToMarketInputs_adds_the_rate_and_sets_cash()
    {
        var posterior = new BlackLittermanPosterior([0.01, 0.02, 0.03], Covariance);
        var withRate = posterior.ToMarketInputs(0.02);
        Assert.Equal(0.02, withRate.RiskFree);
        for (var i = 0; i < 3; i++)
        {
            Assert.Equal(posterior.Mu[i] + 0.02, withRate.Mu![i]);
        }
        Assert.NotSame(posterior.Mu, withRate.Mu);
        Assert.NotSame(posterior.Covariance, withRate.Covariance);
        using var data = new MarketData(withRate);
        Assert.Equal(4, data.WeightCount);
        var without = posterior.ToMarketInputs();
        Assert.Null(without.RiskFree);
        Assert.Equal(posterior.Mu, without.Mu);
        Assert.Equal(posterior.Covariance, without.Covariance);
    }

    [Fact]
    public void Redundant_views_throw_PortOptException()
    {
        double[,] p = { { 1.0, -1.0, 0.0 }, { 2.0, -2.0, 0.0 } };
        var error = Assert.Throws<PortOptException>(() => Estimation.OmegaHeLitterman(Covariance, p, Tau));
        Assert.Equal(Status.InvalidArgument, error.Status);
    }

    [Fact]
    public void Views_with_a_wrong_shape_throw_ArgumentException()
    {
        double[,] omega = { { 0.0 } };
        Assert.Throws<ArgumentException>(() => Estimation.Posterior(Covariance, Pi, Tau, new Views(RelativeP, [0.02, 0.03], omega)));
        Assert.Throws<ArgumentException>(() => Estimation.Posterior(Covariance, Pi, Tau, new Views(new double[1, 2], [0.02], omega)));
        Assert.Throws<ArgumentException>(() => Estimation.Posterior(Covariance, Pi, Tau, new Views(RelativeP, [0.02], new double[2, 2])));
        Assert.Throws<ArgumentException>(() => Estimation.Posterior(Covariance, Pi, Tau, new Views(RelativeP, [0.02], null)));
        Assert.Throws<ArgumentException>(() => Estimation.Equilibrium(Covariance, [0.5, 0.5], 2.5));
        Assert.Throws<ArgumentException>(() => Estimation.OmegaIdzorek(Covariance, RelativeP, Tau, [0.5, 0.5]));
        Assert.Throws<ArgumentException>(() => Estimation.OmegaIdzorek(Covariance, new double[1, 2], Tau, [0.5]));
        Assert.Throws<ArgumentException>(() => Estimation.OmegaHeLitterman(Covariance, new double[1, 4], Tau));
        Assert.Throws<ArgumentException>(() => Estimation.Posterior(Covariance, [0.01, 0.02], Tau, new Views(RelativeP, [0.02], omega)));
    }

    [Fact]
    public void BayesStein_on_MarketData_is_the_bits_of_the_arrays()
    {
        var covariance = (double[,])Covariance.Clone();
        for (var i = 0; i < 9; i++)
        {
            covariance[i / 3, i % 3] /= 100;
        }
        double[] muHat = [0.001, 0.002, 0.0005];
        var expected = Estimation.BayesStein(muHat, covariance, 50);
        using var data = new MarketData(new MarketInputs { Mu = muHat, Covariance = covariance });
        var estimate = Estimation.BayesStein(data, 50);
        Assert.Equal(expected.Mu, estimate.Mu);
        Assert.Equal(expected.GrandMean, estimate.GrandMean);
        Assert.Equal(expected.Weight, estimate.Weight);
        Assert.Throws<ArgumentOutOfRangeException>(() => Estimation.BayesStein(data, -1));
    }

    [Fact]
    public void One_MarketData_gives_the_bits_of_the_covariance()
    {
        double[] q = [0.02, 0.05];
        // The equilibrium and the Omegas do not read the mean, so one data with mean Pi serves the whole model.
        using var data = new MarketData(new MarketInputs { Mu = Pi, Covariance = Covariance });
        Assert.Equal(Estimation.Equilibrium(Covariance, MarketWeights, 2.5), Estimation.Equilibrium(data, MarketWeights, 2.5));
        var omega = Estimation.OmegaHeLitterman(Covariance, TwoViewsP, Tau);
        Assert.Equal(omega, Estimation.OmegaHeLitterman(data, TwoViewsP, Tau));
        Assert.Equal(
            Estimation.OmegaIdzorek(Covariance, TwoViewsP, Tau, [0.5, 1.0]),
            Estimation.OmegaIdzorek(data, TwoViewsP, Tau, [0.5, 1.0]));
        var views = new Views(TwoViewsP, q, omega);
        var expected = Estimation.Posterior(Covariance, Pi, Tau, views);
        var posterior = Estimation.Posterior(data, Tau, views);
        Assert.Equal(expected.Mu, posterior.Mu);
        Assert.Equal(expected.Covariance, posterior.Covariance);
        Assert.Equal(3, data.AssetCount);
    }

    [Fact]
    public void MarketData_overloads_check_the_data()
    {
        double[,] omega = { { 0.0 } };
        using var covarianceOnly = new MarketData(new MarketInputs { Covariance = Covariance });
        Assert.Throws<ArgumentException>(() => Estimation.Equilibrium(covarianceOnly, [0.5, 0.5], 2.5));
        Assert.Throws<ArgumentException>(() => Estimation.OmegaHeLitterman(covarianceOnly, new double[1, 4], Tau));
        Assert.Throws<ArgumentException>(() => Estimation.OmegaIdzorek(covarianceOnly, RelativeP, Tau, [0.5, 0.5]));
        Assert.Throws<ArgumentException>(() => Estimation.Posterior(covarianceOnly, Tau, new Views(new double[1, 2], [0.02], omega)));
        var noMean = Assert.Throws<PortOptException>(() => Estimation.Posterior(covarianceOnly, Tau, new Views(RelativeP, [0.02], omega)));
        Assert.Equal(Status.InvalidArgument, noMean.Status);
        Assert.Throws<PortOptException>(() => Estimation.BayesStein(covarianceOnly, 50));
        using var withRate = new MarketData(new MarketInputs { Mu = Pi, Covariance = Covariance, RiskFree = 0.02 });
        Assert.Throws<PortOptException>(() => Estimation.Equilibrium(withRate, MarketWeights, 2.5));
        Assert.Throws<PortOptException>(() => Estimation.Posterior(withRate, Tau, new Views(RelativeP, [0.02], omega)));
        Assert.Throws<ArgumentNullException>(() => Estimation.Equilibrium((MarketData)null!, MarketWeights, 2.5));
        var disposed = new MarketData(new MarketInputs { Covariance = Covariance });
        disposed.Dispose();
        Assert.Throws<ObjectDisposedException>(() => Estimation.Equilibrium(disposed, MarketWeights, 2.5));
        Assert.Throws<ObjectDisposedException>(() => Estimation.OmegaHeLitterman(disposed, RelativeP, Tau));
    }

    private static void Close(double[] expected, double[] actual)
    {
        Assert.Equal(expected.Length, actual.Length);
        for (var i = 0; i < expected.Length; i++)
        {
            Assert.InRange(Math.Abs(expected[i] - actual[i]), 0, Tolerance);
        }
    }

    private static void Close(double[] expected, double[,] actual)
    {
        Close(expected, actual.Cast<double>().ToArray());
    }
}
