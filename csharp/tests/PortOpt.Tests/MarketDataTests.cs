using PortOpt.Tests.Support;

namespace PortOpt.Tests;

public sealed class MarketDataTests
{
    private static readonly double[] Equal = [1.0 / 3, 1.0 / 3, 1.0 / 3];

    [Fact]
    public void Counts_follow_the_inputs_and_the_cash()
    {
        using var risky = Samples.Data();
        using var withCash = new MarketData(new MarketInputs { Mu = Samples.Mu, Covariance = Samples.Covariance, RiskFree = 0.02 });
        Assert.Equal((3, 3), (risky.AssetCount, risky.WeightCount));
        Assert.Equal((3, 4), (withCash.AssetCount, withCash.WeightCount));
    }

    [Fact]
    public void Metrics_are_the_formulas_of_the_header()
    {
        using var data = Samples.Data();
        var covariance = Samples.Covariance;
        var variance = 0.0;
        for (var i = 0; i < 3; i++)
        {
            for (var j = 0; j < 3; j++)
            {
                variance += Equal[i] * covariance[i, j] * Equal[j];
            }
        }
        Assert.Equal(0.1, data.PortfolioReturn(Equal), 15);
        Assert.Equal(Math.Sqrt(variance), data.Volatility(Equal), 15);
        Assert.Equal(0.1 - (2.5 * variance), data.Utility(Equal, 5), 15);
        Assert.Equal(0.1 / Math.Sqrt(variance), data.SharpeRatio(Equal), 13);
        Assert.Equal((0.2 + 0.3 + Math.Sqrt(0.06)) / 3 / Math.Sqrt(variance), data.DiversificationRatio(Equal), 13);
        Assert.Equal(1.0, data.RiskContributions(Equal).Sum(), 15);
    }

    [Fact]
    public void Tracking_error_of_the_benchmark_itself_is_zero()
    {
        double[,] returns = { { 0.01, 0.02 }, { -0.01, 0.03 }, { 0.02, -0.02 } };
        double[] benchmark = [0.01, -0.01, 0.02];
        using var data = new MarketData(new MarketInputs { Returns = returns, Benchmark = benchmark });
        Assert.Equal(0.0, data.TrackingError([1.0, 0.0]));
        Assert.True(data.TrackingError([0.0, 1.0]) > 0);
    }

    [Fact]
    public void A_wrong_weight_count_is_a_native_error()
    {
        using var data = Samples.Data();
        var error = Assert.Throws<PortOptException>(() => data.Volatility([0.5, 0.5]));
        Assert.Equal(Status.InvalidArgument, error.Status);
    }

    [Fact]
    public void Shapes_the_native_side_cannot_see_are_argument_errors()
    {
        Assert.Throws<ArgumentNullException>(() => new MarketData(null!));
        Assert.Throws<ArgumentException>(() => new MarketData(new MarketInputs()));
        Assert.Throws<ArgumentException>(() => new MarketData(new MarketInputs { Mu = [0.1, 0.2], Covariance = Samples.Covariance }));
        Assert.Throws<ArgumentException>(() => new MarketData(new MarketInputs { Covariance = new double[2, 3] }));
        Assert.Throws<ArgumentException>(() => new MarketData(new MarketInputs { Returns = new double[4, 3] }));
        Assert.Throws<ArgumentException>(() => new MarketData(new MarketInputs { Returns = new double[4, 3], Benchmark = new double[3] }));
    }

    [Fact]
    public void A_covariance_that_is_not_positive_definite_is_refused()
    {
        var error = Assert.Throws<PortOptException>(() => new MarketData(new MarketInputs { Covariance = new[,] { { 1.0, 2.0 }, { 2.0, 1.0 } } }));
        Assert.Equal(Status.InvalidArgument, error.Status);
    }

    [Fact]
    public void A_disposed_data_refuses_every_call()
    {
        var data = Samples.Data();
        data.Dispose();
        data.Dispose();
        Assert.Throws<ObjectDisposedException>(() => data.Volatility(Equal));
        Assert.Throws<ObjectDisposedException>(() => data.RiskContributions(Equal));
        Assert.Throws<ObjectDisposedException>(() => Optimizer.Solve(new Spec(), data));
    }
}
