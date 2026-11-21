using PortOpt.Tests.Support;

namespace PortOpt.Tests;

public sealed class SolveTests
{
    public static TheoryData<Spec> EveryObjective => new()
    {
        new Spec { Objective = Objective.MinVariance, Bounds = new Bounds(0, 1) },
        new Spec { Objective = Objective.MaxReturn, Bounds = new Bounds(0, 1), MaxVolatility = 0.22 },
        new Spec { Objective = Objective.MaxSharpe, Bounds = new Bounds(0, 1) },
        new Spec { Objective = Objective.Utility, Bounds = new Bounds(0, 1), RiskAversion = 3 },
        new Spec { Objective = Objective.Importance, Bounds = new Bounds(0, 1), Importance = 0.3 },
        new Spec { Objective = Objective.Kelly, Bounds = new Bounds(-1, 2), Leverage = 2 },
        new Spec { Objective = Objective.MaxDiversification, Bounds = new Bounds(0, 1) },
        new Spec { Objective = Objective.RiskParity, RiskBudget = [0.5, 0.25, 0.25] },
    };

    [Theory]
    [MemberData(nameof(EveryObjective))]
    public void Every_objective_solves_to_weights_that_sum_to_one(Spec spec)
    {
        using var data = Samples.Data();
        var result = Optimizer.Solve(spec, data);
        Assert.Equal(SolveStatus.Optimal, result.Status);
        Assert.Equal(3, result.Weights.Length);
        Assert.Equal(1.0, result.Weights.Sum(), 12);
        Assert.Equal(3, result.BoundMultipliers.Length);
    }

    [Fact]
    public void Minimum_variance_matches_the_c_interface_bit_for_bit()
    {
        using var data = Samples.Data();
        var result = Optimizer.Solve(new Spec { Bounds = new Bounds(0, 1) }, data);
        var (weights, multipliers) = SolveThroughC();
        Assert.Equal(BitsOf(weights), BitsOf(result.Weights));
        double[] named = [result.BudgetMultiplier, result.TargetReturnMultiplier, result.VolatilityMultiplier, result.LeverageMultiplier];
        Assert.Equal(BitsOf(multipliers), BitsOf(named));
    }

    [Fact]
    public void Risk_parity_meets_its_budget()
    {
        using var data = Samples.Data();
        double[] budget = [0.5, 0.25, 0.25];
        var result = Optimizer.Solve(new Spec { Objective = Objective.RiskParity, RiskBudget = budget }, data);
        var contributions = data.RiskContributions(result.Weights);
        for (var i = 0; i < budget.Length; i++)
        {
            Assert.Equal(budget[i], contributions[i], 10);
        }
        Assert.True(result.Residual < 1e-20);
    }

    [Fact]
    public void A_target_return_is_met_and_has_a_multiplier()
    {
        using var data = Samples.Data();
        var result = Optimizer.Solve(new Spec { Bounds = new Bounds(0, 1), TargetReturn = 0.11 }, data);
        Assert.Equal(0.11, data.PortfolioReturn(result.Weights), 12);
        Assert.False(double.IsNaN(result.TargetReturnMultiplier));
        Assert.True(double.IsNaN(result.LeverageMultiplier));
    }

    [Fact]
    public void Tracking_error_with_cash_puts_the_cash_weight_last()
    {
        double[,] returns = { { 0.01, 0.02 }, { -0.01, 0.03 }, { 0.02, -0.02 }, { 0.00, 0.01 } };
        double[] benchmark = [0.012, 0.0, 0.01, 0.004];
        using var data = new MarketData(new MarketInputs { Returns = returns, Benchmark = benchmark, RiskFree = 0.001 });
        var result = Optimizer.Solve(new Spec { Objective = Objective.TrackingError, Bounds = new Bounds(0, 1) }, data);
        Assert.Equal(SolveStatus.Optimal, result.Status);
        Assert.Equal(3, result.Weights.Length);
    }

    [Fact]
    public void An_infeasible_problem_is_a_status_with_nan_weights()
    {
        using var data = Samples.Data();
        var result = Optimizer.Solve(new Spec { Bounds = new Bounds(0, 0.2) }, data);
        Assert.Equal(SolveStatus.Infeasible, result.Status);
        Assert.All(result.Weights, weight => Assert.True(double.IsNaN(weight)));
    }

    [Fact]
    public void A_setter_the_objective_does_not_take_is_refused_with_its_reason()
    {
        using var data = Samples.Data();
        var error = Assert.Throws<PortOptException>(() => Optimizer.Solve(new Spec { TargetVolatility = 0.2 }, data));
        Assert.Equal(Status.InvalidArgument, error.Status);
        Assert.Contains("target_volatility", error.Message, StringComparison.Ordinal);
    }

    [Fact]
    public void Missing_data_is_refused()
    {
        using var data = new MarketData(new MarketInputs { Covariance = Samples.Covariance });
        var error = Assert.Throws<PortOptException>(() => Optimizer.Solve(new Spec { Objective = Objective.MaxSharpe }, data));
        Assert.Equal(Status.InvalidArgument, error.Status);
    }

    [Fact]
    public void Solves_run_in_parallel_on_one_data()
    {
        using var data = Samples.Data();
        var spec = new Spec { Objective = Objective.MaxSharpe, Bounds = new Bounds(0, 1) };
        var expected = BitsOf(Optimizer.Solve(spec, data).Weights);
        var results = new long[64][];
        Parallel.For(0, results.Length, i => results[i] = BitsOf(Optimizer.Solve(spec, data).Weights));
        Assert.All(results, bits => Assert.Equal(expected, bits));
    }

    private static long[] BitsOf(double[] values)
    {
        return [.. values.Select(BitConverter.DoubleToInt64Bits)];
    }

    // The same problem through the C interface alone, as cpp/tests/wrapper_test.cpp does for the C++ wrapper.
    private static (double[] Weights, double[] Multipliers) SolveThroughC()
    {
        Library.EnsureLoaded();
        using var data = NativeData();
        Calls.Check(NativeMethods.po_spec_create(0, out var spec));
        using (spec)
        {
            Calls.Check(NativeMethods.po_spec_set_bounds(spec, 0, 1));
            Calls.Check(NativeMethods.po_solve(spec, data, out var result));
            using (result)
            {
                return (NativeWeights(result), NativeMultipliers(result));
            }
        }
    }

    private static unsafe DataHandle NativeData()
    {
        fixed (double* expected = Samples.Mu)
        fixed (double* matrix = Samples.Covariance)
        {
            Calls.Check(NativeMethods.po_data_create(3, expected, matrix, out var data));
            return data;
        }
    }

    private static unsafe double[] NativeWeights(ResultHandle result)
    {
        var weights = new double[3];
        fixed (double* output = weights)
        {
            Calls.Check(NativeMethods.po_result_weights(result, output, 3));
        }
        return weights;
    }

    private static double[] NativeMultipliers(ResultHandle result)
    {
        var multipliers = new double[4];
        for (var constraint = 0; constraint < 4; constraint++)
        {
            Calls.Check(NativeMethods.po_result_multiplier(result, constraint, out multipliers[constraint]));
        }
        return multipliers;
    }
}
