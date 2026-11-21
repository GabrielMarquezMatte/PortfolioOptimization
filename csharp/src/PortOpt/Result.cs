using PortOpt.Interop;

namespace PortOpt;

/// <summary>
/// The result of a solve, copied out of the native result, so it holds no handle. The weights, the multipliers and
/// the residual are NaN when the status provides no point (<see cref="SolveStatus.Infeasible"/>,
/// <see cref="SolveStatus.Unbounded"/>, <see cref="SolveStatus.NotConverged"/>); the multiplier of a constraint the
/// problem does not have is NaN. portopt.h has the multiplier convention.
/// </summary>
public sealed class Result
{
    private Result()
    {
    }

    /// <summary>How the solve ended.</summary>
    public SolveStatus Status { get; private init; }

    /// <summary>One weight per weight of the data, cash last; they sum to one.</summary>
    public double[] Weights { get; private init; } = [];

    /// <summary>The iterations the solver took.</summary>
    public int Iterations { get; private init; }

    /// <summary>The multiplier of the budget 1'w = 1.</summary>
    public double BudgetMultiplier { get; private init; }

    /// <summary>The multiplier of the target return.</summary>
    public double TargetReturnMultiplier { get; private init; }

    /// <summary>The multiplier of the volatility cap or target.</summary>
    public double VolatilityMultiplier { get; private init; }

    /// <summary>The multiplier of the gross exposure limit.</summary>
    public double LeverageMultiplier { get; private init; }

    /// <summary>One bound multiplier per weight: lower-bound multiplier minus upper-bound multiplier.</summary>
    public double[] BoundMultipliers { get; private init; } = [];

    /// <summary>For <see cref="Objective.RiskParity"/>, sum_i (c_i - b_i)^2 at the returned weights (c the risk
    /// contributions, b the budget); 0 for every other objective.</summary>
    public double Residual { get; private init; }

    internal static unsafe Result Read(ResultHandle handle)
    {
        Calls.Check(NativeMethods.po_result_status(handle, out var status));
        Calls.Check(NativeMethods.po_result_weight_count(handle, out var count));
        Calls.Check(NativeMethods.po_result_iterations(handle, out var iterations));
        Calls.Check(NativeMethods.po_result_residual(handle, out var residual));
        var weights = new double[checked((int)count)];
        var bounds = new double[weights.Length];
        fixed (double* output = weights)
        {
            Calls.Check(NativeMethods.po_result_weights(handle, output, count));
        }
        fixed (double* output = bounds)
        {
            Calls.Check(NativeMethods.po_result_bound_multipliers(handle, output, count));
        }
        return new Result
        {
            Status = (SolveStatus)status,
            Weights = weights,
            Iterations = iterations,
            BudgetMultiplier = Multiplier(handle, Constraint.Budget),
            TargetReturnMultiplier = Multiplier(handle, Constraint.TargetReturn),
            VolatilityMultiplier = Multiplier(handle, Constraint.Volatility),
            LeverageMultiplier = Multiplier(handle, Constraint.Leverage),
            BoundMultipliers = bounds,
            Residual = residual,
        };
    }

    private static double Multiplier(ResultHandle handle, Constraint constraint)
    {
        Calls.Check(NativeMethods.po_result_multiplier(handle, (int)constraint, out var multiplier));
        return multiplier;
    }
}
