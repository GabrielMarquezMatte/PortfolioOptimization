using PortOpt.Interop;

namespace PortOpt;

/// <summary>
/// What to solve, as values. Each property that is set goes to its <c>po_spec_set_*</c> setter, whose comment in
/// portopt.h gives its domain and the objectives that take it; one an objective does not take makes the solve fail
/// with <see cref="Status.InvalidArgument"/>. The setters run at <see cref="Optimizer.Solve"/>, or earlier with
/// <see cref="Validate"/>. Equality compares every value, the arrays included.
/// </summary>
public sealed record Spec
{
    private readonly double[]? _riskBudget;

    /// <summary>What to optimize; <see cref="Objective.MinVariance"/> by default.</summary>
    public Objective Objective { get; init; }

    /// <summary>One pair of bounds for every weight; at most one of this and <see cref="AssetBounds"/>. Without
    /// either every weight is free.</summary>
    public Bounds? Bounds { get; init; }

    /// <summary>One pair of bounds per weight, cash last; at most one of this and <see cref="Bounds"/>.</summary>
    public AssetBounds? AssetBounds { get; init; }

    /// <summary>The required portfolio return, finite.</summary>
    public double? TargetReturn { get; init; }

    /// <summary>The risk aversion delta of <see cref="Objective.Utility"/>, positive and finite.</summary>
    public double? RiskAversion { get; init; }

    /// <summary>The weight lambda of the return in <see cref="Objective.Importance"/>, 0 &lt;= lambda &lt; 1.</summary>
    public double? Importance { get; init; }

    /// <summary>The volatility cap of <see cref="Objective.MaxReturn"/>, positive and finite.</summary>
    public double? MaxVolatility { get; init; }

    /// <summary>The gross exposure limit sum_i |w_i| over the risky weights, positive and finite.</summary>
    public double? Leverage { get; init; }

    /// <summary>The volatility target of <see cref="Objective.RiskParity"/>, positive and finite.</summary>
    public double? TargetVolatility { get; init; }

    /// <summary>The risk budget of <see cref="Objective.RiskParity"/>: one positive entry per asset, summing to one;
    /// copied.</summary>
    public IReadOnlyList<double>? RiskBudget
    {
        get => _riskBudget;
        init => _riskBudget = value is null ? null : [.. value];
    }

    /// <summary>Runs every setter on a native spec, as the solve does, without solving.</summary>
    /// <exception cref="InvalidOperationException">Both <see cref="Bounds"/> and <see cref="AssetBounds"/> are set.</exception>
    /// <exception cref="PortOptException">A value is outside its domain.</exception>
    public void Validate()
    {
        using var handle = CreateHandle();
    }

    /// <inheritdoc/>
    public bool Equals(Spec? other)
    {
        return other is not null && Objective == other.Objective && Bounds == other.Bounds
               && Equals(AssetBounds, other.AssetBounds) && TargetReturn == other.TargetReturn
               && RiskAversion == other.RiskAversion && Importance == other.Importance
               && MaxVolatility == other.MaxVolatility && Leverage == other.Leverage
               && TargetVolatility == other.TargetVolatility && Values.Same(_riskBudget, other._riskBudget);
    }

    /// <inheritdoc/>
    public override int GetHashCode()
    {
        var scalars = HashCode.Combine(TargetReturn, RiskAversion, Importance, MaxVolatility, Leverage, TargetVolatility);
        return HashCode.Combine(Objective, Bounds, AssetBounds, scalars, Values.HashOf(_riskBudget));
    }

    internal SpecHandle CreateHandle()
    {
        if (Bounds is not null && AssetBounds is not null)
        {
            throw new InvalidOperationException("a spec takes Bounds or AssetBounds, not both");
        }
        Library.EnsureLoaded();
        Calls.Check(NativeMethods.po_spec_create((int)Objective, out var handle));
        try
        {
            SetBounds(handle);
            SetScalars(handle);
            SetRiskBudget(handle);
        }
        catch
        {
            handle.Dispose();
            throw;
        }
        return handle;
    }

    private unsafe void SetBounds(SpecHandle handle)
    {
        if (Bounds is { } bounds)
        {
            Calls.Check(NativeMethods.po_spec_set_bounds(handle, bounds.Lower, bounds.Upper));
        }
        if (AssetBounds is null)
        {
            return;
        }
        fixed (double* lower = AssetBounds.LowerSpan)
        fixed (double* upper = AssetBounds.UpperSpan)
        {
            Calls.Check(NativeMethods.po_spec_set_asset_bounds(handle, lower, upper, (nuint)AssetBounds.Lower.Count));
        }
    }

    private void SetScalars(SpecHandle handle)
    {
        Set(handle, TargetReturn, NativeMethods.po_spec_set_target_return);
        Set(handle, RiskAversion, NativeMethods.po_spec_set_risk_aversion);
        Set(handle, Importance, NativeMethods.po_spec_set_importance);
        Set(handle, MaxVolatility, NativeMethods.po_spec_set_max_volatility);
        Set(handle, Leverage, NativeMethods.po_spec_set_leverage);
        Set(handle, TargetVolatility, NativeMethods.po_spec_set_target_volatility);
    }

    private static void Set(SpecHandle handle, double? value, Func<SpecHandle, double, int> setter)
    {
        if (value is { } number)
        {
            Calls.Check(setter(handle, number));
        }
    }

    private unsafe void SetRiskBudget(SpecHandle handle)
    {
        if (_riskBudget is null)
        {
            return;
        }
        fixed (double* budget = _riskBudget)
        {
            Calls.Check(NativeMethods.po_spec_set_risk_budget(handle, budget, (nuint)_riskBudget.Length));
        }
    }
}
