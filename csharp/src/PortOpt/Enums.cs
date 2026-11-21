namespace PortOpt;

/// <summary>The status of a native call (<c>PO_OK</c>, <c>PO_INVALID_ARGUMENT</c>, ...).</summary>
public enum Status
{
    /// <summary>The call succeeded.</summary>
    Ok = 0,

    /// <summary>An argument is outside its domain, or the data does not have what the objective needs.</summary>
    InvalidArgument = 1,

    /// <summary>A handle is null.</summary>
    InvalidHandle = 2,

    /// <summary>The library does not support the request, or has another ABI version.</summary>
    Unsupported = 3,

    /// <summary>The library could not allocate.</summary>
    OutOfMemory = 4,

    /// <summary>An unexpected failure inside the library.</summary>
    InternalError = 5,
}

/// <summary>What a solve optimizes (<c>PO_OBJECTIVE_*</c>); portopt.h has the formula and the inputs of each.</summary>
public enum Objective
{
    /// <summary>Minimum variance w' Sigma w; needs the covariance, and the expected returns with a target
    /// return.</summary>
    MinVariance = 0,

    /// <summary>Maximum expected return under <see cref="Spec.MaxVolatility"/>; needs the expected returns and the
    /// covariance.</summary>
    MaxReturn = 1,

    /// <summary>Maximum Sharpe ratio; needs the expected returns and the covariance.</summary>
    MaxSharpe = 2,

    /// <summary>Minimum (delta / 2) w' Sigma w - mu~' w, delta = <see cref="Spec.RiskAversion"/> (5 by default).</summary>
    Utility = 3,

    /// <summary>Minimum (1 - lambda) w' Sigma w - lambda mu~' w, lambda = <see cref="Spec.Importance"/> (0.5 by
    /// default).</summary>
    Importance = 4,

    /// <summary>The Kelly portfolio: utility with delta = 1 / (1 + rate); needs 1 + rate &gt; 0.</summary>
    Kelly = 5,

    /// <summary>Minimum tracking error against the benchmark; needs return series of full column rank.</summary>
    TrackingError = 6,

    /// <summary>Maximum diversification ratio; needs the covariance.</summary>
    MaxDiversification = 7,

    /// <summary>Risk contributions as close as possible to <see cref="Spec.RiskBudget"/> (1 / n by default); no
    /// cash and no leverage.</summary>
    RiskParity = 8,
}

/// <summary>How a solve ended (<c>PO_SOLVE_*</c>).</summary>
public enum SolveStatus
{
    /// <summary>A solution: the weights and multipliers are set.</summary>
    Optimal = 0,

    /// <summary>Certified infeasible; the weights are NaN.</summary>
    Infeasible = 1,

    /// <summary>The supremum exists but is not attained; the weights are NaN.</summary>
    Unbounded = 2,

    /// <summary>The solver did not end at a point it can vouch for; the weights are NaN.</summary>
    NotConverged = 3,
}

// The named constraints of po_result_multiplier (PO_CONSTRAINT_*).
internal enum Constraint
{
    Budget = 0,
    TargetReturn = 1,
    Volatility = 2,
    Leverage = 3,
}
