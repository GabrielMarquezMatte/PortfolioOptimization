using PortOpt.Interop;

namespace PortOpt;

/// <summary>The solve entry point (<c>po_solve</c>).</summary>
public static class Optimizer
{
    /// <summary>
    /// Solves the spec on the data. An infeasible, unbounded or unfinished problem is a result with that
    /// <see cref="Result.Status"/>, not an exception. Safe in parallel on the same spec and data.
    /// </summary>
    /// <param name="spec">What to solve.</param>
    /// <param name="data">The market data.</param>
    /// <returns>The result, which holds no native handle.</returns>
    /// <exception cref="InvalidOperationException">Both bounds forms are set on the spec.</exception>
    /// <exception cref="PortOptException">A value outside its domain, a setter the objective does not take, data
    /// the objective needs and does not have, or bounds or a risk budget of the wrong length.</exception>
    /// <exception cref="ObjectDisposedException">The data was disposed.</exception>
    public static Result Solve(Spec spec, MarketData data)
    {
        using var specHandle = spec.CreateHandle();
        Calls.Check(NativeMethods.po_solve(specHandle, data.Handle, out var resultHandle));
        using (resultHandle)
        {
            return Result.Read(resultHandle);
        }
    }
}
