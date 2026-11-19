using System.Runtime.InteropServices;

namespace PortOpt.Interop;

internal static unsafe partial class NativeMethods
{
    [LibraryImport(LibraryName)]
    internal static partial int po_solve(SpecHandle spec, DataHandle data, out ResultHandle result);

    [LibraryImport(LibraryName)] internal static partial int po_result_status(ResultHandle result, out int status);
    [LibraryImport(LibraryName)] internal static partial int po_result_weight_count(ResultHandle result, out nuint count);

    [LibraryImport(LibraryName)]
    internal static partial int po_result_weights(ResultHandle result, double* weights, nuint count);

    [LibraryImport(LibraryName)] internal static partial int po_result_iterations(ResultHandle result, out int iterations);

    [LibraryImport(LibraryName)]
    internal static partial int po_result_multiplier(ResultHandle result, int constraint, out double multiplier);

    [LibraryImport(LibraryName)]
    internal static partial int po_result_bound_multipliers(ResultHandle result, double* multipliers, nuint count);

    [LibraryImport(LibraryName)] internal static partial int po_result_residual(ResultHandle result, out double residual);
    [LibraryImport(LibraryName)] internal static partial void po_result_free(nint result);

    [LibraryImport(LibraryName)]
    internal static partial int po_portfolio_return(DataHandle data, double* weights, nuint count, out double value);

    [LibraryImport(LibraryName)]
    internal static partial int po_portfolio_volatility(DataHandle data, double* weights, nuint count, out double value);

    [LibraryImport(LibraryName)]
    internal static partial int po_utility(DataHandle data, double* weights, nuint count, double riskAversion,
                                           out double value);

    [LibraryImport(LibraryName)]
    internal static partial int po_tracking_error(DataHandle data, double* weights, nuint count, out double value);

    [LibraryImport(LibraryName)]
    internal static partial int po_sharpe_ratio(DataHandle data, double* weights, nuint count, out double value);

    [LibraryImport(LibraryName)]
    internal static partial int po_diversification_ratio(DataHandle data, double* weights, nuint count, out double value);

    [LibraryImport(LibraryName)]
    internal static partial int po_risk_contributions(DataHandle data, double* weights, nuint count, double* contributions);
}
