using System.Runtime.InteropServices;

namespace PortOpt.Interop;

// The C interface, one declaration per function of include/portopt/portopt.abi, under its C name. ContractTests checks
// every declaration and constant against the contract.
internal static unsafe partial class NativeMethods
{
    // Without a lib prefix on Windows and with only the extension appended elsewhere, this one name finds
    // libportopt.dll, libportopt.so and libportopt.dylib.
    internal const string LibraryName = "libportopt";

    internal const int AbiVersion = 1;

    [LibraryImport(LibraryName)] internal static partial int po_abi_version();
    [LibraryImport(LibraryName)] internal static partial nint po_version();
    [LibraryImport(LibraryName)] internal static partial nint po_last_error();

    [LibraryImport(LibraryName)] internal static partial int po_spec_create(int objective, out SpecHandle spec);
    [LibraryImport(LibraryName)] internal static partial int po_spec_set_bounds(SpecHandle spec, double lower, double upper);

    [LibraryImport(LibraryName)]
    internal static partial int po_spec_set_weight_bounds(SpecHandle spec, double* lower, double* upper, nuint count);

    [LibraryImport(LibraryName)] internal static partial int po_spec_set_target_return(SpecHandle spec, double value);
    [LibraryImport(LibraryName)] internal static partial int po_spec_set_risk_aversion(SpecHandle spec, double value);
    [LibraryImport(LibraryName)] internal static partial int po_spec_set_importance(SpecHandle spec, double value);
    [LibraryImport(LibraryName)] internal static partial int po_spec_set_max_volatility(SpecHandle spec, double value);
    [LibraryImport(LibraryName)] internal static partial int po_spec_set_target_volatility(SpecHandle spec, double value);

    [LibraryImport(LibraryName)]
    internal static partial int po_spec_set_risk_budget(SpecHandle spec, double* budget, nuint count);

    [LibraryImport(LibraryName)] internal static partial int po_spec_set_leverage(SpecHandle spec, double value);
    [LibraryImport(LibraryName)] internal static partial void po_spec_free(nint spec);

    [LibraryImport(LibraryName)]
    internal static partial int po_data_create(nuint assets, double* mu, double* covariance, out DataHandle data);

    [LibraryImport(LibraryName)] internal static partial int po_data_set_risk_free(DataHandle data, double rate);

    [LibraryImport(LibraryName)]
    internal static partial int po_data_set_returns(DataHandle data, double* returns, double* benchmark, nuint periods);

    [LibraryImport(LibraryName)] internal static partial int po_data_asset_count(DataHandle data, out nuint count);
    [LibraryImport(LibraryName)] internal static partial int po_data_weight_count(DataHandle data, out nuint count);
    [LibraryImport(LibraryName)] internal static partial void po_data_free(nint data);
}
