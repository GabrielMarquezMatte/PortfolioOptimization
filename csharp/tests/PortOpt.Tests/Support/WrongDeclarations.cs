using System.Runtime.InteropServices;

namespace PortOpt.Tests.Support;

// Declarations the contract check must reject. They are only reflected, never called.
internal static unsafe partial class WrongDeclarations
{
    // const double* is an input: only double* fits, never out double.
    [DefaultDllImportSearchPaths(DllImportSearchPath.SafeDirectories)]
    [LibraryImport("libportopt")]
    internal static partial int po_spec_set_risk_budget(SpecHandle spec, out double budget, nuint count);

    // int32_t* is an output: out int or int*, never ref int.
    [DefaultDllImportSearchPaths(DllImportSearchPath.SafeDirectories)]
    [LibraryImport("libportopt")]
    internal static partial int po_result_status(ResultHandle result, ref int status);

    // po_spec** is a created handle: out SpecHandle, never the bare handle.
    [DefaultDllImportSearchPaths(DllImportSearchPath.SafeDirectories)]
    [LibraryImport("libportopt")]
    internal static partial int po_spec_create(int objective, SpecHandle spec);

    // A free takes the raw pointer the SafeHandle releases.
    [DefaultDllImportSearchPaths(DllImportSearchPath.SafeDirectories)]
    [LibraryImport("libportopt")]
    internal static partial void po_data_free(DataHandle data);
}
