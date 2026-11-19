using Microsoft.Win32.SafeHandles;

namespace PortOpt.Interop;

// One SafeHandle per opaque type of the C interface; the parameterless constructor is the one the LibraryImport
// generator uses for out parameters.

internal sealed class SpecHandle : SafeHandleZeroOrMinusOneIsInvalid
{
    public SpecHandle()
        : base(true)
    {
    }

    protected override bool ReleaseHandle()
    {
        NativeMethods.po_spec_free(handle);
        return true;
    }
}

internal sealed class DataHandle : SafeHandleZeroOrMinusOneIsInvalid
{
    public DataHandle()
        : base(true)
    {
    }

    protected override bool ReleaseHandle()
    {
        NativeMethods.po_data_free(handle);
        return true;
    }
}

internal sealed class ResultHandle : SafeHandleZeroOrMinusOneIsInvalid
{
    public ResultHandle()
        : base(true)
    {
    }

    protected override bool ReleaseHandle()
    {
        NativeMethods.po_result_free(handle);
        return true;
    }
}
