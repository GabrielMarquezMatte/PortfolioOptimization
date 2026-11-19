using System.Runtime.InteropServices;
using PortOpt.Interop;

namespace PortOpt;

/// <summary>The native library this process loaded.</summary>
public static class NativeAbi
{
    /// <summary>The ABI version of the loaded library (<c>po_abi_version</c>); loading checks it is the one this
    /// package was built for.</summary>
    /// <exception cref="DllNotFoundException">The native library cannot be found.</exception>
    /// <exception cref="PortOptException">The library has another ABI version.</exception>
    public static int Version
    {
        get
        {
            Library.EnsureLoaded();
            return NativeMethods.po_abi_version();
        }
    }

    /// <summary>The version of the loaded library, "major.minor.patch" (<c>po_version</c>).</summary>
    /// <exception cref="DllNotFoundException">The native library cannot be found.</exception>
    /// <exception cref="PortOptException">The library has another ABI version.</exception>
    public static string LibraryVersion
    {
        get
        {
            Library.EnsureLoaded();
            return Marshal.PtrToStringUTF8(NativeMethods.po_version()) ?? "";
        }
    }
}
