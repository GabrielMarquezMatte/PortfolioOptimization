using System.Runtime.InteropServices;

namespace PortOpt.Interop;

internal static class Library
{
    private const string OverrideVariable = "PORTOPT_LIB";

    // A resolver can be registered once per assembly, so a failed load stays failed for the process.
    private static readonly Lazy<bool> _loaded = new(Load, LazyThreadSafetyMode.ExecutionAndPublication);

    public static void EnsureLoaded()
    {
        _ = _loaded.Value;
    }

    public static string? ResolvePath(string? environmentValue)
    {
        if (string.IsNullOrEmpty(environmentValue))
        {
            return null;
        }
        if (!File.Exists(environmentValue))
        {
            throw new DllNotFoundException($"{OverrideVariable} points to {environmentValue}, which is not a file");
        }
        return Path.GetFullPath(environmentValue);
    }

    private static bool Load()
    {
        var path = ResolvePath(Environment.GetEnvironmentVariable(OverrideVariable));
        if (path is not null)
        {
            NativeLibrary.SetDllImportResolver(typeof(NativeMethods).Assembly, (name, _, _) => Resolve(name, path));
        }
        var abi = NativeMethods.po_abi_version();
        if (abi != NativeMethods.AbiVersion)
        {
            var library = path ?? NativeMethods.LibraryName;
            throw new PortOptException(Status.Unsupported,
                $"{library} has ABI version {abi}, this package needs {NativeMethods.AbiVersion}");
        }
        return true;
    }

    private static nint Resolve(string name, string path)
    {
        if (!string.Equals(name, NativeMethods.LibraryName, StringComparison.Ordinal))
        {
            return 0;
        }
        return NativeLibrary.Load(path);
    }
}
