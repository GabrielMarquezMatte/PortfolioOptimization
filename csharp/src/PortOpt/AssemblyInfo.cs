using System.Diagnostics.CodeAnalysis;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

// Every native signature is blittable (handles, integers, doubles, pointers), so the runtime never needs to marshal
// and the library stays trimming- and AOT-safe.
[assembly: DisableRuntimeMarshalling]

// The NuGet package puts the library where the host's own probing finds it, and PORTOPT_LIB loads it by full path. A
// Native AOT app has no deps.json: its publish puts the library next to the executable, which only AssemblyDirectory
// searches on Linux and macOS. The current directory of the process is never searched.
[assembly: DefaultDllImportSearchPaths(DllImportSearchPath.AssemblyDirectory | DllImportSearchPath.SafeDirectories)]
[assembly: SuppressMessage("Security", "CA5393:Do not use unsafe DllImportSearchPath value", Scope = "module",
    Justification = "The assembly directory is the application directory, which SafeDirectories already searches on Windows; on Linux and macOS it is the only place a Native AOT app finds the packaged library.")]
[assembly: InternalsVisibleTo("PortOpt.Tests")]
[assembly: InternalsVisibleTo("PortOpt.Benchmarks")]
