using System.Runtime.InteropServices;

namespace PortOpt.Interop;

internal static unsafe partial class NativeMethods
{
    [LibraryImport(LibraryName)]
    internal static partial int po_sample_moments(nuint t, nuint n, double* returns, double* outMu, double* outCov);

    [LibraryImport(LibraryName)]
    internal static partial int po_bayes_stein(DataHandle data, nuint t, double* outMu, out double outGrandMean, out double outWeight);

    [LibraryImport(LibraryName)]
    internal static partial int po_bl_equilibrium(DataHandle data, double* weights, nuint n, double delta, double* outPi);

    [LibraryImport(LibraryName)]
    internal static partial int po_views_create(nuint k, nuint n, double* p, double* q, out ViewsHandle views);

    [LibraryImport(LibraryName)] internal static partial int po_views_set_omega(ViewsHandle views, double* omega, nuint k);
    [LibraryImport(LibraryName)] internal static partial void po_views_free(nint views);

    [LibraryImport(LibraryName)]
    internal static partial int po_bl_omega_he_litterman(DataHandle data, ViewsHandle views, double tau, double* outOmega);

    [LibraryImport(LibraryName)]
    internal static partial int po_bl_omega_idzorek(
        DataHandle data, ViewsHandle views, double tau, double* confidence, double* outOmega);

    [LibraryImport(LibraryName)]
    internal static partial int po_bl_posterior(DataHandle data, ViewsHandle views, double tau, double* outMu, double* outCov);
}
