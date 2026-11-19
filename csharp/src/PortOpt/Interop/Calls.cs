using System.Runtime.InteropServices;

namespace PortOpt.Interop;

internal static class Calls
{
    // The message of po_last_error belongs to the calling thread and lives only until its next PortOpt call, so it is
    // read here, before anything else can call into the library.
    public static void Check(int status)
    {
        if (status == 0)
        {
            return;
        }
        throw new PortOptException((Status)status, LastError());
    }

    public static string LastError()
    {
        return Marshal.PtrToStringUTF8(NativeMethods.po_last_error()) ?? "";
    }
}
