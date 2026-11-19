namespace PortOpt;

/// <summary>
/// A call into the native library failed. <see cref="Status"/> is the status it returned and the message is the
/// library's own explanation, read on the thread that made the call.
/// </summary>
public sealed class PortOptException : Exception
{
    /// <summary>An exception with status <see cref="Status.InternalError"/> and no message.</summary>
    public PortOptException()
    {
    }

    /// <summary>An exception with status <see cref="Status.InternalError"/>.</summary>
    /// <param name="message">What failed.</param>
    public PortOptException(string message)
        : base(message)
    {
    }

    /// <summary>An exception with status <see cref="Status.InternalError"/> and the exception that caused it.</summary>
    /// <param name="message">What failed.</param>
    /// <param name="innerException">The cause.</param>
    public PortOptException(string message, Exception innerException)
        : base(message, innerException)
    {
    }

    /// <summary>An exception for a native status.</summary>
    /// <param name="status">The status the native call returned.</param>
    /// <param name="message">The native message.</param>
    public PortOptException(Status status, string message)
        : base(message)
    {
        Status = status;
    }

    /// <summary>The status the native call returned.</summary>
    public Status Status { get; } = Status.InternalError;
}
