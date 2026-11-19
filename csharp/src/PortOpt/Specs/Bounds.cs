using System.Runtime.InteropServices;

namespace PortOpt;

/// <summary>One pair of bounds for every weight, cash included (<c>po_spec_set_bounds</c>).</summary>
/// <param name="Lower">The lower bound; may be -inf.</param>
/// <param name="Upper">The upper bound; may be +inf.</param>
[StructLayout(LayoutKind.Auto)]
public readonly record struct Bounds(double Lower, double Upper);

/// <summary>
/// One pair of bounds per weight, cash last (<c>po_spec_set_asset_bounds</c>): as many pairs as the weight count of
/// the data. The values are copied; equality compares them.
/// </summary>
public sealed record AssetBounds
{
    private readonly double[] _lower;
    private readonly double[] _upper;

    /// <summary>Bounds per weight.</summary>
    /// <param name="lower">The lower bound of each weight.</param>
    /// <param name="upper">The upper bound of each weight, as many as <paramref name="lower"/>.</param>
    /// <exception cref="ArgumentException">The two sides have different lengths.</exception>
    public AssetBounds(IReadOnlyList<double> lower, IReadOnlyList<double> upper)
    {
        if (lower.Count != upper.Count)
        {
            throw new ArgumentException($"lower has {lower.Count} bounds and upper {upper.Count}", nameof(upper));
        }
        _lower = [.. lower];
        _upper = [.. upper];
    }

    /// <summary>The lower bound of each weight.</summary>
    public IReadOnlyList<double> Lower => _lower;

    /// <summary>The upper bound of each weight.</summary>
    public IReadOnlyList<double> Upper => _upper;

    /// <inheritdoc/>
    public bool Equals(AssetBounds? other)
    {
        return other is not null && Values.Same(_lower, other._lower) && Values.Same(_upper, other._upper);
    }

    /// <inheritdoc/>
    public override int GetHashCode()
    {
        return HashCode.Combine(Values.HashOf(_lower), Values.HashOf(_upper));
    }

    internal ReadOnlySpan<double> LowerSpan => _lower;

    internal ReadOnlySpan<double> UpperSpan => _upper;
}

internal static class Values
{
    public static bool Same(double[]? left, double[]? right)
    {
        if (left is null || right is null)
        {
            return left is null && right is null;
        }
        return left.SequenceEqual(right);
    }

    public static int HashOf(IReadOnlyList<double>? values)
    {
        var hash = new HashCode();
        foreach (var value in values ?? [])
        {
            hash.Add(value);
        }
        return hash.ToHashCode();
    }
}
