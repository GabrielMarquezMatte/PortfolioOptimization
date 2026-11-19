using BenchmarkDotNet.Attributes;
using PortOpt.Interop;

namespace PortOpt.Benchmarks;

// The cost of the binding: a long-only minimum-variance solve through Optimizer.Solve against the same calls made
// directly on NativeMethods (create the spec, solve, read status and weights, free). The difference is what the
// binding adds: the record-to-spec translation, the four multipliers and the bound multipliers it always reads, and
// the managed result.
[MemoryDiagnoser]
public class BindingOverhead
{
    private MarketData _data = null!;
    private Spec _spec = null!;
    private double[] _weights = [];

    [Params(10, 100)]
    public int Assets { get; set; }

    [GlobalSetup]
    public void Setup()
    {
        var mu = new double[Assets];
        var covariance = new double[Assets, Assets];
        for (var i = 0; i < Assets; i++)
        {
            mu[i] = 0.05 + (0.001 * i);
            for (var j = 0; j < Assets; j++)
            {
                covariance[i, j] = (i == j ? 0.04 : 0.0) + (0.01 / (1 + Math.Abs(i - j)));
            }
        }
        _data = new MarketData(mu, covariance);
        _spec = new Spec { Bounds = new Bounds(0, 1) };
        _weights = new double[Assets];
    }

    [GlobalCleanup]
    public void Cleanup()
    {
        _data.Dispose();
    }

    [Benchmark(Baseline = true)]
    public unsafe double Native()
    {
        Calls.Check(NativeMethods.po_spec_create(0, out var spec));
        using (spec)
        {
            Calls.Check(NativeMethods.po_spec_set_bounds(spec, 0, 1));
            Calls.Check(NativeMethods.po_solve(spec, _data.Handle, out var result));
            using (result)
            {
                Calls.Check(NativeMethods.po_result_status(result, out _));
                fixed (double* output = _weights)
                {
                    Calls.Check(NativeMethods.po_result_weights(result, output, (nuint)_weights.Length));
                }
                return _weights[0];
            }
        }
    }

    [Benchmark]
    public double Binding()
    {
        return Optimizer.Solve(_spec, _data).Weights[0];
    }
}
