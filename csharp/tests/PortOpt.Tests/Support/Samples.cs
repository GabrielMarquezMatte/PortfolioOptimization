namespace PortOpt.Tests.Support;

// The three-asset problem of cpp/tests/wrapper_test.cpp.
internal static class Samples
{
    public static double[] Mu => [0.08, 0.12, 0.10];

    public static double[,] Covariance => new[,]
    {
        { 0.040, 0.006, 0.010 },
        { 0.006, 0.090, 0.012 },
        { 0.010, 0.012, 0.060 },
    };

    public static MarketData Data()
    {
        return new MarketData(Mu, Covariance);
    }
}
