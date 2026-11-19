namespace PortOpt.Tests.Support;

// The repository checkout the tests run from: the contract lives there, not in the build output.
internal static class Repository
{
    public static string Root { get; } = FindRoot();

    public static string[] ReadContract()
    {
        return File.ReadAllLines(Path.Combine(Root, "include", "portopt", "portopt.abi"));
    }

    private static string FindRoot()
    {
        for (var directory = new DirectoryInfo(AppContext.BaseDirectory); directory is not null; directory = directory.Parent)
        {
            if (File.Exists(Path.Combine(directory.FullName, "include", "portopt", "portopt.abi")))
            {
                return directory.FullName;
            }
        }
        throw new InvalidOperationException($"no include/portopt/portopt.abi above {AppContext.BaseDirectory}");
    }
}
