using System.Globalization;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text.RegularExpressions;

namespace PortOpt.Tests.Support;

// Compares the hand-written declarations of NativeMethods and the C# enums with the lines of
// include/portopt/portopt.abi. A declaration may use any C# type the table below accepts for its C type.
internal static partial class ContractCheck
{
    private static readonly Dictionary<string, Type> Scalars = new(StringComparer.Ordinal)
    {
        ["po_status"] = typeof(int),
        ["int32_t"] = typeof(int),
        ["size_t"] = typeof(nuint),
        ["double"] = typeof(double),
    };

    private static readonly string[] StatusNames =
        ["PO_OK", "PO_INVALID_ARGUMENT", "PO_INVALID_HANDLE", "PO_UNSUPPORTED", "PO_OUT_OF_MEMORY", "PO_INTERNAL_ERROR"];

    private static readonly Dictionary<string, string> Constants = new(StringComparer.Ordinal)
    {
        ["PO_ABI_VERSION"] = "AbiVersion",
    };

    private static readonly (string Prefix, Type Enum)[] Prefixes =
    [
        ("PO_OBJECTIVE_", typeof(Objective)),
        ("PO_SOLVE_", typeof(SolveStatus)),
        ("PO_CONSTRAINT_", typeof(Constraint)),
    ];

    [GeneratedRegex(@"^function (?<ret>.+?) (?<name>po_\w+)\((?<params>.*)\)$", RegexOptions.None, 1000)]
    private static partial Regex FunctionLine();

    [GeneratedRegex(@"^define (?<name>PO_\w+) (?<value>-?\d+)$", RegexOptions.None, 1000)]
    private static partial Regex DefineLine();

    [GeneratedRegex(@"\bconst\b", RegexOptions.None, 1000)]
    private static partial Regex ConstWord();

    public static List<string> FunctionProblems(IEnumerable<string> lines, Type native)
    {
        var declared = native.GetMethods(BindingFlags.Static | BindingFlags.Public | BindingFlags.NonPublic)
            .Where(method => method.GetCustomAttribute<LibraryImportAttribute>() is not null)
            .ToDictionary(method => method.Name, StringComparer.Ordinal);
        var problems = new List<string>();
        var seen = new HashSet<string>(StringComparer.Ordinal);
        foreach (var match in lines.Select(line => FunctionLine().Match(line)).Where(match => match.Success))
        {
            var name = match.Groups["name"].Value;
            seen.Add(name);
            if (!declared.TryGetValue(name, out var method))
            {
                problems.Add($"{name} is in the contract but not declared");
                continue;
            }
            problems.AddRange(SignatureProblems(method, match));
        }
        problems.AddRange(declared.Keys.Where(name => !seen.Contains(name))
            .Select(name => $"{name} is declared but not in the contract"));
        return problems;
    }

    public static List<string> DefineProblems(IEnumerable<string> lines)
    {
        var problems = new List<string>();
        var matched = new HashSet<(Type, string)>();
        foreach (var match in lines.Select(line => DefineLine().Match(line)).Where(match => match.Success))
        {
            var value = int.Parse(match.Groups["value"].Value, CultureInfo.InvariantCulture);
            var problem = DefineProblem(match.Groups["name"].Value, value, matched);
            if (problem is not null)
            {
                problems.Add(problem);
            }
        }
        foreach (var type in Prefixes.Select(entry => entry.Enum).Append(typeof(Status)))
        {
            problems.AddRange(Enum.GetNames(type).Where(member => !matched.Contains((type, member)))
                .Select(member => $"{type.Name}.{member} has no define"));
        }
        return problems;
    }

    private static string? DefineProblem(string name, int value, HashSet<(Type, string)> matched)
    {
        if (Constants.TryGetValue(name, out var field))
        {
            var constant = typeof(NativeMethods).GetField(field, BindingFlags.Static | BindingFlags.Public | BindingFlags.NonPublic);
            return Equals(constant?.GetRawConstantValue(), value) ? null : $"{name} = {value} but NativeMethods.{field} differs";
        }
        if (StatusNames.Contains(name, StringComparer.Ordinal))
        {
            return MemberProblem(name, name["PO_".Length..], value, typeof(Status), matched);
        }
        foreach (var (prefix, type) in Prefixes.Where(entry => name.StartsWith(entry.Prefix, StringComparison.Ordinal)))
        {
            return MemberProblem(name, name[prefix.Length..], value, type, matched);
        }
        return $"{name} has no C# destination";
    }

    private static string? MemberProblem(string name, string suffix, int value, Type type, HashSet<(Type, string)> matched)
    {
        var wanted = suffix.Replace("_", "", StringComparison.Ordinal);
        var member = Enum.GetNames(type)
            .FirstOrDefault(candidate => string.Equals(candidate, wanted, StringComparison.OrdinalIgnoreCase));
        if (member is null)
        {
            return $"{name} has no member in {type.Name}";
        }
        matched.Add((type, member));
        var actual = Convert.ToInt32(Enum.Parse(type, member), CultureInfo.InvariantCulture);
        return actual == value ? null : $"{name} = {value} but {type.Name}.{member} = {actual}";
    }

    private static IEnumerable<string> SignatureProblems(MethodInfo method, Match match)
    {
        var name = method.Name;
        var returned = match.Groups["ret"].Value;
        if (!Accepted(name, returned, true).Contains(method.ReturnType))
        {
            yield return $"{name} returns {method.ReturnType} for C {returned}";
        }
        var parameters = match.Groups["params"].Value;
        string[] cTypes = string.Equals(parameters, "void", StringComparison.Ordinal) ? [] : parameters.Split(", ");
        var actual = method.GetParameters();
        if (actual.Length != cTypes.Length)
        {
            yield return $"{name} has {actual.Length} parameters, the contract {cTypes.Length}";
            yield break;
        }
        for (var i = 0; i < cTypes.Length; i++)
        {
            var type = actual[i].ParameterType;
            if (!Accepted(name, cTypes[i], false).Contains(type) || (type.IsByRef && !actual[i].IsOut))
            {
                yield return $"{name} parameter {i} is {type}{(actual[i].IsOut ? " (out)" : "")} for C {cTypes[i]}";
            }
        }
    }

    // The C# types a declaration may use for one C type of the contract. A by-reference type is always an output and
    // must be declared out.
    private static Type[] Accepted(string function, string cType, bool isReturn)
    {
        var input = cType.StartsWith("const ", StringComparison.Ordinal);
        var bare = ConstWord().Replace(cType, "").Replace(" ", "", StringComparison.Ordinal);
        var stars = bare.Count(character => character == '*');
        var baseName = bare.TrimEnd('*');
        if (string.Equals(baseName, "void", StringComparison.Ordinal) && stars == 0)
        {
            return [typeof(void)];
        }
        if (string.Equals(baseName, "char", StringComparison.Ordinal))
        {
            return isReturn && stars == 1 ? [typeof(nint), typeof(byte).MakePointerType()] : [];
        }
        if (Scalars.TryGetValue(baseName, out var scalar))
        {
            return Pointers(scalar, stars, input);
        }
        var handle = typeof(NativeMethods).Assembly.GetType("PortOpt.Interop." + HandleName(baseName));
        if (handle is null)
        {
            return [];
        }
        if (stars == 1)
        {
            // po_data_set_risk_free also ends in _free: a free function is exactly the handle type's name plus _free.
            return string.Equals(function, baseName + "_free", StringComparison.Ordinal) ? [typeof(nint)] : [handle];
        }
        return stars == 2 ? [handle.MakeByRefType()] : [];
    }

    // const T* is an input buffer (T* only); T* is an output (out T) or an output buffer (T*).
    private static Type[] Pointers(Type scalar, int stars, bool input)
    {
        if (stars == 0)
        {
            return [scalar];
        }
        if (stars == 1)
        {
            return input ? [scalar.MakePointerType()] : [scalar.MakeByRefType(), scalar.MakePointerType()];
        }
        return [];
    }

    private static string HandleName(string cName)
    {
        var parts = cName["po_".Length..].Split('_').Select(part => char.ToUpperInvariant(part[0]) + part[1..]);
        return string.Concat(parts) + "Handle";
    }
}
