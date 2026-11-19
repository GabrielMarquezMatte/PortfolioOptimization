using PortOpt.Tests.Support;

namespace PortOpt.Tests;

public sealed class ContractTests
{
    [Fact]
    public void Every_contract_function_has_a_matching_declaration()
    {
        Assert.Empty(ContractCheck.FunctionProblems(Repository.ReadContract(), typeof(NativeMethods)));
    }

    [Fact]
    public void Every_contract_define_matches_a_member()
    {
        Assert.Empty(ContractCheck.DefineProblems(Repository.ReadContract()));
    }

    [Fact]
    public void A_missing_declaration_is_reported()
    {
        var problems = ContractCheck.FunctionProblems(["function po_status po_not_declared(int32_t)"], typeof(NativeMethods));
        Assert.Contains(problems, problem => problem.Contains("po_not_declared", StringComparison.Ordinal));
    }

    [Fact]
    public void A_declaration_outside_the_contract_is_reported()
    {
        var problems = ContractCheck.FunctionProblems(["function int32_t po_abi_version(void)"], typeof(NativeMethods));
        Assert.Contains(problems, problem => problem.Contains("po_solve", StringComparison.Ordinal));
    }

    [Fact]
    public void A_wrong_parameter_type_is_reported()
    {
        var line = "function po_status po_data_asset_count(const po_data*, double*)";
        var problems = ContractCheck.FunctionProblems([line], typeof(NativeMethods));
        Assert.Contains(problems, problem => problem.StartsWith("po_data_asset_count parameter 1", StringComparison.Ordinal));
    }

    [Fact]
    public void An_unmapped_define_is_reported()
    {
        var problems = ContractCheck.DefineProblems(["define PO_NEW_THING 3"]);
        Assert.Contains(problems, problem => problem.Contains("PO_NEW_THING", StringComparison.Ordinal));
    }

    [Fact]
    public void A_define_with_another_value_is_reported()
    {
        var problems = ContractCheck.DefineProblems(["define PO_OBJECTIVE_KELLY 9"]);
        Assert.Contains(problems, problem => problem.StartsWith("PO_OBJECTIVE_KELLY = 9", StringComparison.Ordinal));
    }

    [Theory]
    [InlineData("function po_status po_spec_set_risk_budget(po_spec*, const double*, size_t)", "po_spec_set_risk_budget parameter 1")]
    [InlineData("function po_status po_result_status(const po_result*, int32_t*)", "po_result_status parameter 1")]
    [InlineData("function po_status po_spec_create(int32_t, po_spec**)", "po_spec_create parameter 1")]
    [InlineData("function void po_data_free(po_data*)", "po_data_free parameter 0")]
    public void A_wrong_declaration_is_reported(string line, string expected)
    {
        var problems = ContractCheck.FunctionProblems([line], typeof(WrongDeclarations));
        Assert.Contains(problems, problem => problem.StartsWith(expected, StringComparison.Ordinal));
    }
}
