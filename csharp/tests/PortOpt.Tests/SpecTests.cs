namespace PortOpt.Tests;

public sealed class SpecTests
{
    [Fact]
    public void A_new_spec_is_minimum_variance_with_nothing_set()
    {
        var spec = new Spec();
        Assert.Equal(Objective.MinVariance, spec.Objective);
        Assert.Null(spec.Bounds);
        Assert.Null(spec.RiskBudget);
        spec.Validate();
    }

    [Fact]
    public void Equality_compares_the_arrays_by_value()
    {
        var left = new Spec { Objective = Objective.RiskParity, RiskBudget = [0.5, 0.25, 0.25], AssetBounds = new([0, 0, 0], [1, 1, 1]) };
        var right = new Spec { Objective = Objective.RiskParity, RiskBudget = [0.5, 0.25, 0.25], AssetBounds = new([0, 0, 0], [1, 1, 1]) };
        Assert.Equal(left, right);
        Assert.Equal(left.GetHashCode(), right.GetHashCode());
        Assert.NotEqual(left, right with { RiskBudget = [0.4, 0.3, 0.3] });
        Assert.NotEqual(left, right with { RiskBudget = null });
    }

    [Fact]
    public void The_risk_budget_is_copied()
    {
        double[] budget = [0.5, 0.25, 0.25];
        var spec = new Spec { Objective = Objective.RiskParity, RiskBudget = budget };
        budget[0] = 0.9;
        Assert.Equal(0.5, spec.RiskBudget![0]);
    }

    [Fact]
    public void A_value_outside_its_domain_fails_validation()
    {
        var error = Assert.Throws<PortOptException>(() => new Spec { Objective = Objective.Utility, RiskAversion = -1 }.Validate());
        Assert.Equal(Status.InvalidArgument, error.Status);
        Assert.NotEmpty(error.Message);
    }

    [Fact]
    public void Bounds_with_lower_above_upper_fail_validation()
    {
        var error = Assert.Throws<PortOptException>(() => new Spec { Bounds = new Bounds(1, 0) }.Validate());
        Assert.Equal(Status.InvalidArgument, error.Status);
    }

    [Fact]
    public void Both_forms_of_bounds_are_refused()
    {
        var spec = new Spec { Bounds = new Bounds(0, 1), AssetBounds = new([0, 0, 0], [1, 1, 1]) };
        Assert.Throws<InvalidOperationException>(spec.Validate);
    }

    [Fact]
    public void Asset_bounds_need_two_sides_of_one_length()
    {
        Assert.Throws<ArgumentException>(() => new AssetBounds([0, 0], [1, 1, 1]));
    }

    [Fact]
    public void An_unknown_objective_is_refused()
    {
        var error = Assert.Throws<PortOptException>(() => new Spec { Objective = (Objective)42 }.Validate());
        Assert.Equal(Status.InvalidArgument, error.Status);
    }
}
