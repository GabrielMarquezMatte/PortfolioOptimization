import numpy as np
import pytest

import portopt as po

COV = np.array([[0.04, 0.006, 0.002], [0.006, 0.09, 0.012], [0.002, 0.012, 0.16]])
MU = np.array([0.05, 0.08, 0.11])


def run(data: po.Data, **kwargs: object) -> po.Result:
    return po.solve(po.Spec(po.Objective.RISK_PARITY, **kwargs), data)  # type: ignore[arg-type]


def numpy_contributions(weights: np.ndarray, cov: np.ndarray) -> np.ndarray:
    marginal = cov @ weights
    contributions: np.ndarray = weights * marginal / (weights @ marginal)
    return contributions


def test_equal_risk_contribution_diagonal() -> None:
    variances = np.array([0.04, 0.09, 0.16, 0.25])
    result = run(po.Data(None, np.diag(variances)))
    assert result.status == po.SolveStatus.OPTIMAL
    expected = 1.0 / np.sqrt(variances)
    np.testing.assert_allclose(result.weights, expected / expected.sum(), rtol=1e-12)


def test_risk_budget_matches_contributions() -> None:
    budget = [0.5, 0.3, 0.2]
    data = po.Data(None, COV)
    result = run(data, risk_budget=budget)
    assert result.status == po.SolveStatus.OPTIMAL
    np.testing.assert_allclose(data.risk_contributions(result.weights), budget, rtol=1e-9)
    assert result.weights.sum() == pytest.approx(1.0, abs=1e-12)


def test_residual_zero_unconstrained() -> None:
    result = run(po.Data(None, COV))
    assert result.status == po.SolveStatus.OPTIMAL
    assert result.residual <= 1e-16
    assert result.budget_multiplier == 0.0
    assert np.isnan(result.target_return_multiplier)
    assert np.isnan(result.volatility_multiplier)
    assert np.isnan(result.leverage_multiplier)


def test_target_return_residual_positive_and_met() -> None:
    data = po.Data(MU, COV)
    free = run(data)
    target = data.portfolio_return(free.weights) + 0.015
    result = run(data, bounds=(0.0, 1.0), target_return=target)
    assert result.status == po.SolveStatus.OPTIMAL
    assert result.residual > 1e-8
    assert data.portfolio_return(result.weights) == pytest.approx(target, abs=1e-12)
    assert np.isfinite(result.target_return_multiplier)
    expected = np.sum((data.risk_contributions(result.weights) - 1.0 / 3.0) ** 2)
    assert result.residual == pytest.approx(expected, rel=1e-9)


def test_target_volatility_met() -> None:
    data = po.Data(None, COV)
    free = run(data)
    target = 1.2 * data.volatility(free.weights)
    result = run(data, target_volatility=target)
    assert result.status == po.SolveStatus.OPTIMAL
    assert data.volatility(result.weights) == pytest.approx(target, rel=1e-9)
    assert result.weights.sum() == pytest.approx(1.0, abs=1e-12)
    assert np.isfinite(result.volatility_multiplier)


def test_risk_contributions_sum_to_one_and_match_numpy() -> None:
    data = po.Data(MU, COV)
    weights = np.array([0.5, 0.3, 0.2])
    contributions = data.risk_contributions(weights)
    np.testing.assert_allclose(contributions, numpy_contributions(weights, COV), rtol=1e-14)
    assert contributions.sum() == pytest.approx(1.0, abs=1e-14)
    cash = po.Data(MU, COV, risk_free=0.02)
    with_cash = cash.risk_contributions([0.4, 0.3, 0.2, 0.1])
    assert with_cash[3] == 0.0
    assert with_cash.sum() == pytest.approx(1.0, abs=1e-14)
    assert np.isnan(data.risk_contributions(np.zeros(3))).all()
    with pytest.raises(ValueError, match="n must be 3"):
        data.risk_contributions([0.5, 0.5])


def test_cash_and_leverage_raise_value_error() -> None:
    with pytest.raises(ValueError, match="risk_free does not apply to RISK_PARITY"):
        run(po.Data(MU, COV, risk_free=0.02))
    with pytest.raises(ValueError, match="leverage does not apply to RISK_PARITY"):
        run(po.Data(MU, COV), leverage=2.0)
    with pytest.raises(ValueError, match="risk_budget has 2 entries; expected 3"):
        run(po.Data(MU, COV), risk_budget=[0.5, 0.5])
    with pytest.raises(ValueError, match="target_volatility does not apply to MIN_VARIANCE"):
        po.solve(po.Spec(po.Objective.MIN_VARIANCE, target_volatility=0.1), po.Data(MU, COV))
    with pytest.raises(ValueError, match="sigma must be positive and finite"):
        po.Spec(po.Objective.RISK_PARITY, target_volatility=0.0)
    with pytest.raises(ValueError, match="risk budget entries"):
        po.Spec(po.Objective.RISK_PARITY, risk_budget=[0.5, 0.6])


def test_unreachable_volatility_is_max_iter_with_nan() -> None:
    data = po.Data(None, np.diag([0.04, 0.09]))
    result = run(data, bounds=(0.0, 1.0), target_volatility=0.5)
    assert result.status == po.SolveStatus.NOT_CONVERGED
    assert np.isnan(result.weights).all()
    assert np.isnan(result.bound_multipliers).all()
    assert np.isnan(result.residual)
    assert np.isnan(result.budget_multiplier)
    assert np.isnan(result.volatility_multiplier)


def test_spec_with_risk_budget_array_hashes() -> None:
    first = po.Spec(po.Objective.RISK_PARITY, risk_budget=np.array([0.5, 0.3, 0.2]), target_volatility=0.2)
    same = po.Spec(po.Objective.RISK_PARITY, risk_budget=[0.5, 0.3, 0.2], target_volatility=0.2)
    other = po.Spec(po.Objective.RISK_PARITY, risk_budget=(0.5, 0.2, 0.3), target_volatility=0.2)
    assert first.risk_budget == (0.5, 0.3, 0.2)
    assert first == same
    assert hash(first) == hash(same)
    assert first != other
    assert len({first, same, other}) == 2
