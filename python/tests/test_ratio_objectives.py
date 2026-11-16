from dataclasses import replace

import numpy as np
import pytest
from numpy.typing import NDArray

import portopt as po

COV = np.array([[0.04, 0.006, 0.002], [0.006, 0.09, 0.012], [0.002, 0.012, 0.16]])
MU = np.array([0.05, 0.08, 0.11])
RATE = 0.02


def run(objective: po.Objective, data: po.Data, **kwargs: object) -> po.Result:
    return po.solve(po.Spec(objective, **kwargs), data)  # type: ignore[arg-type]


def normalized(direction: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.asarray(direction / direction.sum(), dtype=np.float64)


def test_tangency_closed_form() -> None:
    expected = normalized(np.linalg.solve(COV, MU))
    result = run(po.Objective.MAX_SHARPE, po.Data(MU, COV))
    assert result.status == po.SolveStatus.OPTIMAL
    np.testing.assert_allclose(result.weights, expected, rtol=0, atol=1e-12)
    assert result.weights.sum() == pytest.approx(1.0, abs=1e-12)

    excess = MU - RATE
    with_cash = run(po.Objective.MAX_SHARPE, po.Data(MU, COV, risk_free=RATE))
    np.testing.assert_allclose(with_cash.weights[:3], normalized(np.linalg.solve(COV, excess)), rtol=0, atol=1e-12)


def test_sharpe_cash_weight_zero() -> None:
    data = po.Data(MU, COV, risk_free=RATE)
    result = run(po.Objective.MAX_SHARPE, data)
    assert result.status == po.SolveStatus.OPTIMAL
    assert result.weights.size == 4
    assert result.weights[3] == 0.0
    assert result.weights.sum() == pytest.approx(1.0, abs=1e-12)
    assert np.isnan(result.budget_multiplier)
    assert result.bound_multipliers.size == 4


def test_sharpe_unbounded_status() -> None:
    result = run(po.Objective.MAX_SHARPE, po.Data(MU, COV, risk_free=0.5))
    assert result.status == po.SolveStatus.UNBOUNDED
    assert np.isnan(result.weights).all()
    assert np.isnan(result.budget_multiplier)
    assert np.isnan(result.leverage_multiplier)
    assert np.isnan(result.bound_multipliers).all()


def test_sharpe_long_only_with_rate_above_every_mean_is_infeasible() -> None:
    result = run(po.Objective.MAX_SHARPE, po.Data(MU, COV, risk_free=0.5), bounds=(0.0, 1.0))
    assert result.status == po.SolveStatus.INFEASIBLE
    assert np.isnan(result.weights).all()


def test_diversification_closed_form() -> None:
    sigma = np.sqrt(np.diag(COV))
    expected = normalized(np.linalg.solve(COV, sigma))
    result = run(po.Objective.MAX_DIVERSIFICATION, po.Data(None, COV))
    assert result.status == po.SolveStatus.OPTIMAL
    np.testing.assert_allclose(result.weights, expected, rtol=0, atol=1e-12)
    assert po.Data(None, COV).diversification_ratio(result.weights) == pytest.approx(
        np.sqrt(sigma @ np.linalg.solve(COV, sigma)), abs=1e-12
    )


def test_leverage_caps_gross_exposure() -> None:
    cov = np.array([[0.04, 0.057, 0.0], [0.057, 0.09, 0.0], [0.0, 0.0, 0.16]])
    data = po.Data(np.array([0.10, 0.04, 0.08]), cov, risk_free=RATE)
    free = run(po.Objective.MAX_SHARPE, data)
    assert np.abs(free.weights).sum() > 4.0
    capped = run(po.Objective.MAX_SHARPE, data, leverage=3.0)
    assert capped.status == po.SolveStatus.OPTIMAL
    assert np.abs(capped.weights[:3]).sum() == pytest.approx(3.0, abs=1e-12)
    assert capped.weights.sum() == pytest.approx(1.0, abs=1e-12)
    assert capped.leverage_multiplier < 0.0
    assert data.sharpe_ratio(capped.weights) < data.sharpe_ratio(free.weights)
    loose = run(po.Objective.MAX_SHARPE, data, leverage=10.0)
    np.testing.assert_allclose(loose.weights, free.weights, rtol=0, atol=1e-12)
    assert loose.leverage_multiplier == 0.0


def test_leverage_below_one_without_cash_is_infeasible() -> None:
    result = run(po.Objective.MIN_VARIANCE, po.Data(MU, COV), leverage=0.9)
    assert result.status == po.SolveStatus.INFEASIBLE
    assert np.isnan(result.leverage_multiplier)


@pytest.mark.parametrize("value", [0.0, -1.0, float("nan"), float("inf")])
def test_invalid_leverage_raises_at_construction(value: float) -> None:
    with pytest.raises(ValueError, match="max_gross must be positive and finite"):
        po.Spec(po.Objective.MIN_VARIANCE, leverage=value)


def test_leverage_multiplier_nan_without_leverage() -> None:
    result = run(po.Objective.MIN_VARIANCE, po.Data(MU, COV))
    assert result.status == po.SolveStatus.OPTIMAL
    assert np.isnan(result.leverage_multiplier)


def test_spec_with_leverage_and_array_bounds_hashes() -> None:
    upper = np.array([0.5, 1.0, 1.0])
    first = po.Spec(po.Objective.MAX_SHARPE, bounds=(np.zeros(3), upper), leverage=1.5)
    same = po.Spec(po.Objective.MAX_SHARPE, bounds=([0.0, 0.0, 0.0], (0.5, 1.0, 1.0)), leverage=1.5)
    assert first == same
    assert hash(first) == hash(same)
    assert first != replace(first, leverage=2.0)
    assert first != po.Spec(po.Objective.MAX_SHARPE, bounds=(np.zeros(3), upper))
    assert "leverage=1.5" in repr(first)
    assert po.solve(first, po.Data(MU, COV)).status == po.SolveStatus.OPTIMAL


def test_ratio_metrics_match_numpy() -> None:
    w = np.array([0.5, 0.3, 0.1, 0.1])
    data = po.Data(MU, COV, risk_free=RATE)
    risky = w[:3]
    volatility = np.sqrt(risky @ COV @ risky)
    expected_return = risky @ MU + w[3] * RATE
    assert data.sharpe_ratio(w) == pytest.approx((expected_return - RATE) / volatility, abs=1e-14)
    assert data.diversification_ratio(w) == pytest.approx(np.sqrt(np.diag(COV)) @ risky / volatility, abs=1e-14)

    plain = po.Data(MU, COV)
    assert plain.sharpe_ratio(risky) == pytest.approx((MU @ risky) / volatility, abs=1e-14)
    with pytest.raises(ValueError, match="mu"):
        po.Data(None, COV).sharpe_ratio(risky)
    assert po.Data(None, COV).diversification_ratio(risky) == pytest.approx(
        np.sqrt(np.diag(COV)) @ risky / volatility, abs=1e-14
    )
