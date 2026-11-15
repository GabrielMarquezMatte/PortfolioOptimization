from dataclasses import FrozenInstanceError, replace

import numpy as np
import pytest

import portopt as po

COV = np.array([[0.04, 0.006, 0.002], [0.006, 0.09, 0.012], [0.002, 0.012, 0.16]])
MU = np.array([0.05, 0.08, 0.11])
RATE = 0.02


def run(objective: po.Objective, data: po.Data, **kwargs: object) -> po.Result:
    return po.solve(po.Spec(objective, **kwargs), data)  # type: ignore[arg-type]


def random_series(periods: int = 60) -> np.ndarray:
    return np.random.default_rng(7).normal(0.001, 0.01, size=(periods, 3))


def test_spec_is_frozen_dataclass() -> None:
    spec = po.Spec(po.Objective.UTILITY, bounds=(0.0, 1.0), risk_aversion=3.0)
    with pytest.raises(FrozenInstanceError):
        spec.target_return = 0.1  # type: ignore[misc]
    assert spec == po.Spec(po.Objective.UTILITY, bounds=(0.0, 1.0), risk_aversion=3.0)
    assert spec != po.Spec(po.Objective.UTILITY, bounds=(0.0, 1.0), risk_aversion=4.0)
    text = repr(spec)
    assert "risk_aversion=3.0" in text
    assert "objective=<Objective.UTILITY: 3>" in text
    assert "handle" not in text


def test_spec_with_array_bounds_compares_and_hashes_by_value() -> None:
    upper = np.array([0.5, 1.0, 1.0])
    first = po.Spec(po.Objective.MIN_VARIANCE, bounds=(np.zeros(3), upper))
    same = po.Spec(po.Objective.MIN_VARIANCE, bounds=([0.0, 0.0, 0.0], (0.5, 1.0, 1.0)))
    mixed = po.Spec(po.Objective.MIN_VARIANCE, bounds=(0.0, upper))
    other = po.Spec(po.Objective.MIN_VARIANCE, bounds=(np.zeros(3), np.array([0.5, 1.0, 0.9])))
    assert first == same
    assert hash(first) == hash(same)
    assert first != other
    targeted = po.Spec(po.Objective.MIN_VARIANCE, bounds=(0.0, upper), target_return=0.07)
    assert replace(mixed, target_return=0.07) == targeted
    assert po.solve(first, po.Data(MU, COV)).status == po.SolveStatus.OPTIMAL


def test_spec_broadcasts_a_scalar_side_against_an_array() -> None:
    mixed = po.Spec(po.Objective.MIN_VARIANCE, bounds=(0.0, [1.0, 1.0]))
    arrays = po.Spec(po.Objective.MIN_VARIANCE, bounds=([0.0, 0.0], [1.0, 1.0]))
    assert mixed == arrays
    assert hash(mixed) == hash(arrays)
    assert mixed.bounds == ((0.0, 0.0), (1.0, 1.0))
    assert po.Spec(po.Objective.MIN_VARIANCE, bounds=(0.0, 1.0)).bounds == (0.0, 1.0)
    assert mixed != po.Spec(po.Objective.MIN_VARIANCE, bounds=(0.0, 1.0))


def test_spec_replace_builds_validated_spec() -> None:
    data = po.Data(MU, COV)
    spec = po.Spec(po.Objective.MIN_VARIANCE)
    targeted = replace(spec, target_return=0.08)
    assert targeted.target_return == 0.08
    assert spec.target_return is None
    result = po.solve(targeted, data)
    assert result.status == po.SolveStatus.OPTIMAL
    assert data.portfolio_return(result.weights) == pytest.approx(0.08, abs=1e-12)
    with pytest.raises(ValueError, match="lambda must satisfy 0 <= lambda < 1"):
        replace(po.Spec(po.Objective.IMPORTANCE, importance=0.5), importance=1.0)


@pytest.mark.parametrize(
    "objective, setter, value, message",
    [
        (po.Objective.UTILITY, "risk_aversion", 0.0, "delta must be positive and finite"),
        (po.Objective.UTILITY, "risk_aversion", float("nan"), "delta must be positive and finite"),
        (po.Objective.IMPORTANCE, "importance", -0.1, "lambda must satisfy 0 <= lambda < 1"),
        (po.Objective.MAX_RETURN, "max_volatility", 0.0, "sigma must be positive and finite"),
        (po.Objective.MIN_VARIANCE, "target_return", float("inf"), "r must be finite"),
    ],
)
def test_invalid_setter_value_raises_at_construction(
    objective: po.Objective, setter: str, value: float, message: str
) -> None:
    with pytest.raises(ValueError, match=message):
        po.Spec(objective, **{setter: value})  # type: ignore[arg-type]


def test_utility_matches_closed_form() -> None:
    delta = 3.0
    ones = np.ones(3)
    inverse_ones, inverse_mu = np.linalg.solve(COV, ones), np.linalg.solve(COV, MU)
    nu = (ones @ inverse_mu - delta) / (ones @ inverse_ones)
    expected = (inverse_mu - nu * inverse_ones) / delta
    result = run(po.Objective.UTILITY, po.Data(MU, COV), risk_aversion=delta)
    assert result.status == po.SolveStatus.OPTIMAL
    np.testing.assert_allclose(result.weights, expected, rtol=0, atol=1e-12)
    assert np.isnan(result.volatility_multiplier)


def test_utility_default_risk_aversion_is_five() -> None:
    data = po.Data(MU, COV)
    np.testing.assert_array_equal(
        run(po.Objective.UTILITY, data).weights, run(po.Objective.UTILITY, data, risk_aversion=5.0).weights
    )


def test_importance_zero_is_min_variance() -> None:
    data = po.Data(MU, COV)
    importance = run(po.Objective.IMPORTANCE, data, importance=0.0, bounds=(0.0, 1.0))
    minimum = run(po.Objective.MIN_VARIANCE, data, bounds=(0.0, 1.0))
    np.testing.assert_allclose(importance.weights, minimum.weights, rtol=0, atol=1e-12)


def test_importance_matches_utility() -> None:
    data = po.Data(MU, COV)
    weight = 0.25
    importance = run(po.Objective.IMPORTANCE, data, importance=weight, bounds=(0.0, 1.0))
    utility = run(po.Objective.UTILITY, data, risk_aversion=2 * (1 - weight) / weight, bounds=(0.0, 1.0))
    np.testing.assert_allclose(importance.weights, utility.weights, rtol=0, atol=1e-12)


def test_kelly_is_utility_of_rate() -> None:
    plain = po.Data(MU, COV)
    np.testing.assert_allclose(
        run(po.Objective.KELLY, plain).weights,
        run(po.Objective.UTILITY, plain, risk_aversion=1.0).weights,
        rtol=0,
        atol=1e-12,
    )
    cash = po.Data(MU, COV, risk_free=RATE)
    kelly = run(po.Objective.KELLY, cash, bounds=(0.0, 1.0))
    utility = run(po.Objective.UTILITY, cash, risk_aversion=1.0 / (1.0 + RATE), bounds=(0.0, 1.0))
    assert kelly.weights.size == 4
    np.testing.assert_allclose(kelly.weights, utility.weights, rtol=0, atol=1e-12)


def test_tracking_error_replication() -> None:
    returns = random_series()
    held = np.array([0.5, 0.3, 0.2])
    data = po.Data(None, None, returns=returns, benchmark=returns @ held)
    assert data.asset_count == 3
    result = run(po.Objective.TRACKING_ERROR, data)
    assert result.status == po.SolveStatus.OPTIMAL
    np.testing.assert_allclose(result.weights, held, rtol=0, atol=1e-9)
    assert data.tracking_error(result.weights) == pytest.approx(0.0, abs=1e-12)


def test_max_return_on_capital_market_line() -> None:
    cap = 0.12
    excess = MU - RATE
    slope = np.sqrt(excess @ np.linalg.solve(COV, excess))
    expected = cap * np.linalg.solve(COV, excess) / slope
    data = po.Data(MU, COV, risk_free=RATE)
    result = run(po.Objective.MAX_RETURN, data, max_volatility=cap)
    assert result.status == po.SolveStatus.OPTIMAL
    np.testing.assert_allclose(result.weights[:3], expected, rtol=0, atol=1e-10)
    assert result.weights.sum() == pytest.approx(1.0, abs=1e-12)
    assert data.volatility(result.weights) == pytest.approx(cap, abs=1e-12)
    assert data.portfolio_return(result.weights) == pytest.approx(RATE + cap * slope, abs=1e-12)
    assert result.volatility_multiplier == pytest.approx(-slope, abs=1e-9)


def test_max_return_below_minimum_is_infeasible() -> None:
    minimum = 1.0 / np.sqrt(np.ones(3) @ np.linalg.solve(COV, np.ones(3)))
    result = run(po.Objective.MAX_RETURN, po.Data(MU, COV), max_volatility=0.99 * minimum)
    assert result.status == po.SolveStatus.INFEASIBLE
    assert np.isnan(result.weights).all()
    assert np.isnan(result.volatility_multiplier)


def test_min_variance_has_no_volatility_multiplier() -> None:
    assert np.isnan(run(po.Objective.MIN_VARIANCE, po.Data(MU, COV)).volatility_multiplier)


def test_setter_outside_matrix_raises_value_error() -> None:
    data = po.Data(MU, COV)
    with pytest.raises(ValueError, match="risk_aversion does not apply to MIN_VARIANCE"):
        run(po.Objective.MIN_VARIANCE, data, risk_aversion=3.0)
    with pytest.raises(ValueError, match="MAX_RETURN requires max_volatility"):
        run(po.Objective.MAX_RETURN, data)
    with pytest.raises(ValueError, match="UTILITY requires mu"):
        run(po.Objective.UTILITY, po.Data(None, COV))


def test_returns_without_benchmark_raises() -> None:
    with pytest.raises(ValueError, match="returns and benchmark must be given together"):
        po.Data(None, None, returns=random_series())
    with pytest.raises(ValueError, match="returns and benchmark must be given together"):
        po.Data(MU, COV, benchmark=np.zeros(60))


def test_returns_wrong_width_raises() -> None:
    returns = random_series()
    with pytest.raises(ValueError, match="returns has 2 columns, expected 3"):
        po.Data(MU, COV, returns=returns[:, :2], benchmark=returns[:, 0])
    with pytest.raises(ValueError, match="benchmark has 59 entries, expected 60"):
        po.Data(None, None, returns=returns, benchmark=returns[1:, 0])
    with pytest.raises(ValueError, match="returns must be 2-D"):
        po.Data(MU, COV, returns=returns[:, 0], benchmark=returns[:, 0])


def test_metrics_utility_and_tracking_error() -> None:
    returns = random_series()
    benchmark = returns[:, 0]
    data = po.Data(MU, COV, risk_free=RATE, returns=returns, benchmark=benchmark)
    w = np.array([0.3, 0.2, 0.1, 0.4])
    delta = 4.0
    expected_return = w[:3] @ MU + w[3] * RATE
    assert data.utility(w, delta) == pytest.approx(expected_return - delta / 2 * w[:3] @ COV @ w[:3], abs=1e-14)
    residual = returns @ w[:3] + w[3] * RATE - benchmark
    assert data.tracking_error(w) == pytest.approx(np.sqrt(residual @ residual / residual.size), abs=1e-14)
    with pytest.raises(ValueError, match="delta"):
        data.utility(w, -1.0)
    with pytest.raises(ValueError, match="returns"):
        po.Data(MU, COV).tracking_error(w[:3])
