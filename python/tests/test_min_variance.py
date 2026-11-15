import numpy as np
import pytest

import portopt as po

COV = np.array([[0.04, 0.006, 0.002], [0.006, 0.09, 0.012], [0.002, 0.012, 0.16]])
MU = np.array([0.05, 0.08, 0.11])


def min_variance(data: po.Data, **kwargs: object) -> po.Result:
    return po.solve(po.Spec(po.Objective.MIN_VARIANCE, **kwargs), data)  # type: ignore[arg-type]


def test_unconstrained_matches_closed_form() -> None:
    expected = np.linalg.solve(COV, np.ones(3))
    expected /= expected.sum()
    result = min_variance(po.Data(MU, COV))
    assert result.status == po.SolveStatus.OPTIMAL
    np.testing.assert_allclose(result.weights, expected, rtol=0, atol=1e-12)
    assert not result.weights.flags.writeable
    assert not result.bound_multipliers.flags.writeable


def test_long_only_weights_nonnegative_and_sum_to_one() -> None:
    result = min_variance(po.Data(MU, COV), bounds=(0.0, 1.0), target_return=0.1)
    assert result.status == po.SolveStatus.OPTIMAL
    assert (result.weights >= -1e-12).all()
    assert result.weights.sum() == pytest.approx(1.0, abs=1e-12)


def test_asset_bounds_vector_and_mixed_scalar() -> None:
    data = po.Data(MU, COV)
    upper = np.array([0.5, 1.0, 1.0])
    vectors = min_variance(data, bounds=(np.zeros(3), upper))
    mixed = min_variance(data, bounds=(0.0, upper))
    assert vectors.status == po.SolveStatus.OPTIMAL
    assert vectors.weights[0] == pytest.approx(0.5, abs=1e-12)
    np.testing.assert_array_equal(mixed.weights, vectors.weights)
    with pytest.raises(ValueError, match="entries"):
        min_variance(data, bounds=(np.zeros(2), np.ones(2)))


def test_frontier_with_list_comprehension() -> None:
    data = po.Data(MU, COV)
    minimum = data.volatility(min_variance(data).weights) ** 2
    right = [data.volatility(min_variance(data, target_return=t).weights) ** 2 for t in (0.07, 0.08, 0.09, 0.10, 0.11)]
    left = [data.volatility(min_variance(data, target_return=t).weights) ** 2 for t in (0.06, 0.05, 0.04, 0.03, 0.02)]
    assert minimum < min(right + left)
    assert right == sorted(right)
    assert left == sorted(left)


def test_cash_weight_is_last() -> None:
    data = po.Data(MU, COV, risk_free=0.02)
    result = min_variance(data, bounds=(0.0, 1.0))
    assert data.weight_count == 4
    assert result.weights.size == 4
    assert result.weights[-1] == pytest.approx(1.0, abs=1e-9)
    assert np.isnan(result.budget_multiplier)


def test_infeasible_returns_status_not_exception() -> None:
    result = min_variance(po.Data(MU, COV), bounds=(0.0, 0.2))
    assert result.status == po.SolveStatus.INFEASIBLE
    assert np.isnan(result.weights).all()


def test_invalid_covariance_raises_value_error() -> None:
    with pytest.raises(ValueError):
        po.Data(MU, np.ones((3, 2)))
    with pytest.raises(ValueError):
        po.Data(MU, np.array([[1.0, 2.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]))
    with pytest.raises(ValueError):
        po.Data(MU[:2], np.array([[1.0, 2.0], [2.0, 1.0]]))


def test_singular_covariance_raises_value_error() -> None:
    duplicated = np.array([[0.04, 0.04, 0.01], [0.04, 0.04, 0.01], [0.01, 0.01, 0.09]])
    with pytest.raises(ValueError, match="covariance is not positive definite"):
        po.Data(MU, duplicated)


def test_unsupported_objective_raises_not_implemented() -> None:
    with pytest.raises(NotImplementedError):
        po.solve(po.Spec(po.Objective.RISK_PARITY), po.Data(MU, COV))


def test_data_accepts_fortran_and_int_arrays() -> None:
    reference = min_variance(po.Data(MU, COV)).weights
    fortran = min_variance(po.Data(MU, np.asfortranarray(COV))).weights
    np.testing.assert_array_equal(fortran, reference)
    integer = min_variance(po.Data(np.array([1, 2, 3]), np.eye(3, dtype=int))).weights
    expected = min_variance(po.Data(np.array([1.0, 2.0, 3.0]), np.eye(3))).weights
    np.testing.assert_array_equal(integer, expected)


def test_metrics_match_numpy() -> None:
    data = po.Data(MU, COV, risk_free=0.02)
    w = np.array([0.3, 0.2, 0.1, 0.4])
    assert data.portfolio_return(w) == pytest.approx(w[:3] @ MU + w[3] * 0.02, abs=1e-14)
    assert data.volatility(w) == pytest.approx(np.sqrt(w[:3] @ COV @ w[:3]), abs=1e-14)
