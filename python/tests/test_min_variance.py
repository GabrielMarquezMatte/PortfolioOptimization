import ctypes
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace

import numpy as np
import pytest
from numpy.typing import NDArray

import portopt as po
from portopt import _native, _result
from portopt._native import lib

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


def test_result_hashes_and_compares_by_identity() -> None:
    data = po.Data(MU, COV)
    first, second = min_variance(data), min_variance(data)
    assert first == first
    assert first != second
    assert not first == second
    assert {first, second} == {first, second}


def test_long_only_weights_nonnegative_and_sum_to_one() -> None:
    result = min_variance(po.Data(MU, COV), bounds=(0.0, 1.0), target_return=0.1)
    assert result.status == po.SolveStatus.OPTIMAL
    assert (result.weights >= -1e-12).all()
    assert result.weights.sum() == pytest.approx(1.0, abs=1e-12)


def test_weight_bounds_vector_and_mixed_scalar() -> None:
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
    spec = po.Spec(po.Objective.MIN_VARIANCE)

    def variance(target: float) -> float:
        return data.volatility(po.solve(replace(spec, target_return=target), data).weights) ** 2

    minimum = data.volatility(po.solve(spec, data).weights) ** 2
    right = [variance(target) for target in (0.07, 0.08, 0.09, 0.10, 0.11)]
    left = [variance(target) for target in (0.06, 0.05, 0.04, 0.03, 0.02)]
    assert minimum < min(right + left)
    assert right == sorted(right)
    assert left == sorted(left)


def test_cash_weight_is_last() -> None:
    data = po.Data(MU, COV, risk_free=0.02)
    result = min_variance(data, bounds=(0.0, 1.0))
    assert data.weight_count == 4
    with pytest.raises(AttributeError):
        data.weight_count = 3  # type: ignore[misc]
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


@pytest.mark.parametrize("objective", [9, 2**32, -(2**31) - 1])
def test_objective_outside_the_enum_raises_value_error(objective: int) -> None:
    # 2**32 must not wrap into the int32 of po_spec_create and solve MIN_VARIANCE.
    with pytest.raises(ValueError, match=r"objective is not one of PO_OBJECTIVE_\*"):
        po.Spec(objective)  # type: ignore[arg-type]


def test_unsupported_status_raises_not_implemented() -> None:
    # Every objective is solved now, so po_solve no longer returns PO_UNSUPPORTED; the mapping stays.
    with pytest.raises(NotImplementedError):
        _native.check(_native.PO_UNSUPPORTED)


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


def native_vector(function: Callable[..., int], pointer: int, count: int) -> NDArray[np.float64]:
    out = (ctypes.c_double * count)()
    _native.check(function(pointer, out, count))
    return np.array(out)


def test_result_is_what_the_native_readers_return() -> None:
    data = po.Data(MU, COV, risk_free=0.02)
    spec = po.Spec(po.Objective.MIN_VARIANCE, bounds=(0.0, 0.7), target_return=0.06)
    with _native.create_handle(lib.po_solve, lib.po_result_free, spec.pointer, data.pointer) as handle:
        pointer = handle.pointer
        result = _result.read_result(pointer, data.weight_count)
        assert _native.read_scalar(ctypes.c_size_t, lib.po_result_weight_count, pointer) == data.weight_count
        assert result.status is po.SolveStatus(_native.read_scalar(ctypes.c_int32, lib.po_result_status, pointer))
        assert result.iterations == _native.read_scalar(ctypes.c_int32, lib.po_result_iterations, pointer)
        assert np.array_equal(result.weights, native_vector(lib.po_result_weights, pointer, 4))
        assert np.array_equal(result.bound_multipliers, native_vector(lib.po_result_bound_multipliers, pointer, 4))
        multipliers = [_native.read_scalar(ctypes.c_double, lib.po_result_multiplier, pointer, c) for c in range(4)]
        named = [result.budget_multiplier, result.target_return_multiplier, result.volatility_multiplier,
                 result.leverage_multiplier]
        np.testing.assert_array_equal(named, multipliers)
        assert result.residual == _native.read_scalar(ctypes.c_double, lib.po_result_residual, pointer)
    assert type(result.status) is po.SolveStatus and type(result.iterations) is int
    for value in (*named, result.residual):
        assert type(value) is float
    for vector in (result.weights, result.bound_multipliers):
        assert vector.dtype == np.float64 and vector.shape == (4,) and not vector.flags.writeable
        with pytest.raises(ValueError):
            vector[0] = 1.0


def test_solves_in_threads_are_independent() -> None:
    data = po.Data(MU, COV, risk_free=0.02)
    specs = [po.Spec(po.Objective.MIN_VARIANCE, bounds=(0.0, 1.0), target_return=r) for r in (0.05, 0.07, 0.09)]
    expected = [po.solve(spec, data) for spec in specs]
    with ThreadPoolExecutor(max_workers=6) as pool:
        results = list(pool.map(lambda index: po.solve(specs[index % 3], data), range(300)))
    for index, result in enumerate(results):
        reference = expected[index % 3]
        assert np.array_equal(result.weights, reference.weights)
        assert np.array_equal(result.bound_multipliers, reference.bound_multipliers, equal_nan=True)
