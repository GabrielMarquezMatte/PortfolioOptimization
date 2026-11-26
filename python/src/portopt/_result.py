"""Solve result and the solve entry point."""

from __future__ import annotations

import ctypes
import enum
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from . import _native
from ._data import Data
from ._native import lib
from ._spec import Spec

PO_CONSTRAINT_BUDGET = 0
PO_CONSTRAINT_TARGET_RETURN = 1
PO_CONSTRAINT_VOLATILITY = 2
PO_CONSTRAINT_LEVERAGE = 3


class SolveStatus(enum.IntEnum):
    OPTIMAL = 0
    INFEASIBLE = 1
    UNBOUNDED = 2
    NOT_CONVERGED = 3


@dataclass(frozen=True, eq=False)
class Result:
    status: SolveStatus
    weights: NDArray[np.float64]
    iterations: int
    budget_multiplier: float
    target_return_multiplier: float
    volatility_multiplier: float
    leverage_multiplier: float
    bound_multipliers: NDArray[np.float64]
    residual: float


# SolveStatus(code) goes through EnumMeta.__call__; the table is one dict lookup. An unknown code still raises there.
STATUSES = {status.value: status for status in SolveStatus}


def solve_status(code: int) -> SolveStatus:
    status = STATUSES.get(code)
    if status is None:
        return SolveStatus(code)
    return status


def read_vectors(pointer: int, count: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """The weights and the bound multipliers, read into one buffer of 2 count doubles made for this call."""
    values = (ctypes.c_double * (2 * count))()
    _native.check(lib.po_result_weights(pointer, values, count))
    upper = (ctypes.c_double * count).from_buffer(values, count * ctypes.sizeof(ctypes.c_double))
    _native.check(lib.po_result_bound_multipliers(pointer, upper, count))
    both = np.frombuffer(values, dtype=np.float64)
    weights, bound_multipliers = both[:count], both[count:]
    weights.flags.writeable = False
    bound_multipliers.flags.writeable = False
    return weights, bound_multipliers


def read_doubles(pointer: int) -> tuple[float, float, float, float, float]:
    """The four named multipliers and the residual, through one output cell made for this call."""
    cell = ctypes.c_double()
    values = []
    for constraint in (PO_CONSTRAINT_BUDGET, PO_CONSTRAINT_TARGET_RETURN, PO_CONSTRAINT_VOLATILITY,
                       PO_CONSTRAINT_LEVERAGE):
        _native.check(lib.po_result_multiplier(pointer, constraint, cell))
        values.append(cell.value)
    _native.check(lib.po_result_residual(pointer, cell))
    return values[0], values[1], values[2], values[3], cell.value


def read_result(pointer: int, count: int) -> Result:
    """count is the weight count of the solved data, which po_result_weight_count would return (portopt.h)."""
    # Every output lives in this call: po_solve releases the GIL, so a shared buffer would race between threads.
    cell = ctypes.c_int32()
    _native.check(lib.po_result_status(pointer, cell))
    status = solve_status(cell.value)
    _native.check(lib.po_result_iterations(pointer, cell))
    weights, bound_multipliers = read_vectors(pointer, count)
    budget, target_return, volatility, leverage, residual = read_doubles(pointer)
    return Result(
        status=status,
        weights=weights,
        iterations=cell.value,
        budget_multiplier=budget,
        target_return_multiplier=target_return,
        volatility_multiplier=volatility,
        leverage_multiplier=leverage,
        bound_multipliers=bound_multipliers,
        residual=residual,
    )


def solve(spec: Spec, data: Data) -> Result:
    with _native.create_handle(lib.po_solve, lib.po_result_free, spec.pointer, data.pointer) as handle:
        return read_result(handle.pointer, data.weight_count)
