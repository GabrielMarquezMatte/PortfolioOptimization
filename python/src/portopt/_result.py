"""Solve result and the solve entry point."""

from __future__ import annotations

import ctypes
import enum
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from . import _native
from ._data import Data
from ._native import lib
from ._spec import Spec, doubles_pointer

PO_CONSTRAINT_BUDGET = 0
PO_CONSTRAINT_TARGET_RETURN = 1
PO_CONSTRAINT_VOLATILITY = 2
PO_CONSTRAINT_LEVERAGE = 3


class SolveStatus(enum.IntEnum):
    OPTIMAL = 0
    INFEASIBLE = 1
    UNBOUNDED = 2
    MAX_ITER = 3


@dataclass(frozen=True)
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


def read_vector(function: Callable[..., int], pointer: int, count: int) -> NDArray[np.float64]:
    out = np.empty(count, dtype=np.float64)
    _native.check(function(pointer, doubles_pointer(out), count))
    out.flags.writeable = False
    return out


def read_multiplier(pointer: int, constraint: int) -> float:
    return float(_native.read_scalar(ctypes.c_double, lib.po_result_multiplier, pointer, constraint))


def read_result(pointer: int) -> Result:
    count = _native.read_scalar(ctypes.c_size_t, lib.po_result_weight_count, pointer)
    return Result(
        status=SolveStatus(_native.read_scalar(ctypes.c_int32, lib.po_result_status, pointer)),
        weights=read_vector(lib.po_result_weights, pointer, count),
        iterations=_native.read_scalar(ctypes.c_int32, lib.po_result_iterations, pointer),
        budget_multiplier=read_multiplier(pointer, PO_CONSTRAINT_BUDGET),
        target_return_multiplier=read_multiplier(pointer, PO_CONSTRAINT_TARGET_RETURN),
        volatility_multiplier=read_multiplier(pointer, PO_CONSTRAINT_VOLATILITY),
        leverage_multiplier=read_multiplier(pointer, PO_CONSTRAINT_LEVERAGE),
        bound_multipliers=read_vector(lib.po_result_bound_multipliers, pointer, count),
        residual=float(_native.read_scalar(ctypes.c_double, lib.po_result_residual, pointer)),
    )


def solve(spec: Spec, data: Data) -> Result:
    with _native.create_handle(lib.po_solve, lib.po_result_free, spec.pointer, data.pointer) as handle:
        return read_result(handle.pointer)
