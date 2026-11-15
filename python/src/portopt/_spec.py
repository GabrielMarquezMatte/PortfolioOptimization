"""Optimization specification."""

from __future__ import annotations

import ctypes
import enum
from dataclasses import KW_ONLY, dataclass, field

import numpy as np
from numpy.typing import ArrayLike, NDArray

from . import _native
from ._native import lib

_DOUBLES = ctypes.POINTER(ctypes.c_double)
SCALAR_SETTERS = (
    ("target_return", "po_spec_set_target_return"),
    ("risk_aversion", "po_spec_set_risk_aversion"),
    ("importance", "po_spec_set_importance"),
    ("max_volatility", "po_spec_set_max_volatility"),
)


class Objective(enum.IntEnum):
    MIN_VARIANCE = 0
    MAX_RETURN = 1
    MAX_SHARPE = 2
    UTILITY = 3
    IMPORTANCE = 4
    KELLY = 5
    TRACKING_ERROR = 6
    MAX_DIVERSIFICATION = 7
    RISK_PARITY = 8


def as_doubles(values: ArrayLike) -> NDArray[np.float64]:
    return np.ascontiguousarray(values, dtype=np.float64)


def doubles_pointer(array: NDArray[np.float64]) -> "ctypes._Pointer[ctypes.c_double]":
    return array.ctypes.data_as(_DOUBLES)


@dataclass(frozen=True)
class Spec:
    objective: Objective
    _: KW_ONLY
    bounds: tuple[float | ArrayLike, float | ArrayLike] | None = None
    target_return: float | None = None
    risk_aversion: float | None = None
    importance: float | None = None
    max_volatility: float | None = None
    _handle: _native.Handle = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        handle = _native.create_handle(lib.po_spec_create, lib.po_spec_free, int(self.objective))
        object.__setattr__(self, "_handle", handle)
        if self.bounds is not None:
            self._set_bounds(*self.bounds)
        for name, setter in SCALAR_SETTERS:
            value = getattr(self, name)
            if value is not None:
                _native.check(getattr(lib, setter)(self.pointer, float(value)))

    @property
    def pointer(self) -> int:
        return self._handle.pointer

    def _set_bounds(self, lower: float | ArrayLike, upper: float | ArrayLike) -> None:
        low, high = np.asarray(lower, dtype=np.float64), np.asarray(upper, dtype=np.float64)
        if low.ndim == 0 and high.ndim == 0:
            _native.check(lib.po_spec_set_bounds(self.pointer, float(low), float(high)))
            return
        low, high = np.broadcast_arrays(low, high)
        if low.ndim != 1:
            raise ValueError("bounds must be scalars or 1-D arrays")
        low, high = as_doubles(low), as_doubles(high)
        _native.check(
            lib.po_spec_set_asset_bounds(self.pointer, doubles_pointer(low), doubles_pointer(high), low.size)
        )
