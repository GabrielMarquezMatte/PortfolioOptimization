"""Market data: expected returns, covariance and an optional risk-free asset."""

from __future__ import annotations

import ctypes
from collections.abc import Callable
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from . import _native
from ._native import lib
from ._spec import as_doubles, doubles_pointer


def checked_mu(mu: ArrayLike | None) -> NDArray[np.float64] | None:
    if mu is None:
        return None
    array = as_doubles(mu)
    if array.ndim != 1:
        raise ValueError(f"mu must be 1-D, got {array.ndim} dimensions")
    return array


def checked_cov(cov: ArrayLike | None) -> NDArray[np.float64] | None:
    if cov is None:
        return None
    array = as_doubles(cov)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError(f"cov must be a square 2-D array, got shape {array.shape}")
    return array


def asset_count(mu: NDArray[np.float64] | None, cov: NDArray[np.float64] | None) -> int:
    sizes = {array.shape[0] for array in (mu, cov) if array is not None}
    if not sizes:
        raise ValueError("Data needs mu or cov")
    if len(sizes) > 1:
        raise ValueError("mu and cov must have the same asset count")
    return int(sizes.pop())


def optional_pointer(array: NDArray[np.float64] | None) -> Any:
    if array is None:
        return None
    return doubles_pointer(array)


class Data:
    def __init__(self, mu: ArrayLike | None, cov: ArrayLike | None, *, risk_free: float | None = None) -> None:
        mu_array, cov_array = checked_mu(mu), checked_cov(cov)
        count = asset_count(mu_array, cov_array)
        self._handle = _native.create_handle(
            lib.po_data_create, lib.po_data_free, count, optional_pointer(mu_array), optional_pointer(cov_array)
        )
        if risk_free is not None:
            _native.check(lib.po_data_set_risk_free(self.pointer, float(risk_free)))
        self.asset_count: int = _native.read_scalar(ctypes.c_size_t, lib.po_data_asset_count, self.pointer)
        self.weight_count: int = _native.read_scalar(ctypes.c_size_t, lib.po_data_weight_count, self.pointer)

    @property
    def pointer(self) -> int:
        return self._handle.pointer

    def portfolio_return(self, weights: ArrayLike) -> float:
        return self._metric(lib.po_portfolio_return, weights)

    def volatility(self, weights: ArrayLike) -> float:
        return self._metric(lib.po_portfolio_volatility, weights)

    def _metric(self, function: Callable[..., int], weights: ArrayLike) -> float:
        w = as_doubles(weights)
        if w.ndim != 1:
            raise ValueError("weights must be 1-D")
        return float(_native.read_scalar(ctypes.c_double, function, self.pointer, doubles_pointer(w), w.size))
