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


Series = tuple[NDArray[np.float64], NDArray[np.float64]]


def checked_returns(returns: ArrayLike | None, benchmark: ArrayLike | None) -> Series | None:
    if returns is None and benchmark is None:
        return None
    if returns is None or benchmark is None:
        raise ValueError("returns and benchmark must be given together")
    series, reference = as_doubles(returns), as_doubles(benchmark)
    if series.ndim != 2:
        raise ValueError(f"returns must be 2-D (periods by assets), got {series.ndim} dimensions")
    if reference.ndim != 1:
        raise ValueError(f"benchmark must be 1-D, got {reference.ndim} dimensions")
    if reference.size != series.shape[0]:
        raise ValueError(f"benchmark has {reference.size} entries, expected {series.shape[0]}")
    return series, reference


def asset_count(mu: NDArray[np.float64] | None, cov: NDArray[np.float64] | None, series: Series | None) -> int:
    sizes = {array.shape[0] for array in (mu, cov) if array is not None}
    if len(sizes) > 1:
        raise ValueError("mu and cov must have the same asset count")
    if series is not None:
        width = series[0].shape[1]
        if sizes and width not in sizes:
            raise ValueError(f"returns has {width} columns, expected {sizes.pop()}")
        sizes.add(width)
    if not sizes:
        raise ValueError("Data needs mu, cov or returns")
    return int(sizes.pop())


def optional_pointer(array: NDArray[np.float64] | None) -> Any:
    if array is None:
        return None
    return doubles_pointer(array)


class Data:
    def __init__(
        self,
        mu: ArrayLike | None,
        cov: ArrayLike | None,
        *,
        risk_free: float | None = None,
        returns: ArrayLike | None = None,
        benchmark: ArrayLike | None = None,
    ) -> None:
        mu_array, cov_array = checked_mu(mu), checked_cov(cov)
        series = checked_returns(returns, benchmark)
        count = asset_count(mu_array, cov_array, series)
        self._handle = _native.create_handle(
            lib.po_data_create, lib.po_data_free, count, optional_pointer(mu_array), optional_pointer(cov_array)
        )
        if risk_free is not None:
            _native.check(lib.po_data_set_risk_free(self.pointer, float(risk_free)))
        if series is not None:
            self._set_returns(*series)
        self.asset_count: int = _native.read_scalar(ctypes.c_size_t, lib.po_data_asset_count, self.pointer)
        self.weight_count: int = _native.read_scalar(ctypes.c_size_t, lib.po_data_weight_count, self.pointer)

    @property
    def pointer(self) -> int:
        return self._handle.pointer

    def _set_returns(self, returns: NDArray[np.float64], benchmark: NDArray[np.float64]) -> None:
        periods = returns.shape[0]
        _native.check(
            lib.po_data_set_returns(self.pointer, doubles_pointer(returns), doubles_pointer(benchmark), periods)
        )

    def portfolio_return(self, weights: ArrayLike) -> float:
        return self._metric(lib.po_portfolio_return, weights)

    def volatility(self, weights: ArrayLike) -> float:
        return self._metric(lib.po_portfolio_volatility, weights)

    def utility(self, weights: ArrayLike, risk_aversion: float) -> float:
        return self._metric(lib.po_utility, weights, float(risk_aversion))

    def tracking_error(self, weights: ArrayLike) -> float:
        return self._metric(lib.po_tracking_error, weights)

    def sharpe_ratio(self, weights: ArrayLike) -> float:
        return self._metric(lib.po_sharpe_ratio, weights)

    def diversification_ratio(self, weights: ArrayLike) -> float:
        return self._metric(lib.po_diversification_ratio, weights)

    def risk_contributions(self, weights: ArrayLike) -> NDArray[np.float64]:
        """c_i = w_i (Sigma w)_i / w' Sigma w, with 0 for the cash entry; all NaN when the volatility is 0."""
        w = as_doubles(weights)
        if w.ndim != 1:
            raise ValueError("weights must be 1-D")
        out = np.empty(w.size, dtype=np.float64)
        _native.check(
            lib.po_risk_contributions(self.pointer, doubles_pointer(w), w.size, doubles_pointer(out))
        )
        return out

    def _metric(self, function: Callable[..., int], weights: ArrayLike, *extra: float) -> float:
        w = as_doubles(weights)
        if w.ndim != 1:
            raise ValueError("weights must be 1-D")
        return float(_native.read_scalar(ctypes.c_double, function, self.pointer, doubles_pointer(w), w.size, *extra))
