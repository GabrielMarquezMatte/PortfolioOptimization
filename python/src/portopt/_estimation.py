"""Estimation utilities: sample moments, Bayes-Stein shrinkage and the Black-Litterman model."""

from __future__ import annotations

import ctypes
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from . import _native
from ._data import Data, checked_cov, checked_mu
from ._native import lib
from ._spec import as_doubles, doubles_pointer

Array = NDArray[np.float64]


class BayesStein(NamedTuple):
    mu: Array
    grand_mean: float
    weight: float


class BlackLittermanPosterior(NamedTuple):
    mu: Array
    cov: Array

    def to_data(self, risk_free: float | None = None) -> Data:
        """The posterior as market data; mu is in excess of the rate, so the rate is added back."""
        if risk_free is None:
            return Data(self.mu, self.cov)
        return Data(self.mu + risk_free, self.cov, risk_free=risk_free)


@dataclass(frozen=True)
class Views:
    p: Array
    q: Array
    omega: Array | None = None

    def __post_init__(self) -> None:
        p, q = as_doubles(self.p), as_doubles(self.q)
        if p.ndim != 2:
            raise ValueError(f"p must be 2-D (views by assets), got {p.ndim} dimensions")
        if q.ndim != 1 or q.size != p.shape[0]:
            raise ValueError(f"q must be 1-D with {p.shape[0]} entries, got shape {q.shape}")
        object.__setattr__(self, "p", p)
        object.__setattr__(self, "q", q)
        object.__setattr__(self, "omega", self._checked_omega(p.shape[0]))

    def _checked_omega(self, k: int) -> Array | None:
        if self.omega is None:
            if k > 0:
                raise ValueError("omega is required when there are views")
            return None
        omega = as_doubles(self.omega)
        if omega.shape != (k, k):
            raise ValueError(f"omega must have shape ({k}, {k}), got {omega.shape}")
        return omega


def checked_p(p: ArrayLike) -> Array:
    array = as_doubles(p)
    if array.ndim != 2:
        raise ValueError(f"p must be 2-D (views by assets), got {array.ndim} dimensions")
    return array


def required_cov(cov: ArrayLike) -> Array:
    array = checked_cov(cov)
    if array is None:
        raise ValueError("cov is required")
    return array


def required_vector(values: ArrayLike, name: str) -> Array:
    array = checked_mu(values)
    if array is None:
        raise ValueError(f"{name} is required")
    return array


def views_handle(p: Array, q: Array, omega: Array | None) -> _native.Handle:
    k, n = p.shape
    handle = _native.create_handle(
        lib.po_views_create, lib.po_views_free, k, n, doubles_pointer(p), doubles_pointer(q)
    )
    if omega is not None and k > 0:
        _native.check(lib.po_views_set_omega(handle.pointer, doubles_pointer(omega), k))
    return handle


def sample_moments(returns: ArrayLike) -> tuple[Array, Array]:
    """Mean and covariance (divisor t - 1) of a periods-by-assets panel."""
    panel = as_doubles(returns)
    if panel.ndim != 2:
        raise ValueError(f"returns must be 2-D (periods by assets), got {panel.ndim} dimensions")
    t, n = panel.shape
    mu, cov = np.empty(n, dtype=np.float64), np.empty((n, n), dtype=np.float64)
    _native.check(lib.po_sample_moments(t, n, doubles_pointer(panel), doubles_pointer(mu), doubles_pointer(cov)))
    return mu, cov


def bayes_stein(mu_hat: ArrayLike, cov: ArrayLike, t: int) -> BayesStein:
    """Jorion shrinkage; cov is the sample covariance with divisor t - 1."""
    mu = required_vector(mu_hat, "mu_hat")
    data = Data(mu, required_cov(cov))
    out = np.empty(mu.size, dtype=np.float64)
    mean, weight = ctypes.c_double(), ctypes.c_double()
    _native.check(
        lib.po_bayes_stein(data.pointer, t, doubles_pointer(out), ctypes.byref(mean), ctypes.byref(weight))
    )
    return BayesStein(out, mean.value, weight.value)


def bl_equilibrium(cov: ArrayLike, w_mkt: ArrayLike, delta: float) -> Array:
    """The implied excess returns Pi = delta Sigma w_mkt."""
    weights = required_vector(w_mkt, "w_mkt")
    data = Data(cov=required_cov(cov))
    out = np.empty(weights.size, dtype=np.float64)
    _native.check(
        lib.po_bl_equilibrium(data.pointer, doubles_pointer(weights), weights.size, float(delta), doubles_pointer(out))
    )
    return out


def bl_omega_he_litterman(cov: ArrayLike, p: ArrayLike, tau: float) -> Array:
    matrix = checked_p(p)
    data = Data(cov=required_cov(cov))
    k = matrix.shape[0]
    with views_handle(matrix, np.zeros(k), None) as views:
        out = np.empty((k, k), dtype=np.float64)
        _native.check(lib.po_bl_omega_he_litterman(data.pointer, views.pointer, float(tau), doubles_pointer(out)))
    return out


def bl_omega_idzorek(cov: ArrayLike, p: ArrayLike, tau: float, confidence: ArrayLike) -> Array:
    matrix = checked_p(p)
    levels = required_vector(confidence, "confidence")
    data = Data(cov=required_cov(cov))
    k = matrix.shape[0]
    with views_handle(matrix, np.zeros(k), None) as views:
        out = np.empty((k, k), dtype=np.float64)
        _native.check(
            lib.po_bl_omega_idzorek(
                data.pointer, views.pointer, float(tau), doubles_pointer(levels), doubles_pointer(out)
            )
        )
    return out


def bl_posterior(cov: ArrayLike, pi: ArrayLike, tau: float, views: Views) -> BlackLittermanPosterior:
    """Posterior mean (in excess of the rate) and covariance."""
    prior = required_vector(pi, "pi")
    data = Data(prior, required_cov(cov))
    mu, out = np.empty(prior.size, dtype=np.float64), np.empty((prior.size,) * 2, dtype=np.float64)
    with views_handle(views.p, views.q, views.omega) as handle:
        _native.check(
            lib.po_bl_posterior(data.pointer, handle.pointer, float(tau), doubles_pointer(mu), doubles_pointer(out))
        )
    return BlackLittermanPosterior(mu, out)
