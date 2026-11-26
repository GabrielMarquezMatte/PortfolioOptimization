"""Estimation utilities: sample moments, Bayes-Stein shrinkage and the Black-Litterman model."""

from __future__ import annotations

import ctypes
import operator
from dataclasses import dataclass
from typing import Any, NamedTuple, overload

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

    # Identity, as Views: tuple equality compares the arrays elementwise and raises, and an array is unhashable.
    def __eq__(self, other: object) -> bool:
        return self is other

    def __ne__(self, other: object) -> bool:
        return self is not other

    def __hash__(self) -> int:
        return object.__hash__(self)


class BlackLittermanPosterior(NamedTuple):
    mu: Array
    cov: Array

    # Identity, as Views: tuple equality compares the arrays elementwise and raises, and an array is unhashable.
    def __eq__(self, other: object) -> bool:
        return self is other

    def __ne__(self, other: object) -> bool:
        return self is not other

    def __hash__(self) -> int:
        return object.__hash__(self)

    def to_data(self, risk_free: float | None = None) -> Data:
        """The posterior as market data; mu is in excess of the rate, so the rate is added back."""
        if risk_free is None:
            return Data(self.mu, self.cov)
        return Data(self.mu + risk_free, self.cov, risk_free=risk_free)


def checked_p(p: ArrayLike) -> Array:
    array = as_doubles(p)
    if array.ndim != 2:
        raise ValueError(f"p must be 2-D (views by assets), got {array.ndim} dimensions")
    return array


@dataclass(frozen=True, eq=False)
class Views:
    p: Array
    q: Array
    omega: Array | None = None

    def __post_init__(self) -> None:
        p, q = checked_p(self.p), as_doubles(self.q)
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


def required_cov(cov: ArrayLike) -> Array:
    array = checked_cov(cov)
    if array is None:
        raise ValueError("cov is required")
    return array


def checked_width(p: Array, n: int) -> Array:
    if p.shape[1] != n:
        raise ValueError(f"p has {p.shape[1]} columns, expected {n} (the asset count of cov)")
    return p


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


@overload
def bayes_stein(data: Data, /, t: int) -> BayesStein: ...


@overload
def bayes_stein(mu_hat: ArrayLike, cov: ArrayLike, t: int) -> BayesStein: ...


def bayes_stein(*args: Any, **kwargs: Any) -> BayesStein:
    """Jorion shrinkage; cov is the sample covariance with divisor t - 1.

    ``bayes_stein(mu_hat, cov, t)`` builds the market data of the call; ``bayes_stein(data, t)`` takes a prebuilt
    `Data` in place of both arrays (its mu is mu_hat, its covariance the sample covariance; no risk-free rate), so
    that the covariance is validated and factored once for every helper that uses it.
    """
    if args and isinstance(args[0], Data):
        return shrunk_mean(*args, **kwargs)
    return shrunk_mean_of_arrays(*args, **kwargs)


def checked_periods(t: int) -> int:
    t = operator.index(t)
    if t < 0:
        raise ValueError(f"t must be non-negative, got {t}")
    return t


def shrunk_mean_of_arrays(mu_hat: ArrayLike, cov: ArrayLike, t: int) -> BayesStein:
    periods = checked_periods(t)
    mu = required_vector(mu_hat, "mu_hat")
    return shrunk_mean(Data(mu, required_cov(cov)), periods)


def shrunk_mean(data: Data, /, t: int) -> BayesStein:
    periods = checked_periods(t)
    out = np.empty(data.asset_count, dtype=np.float64)
    mean, weight = ctypes.c_double(), ctypes.c_double()
    _native.check(lib.po_bayes_stein(data.pointer, periods, doubles_pointer(out), mean, weight))
    return BayesStein(out, mean.value, weight.value)


def covariance_data(cov: ArrayLike | Data) -> Data:
    """A prebuilt Data as it is, or the Data of the covariance alone."""
    if isinstance(cov, Data):
        return cov
    return Data(cov=required_cov(cov))


def bl_equilibrium(cov: ArrayLike | Data, w_mkt: ArrayLike, delta: float) -> Array:
    """The implied excess returns Pi = delta Sigma w_mkt.

    cov is Sigma, or a prebuilt `Data` whose covariance is Sigma (its mu, if any, is not used; no risk-free rate):
    the Data is then reused instead of being built, validated and factored on every call.
    """
    weights = required_vector(w_mkt, "w_mkt")
    data = covariance_data(cov)
    out = np.empty(weights.size, dtype=np.float64)
    _native.check(
        lib.po_bl_equilibrium(data.pointer, doubles_pointer(weights), weights.size, float(delta), doubles_pointer(out))
    )
    return out


def omega_inputs(cov: ArrayLike | Data, p: ArrayLike) -> tuple[Data, Array]:
    if isinstance(cov, Data):
        return cov, checked_width(checked_p(p), cov.asset_count)
    sigma = required_cov(cov)
    matrix = checked_width(checked_p(p), sigma.shape[0])
    return Data(cov=sigma), matrix


def bl_omega_he_litterman(cov: ArrayLike | Data, p: ArrayLike, tau: float) -> Array:
    """He and Litterman's diagonal Omega; cov is Sigma or a prebuilt `Data` of it, as in `bl_equilibrium`."""
    data, matrix = omega_inputs(cov, p)
    k = matrix.shape[0]
    with views_handle(matrix, np.zeros(k), None) as views:
        out = np.empty((k, k), dtype=np.float64)
        _native.check(lib.po_bl_omega_he_litterman(data.pointer, views.pointer, float(tau), doubles_pointer(out)))
    return out


def bl_omega_idzorek(cov: ArrayLike | Data, p: ArrayLike, tau: float, confidence: ArrayLike) -> Array:
    """Idzorek's diagonal Omega; cov is Sigma or a prebuilt `Data` of it, as in `bl_equilibrium`."""
    data, matrix = omega_inputs(cov, p)
    levels = required_vector(confidence, "confidence")
    k = matrix.shape[0]
    if levels.size != k:
        raise ValueError(f"confidence has {levels.size} entries, expected {k} (one per view)")
    with views_handle(matrix, np.zeros(k), None) as views:
        out = np.empty((k, k), dtype=np.float64)
        _native.check(
            lib.po_bl_omega_idzorek(
                data.pointer, views.pointer, float(tau), doubles_pointer(levels), doubles_pointer(out)
            )
        )
    return out


@overload
def bl_posterior(data: Data, /, tau: float, views: Views) -> BlackLittermanPosterior: ...


@overload
def bl_posterior(cov: ArrayLike, pi: ArrayLike, tau: float, views: Views) -> BlackLittermanPosterior: ...


def bl_posterior(*args: Any, **kwargs: Any) -> BlackLittermanPosterior:
    """Posterior mean (in excess of the rate) and covariance.

    ``bl_posterior(cov, pi, tau, views)`` builds the market data of the call; ``bl_posterior(data, tau, views)``
    takes a prebuilt `Data` in place of both arrays (its mu is Pi, its covariance Sigma; no risk-free rate).
    """
    if args and isinstance(args[0], Data):
        return posterior(*args, **kwargs)
    return posterior_of_arrays(*args, **kwargs)


def posterior_of_arrays(cov: ArrayLike, pi: ArrayLike, tau: float, views: Views) -> BlackLittermanPosterior:
    prior = required_vector(pi, "pi")
    sigma = required_cov(cov)
    checked_width(views.p, sigma.shape[0])
    return posterior(Data(prior, sigma), tau, views)


def posterior(data: Data, /, tau: float, views: Views) -> BlackLittermanPosterior:
    n = checked_width(views.p, data.asset_count).shape[1]
    mu, out = np.empty(n, dtype=np.float64), np.empty((n, n), dtype=np.float64)
    with views_handle(views.p, views.q, views.omega) as handle:
        _native.check(
            lib.po_bl_posterior(data.pointer, handle.pointer, float(tau), doubles_pointer(mu), doubles_pointer(out))
        )
    return BlackLittermanPosterior(mu, out)
