import numpy as np
import pytest

from portopt import (
    BlackLittermanPosterior,
    Views,
    bayes_stein,
    bl_equilibrium,
    bl_omega_he_litterman,
    bl_omega_idzorek,
    bl_posterior,
    sample_moments,
)

COV = np.array([[0.04, 0.012, 0.006], [0.012, 0.09, 0.018], [0.006, 0.018, 0.0225]])
W_MKT = np.array([0.5, 0.3, 0.2])


def test_sample_moments_match_numpy() -> None:
    returns = np.random.default_rng(1).normal(0.001, 0.02, size=(60, 4))
    mu, cov = sample_moments(returns)
    np.testing.assert_allclose(mu, returns.mean(axis=0), rtol=1e-14)
    np.testing.assert_allclose(cov, np.cov(returns, rowvar=False, ddof=1), rtol=1e-14)


def test_bayes_stein_matches_numpy() -> None:
    mu_hat, t, n = np.array([0.001, 0.002, 0.0005]), 50, 3
    cov = COV / 100
    result = bayes_stein(mu_hat, cov, t)
    sigma = cov * (t - 1) / (t - n - 2)
    ones = np.ones(n)
    mu0 = ones @ np.linalg.solve(sigma, mu_hat) / (ones @ np.linalg.solve(sigma, ones))
    d = mu_hat - mu0 * ones
    w = (n + 2) / ((n + 2) + t * d @ np.linalg.solve(sigma, d))
    np.testing.assert_allclose(result.grand_mean, mu0, rtol=1e-12)
    np.testing.assert_allclose(result.weight, w, rtol=1e-12)
    np.testing.assert_allclose(result.mu, (1 - w) * mu_hat + w * mu0 * ones, rtol=1e-12)


def test_relative_view_matches_numpy() -> None:
    tau, delta = 0.05, 2.5
    pi = bl_equilibrium(COV, W_MKT, delta)
    np.testing.assert_allclose(pi, delta * COV @ W_MKT, rtol=1e-14)
    p, q = np.array([[1.0, -1.0, 0.0]]), np.array([0.02])
    omega = bl_omega_he_litterman(COV, p, tau)
    np.testing.assert_allclose(omega, tau * p @ COV @ p.T, rtol=1e-14)
    posterior = bl_posterior(COV, pi, tau, Views(p, q, omega))
    a = p @ (tau * COV) @ p.T + omega
    b = tau * COV @ p.T
    np.testing.assert_allclose(posterior.mu, pi + b @ np.linalg.solve(a, q - p @ pi), rtol=1e-12)
    np.testing.assert_allclose(posterior.cov, COV + tau * COV - b @ np.linalg.solve(a, b.T), rtol=1e-12)


def test_full_confidence_reproduces_the_view() -> None:
    p, q, tau = np.array([[1.0, -1.0, 0.0], [0.0, 1.0, 1.0]]), np.array([0.02, 0.05]), 0.05
    pi = bl_equilibrium(COV, W_MKT, 2.5)
    omega = bl_omega_idzorek(COV, p, tau, [1.0, 1.0])
    mu = bl_posterior(COV, pi, tau, Views(p, q, omega)).mu
    np.testing.assert_allclose(p @ mu, q, atol=1e-15)


def test_no_views_returns_the_prior() -> None:
    pi = np.array([0.01, 0.02, 0.03])
    posterior = bl_posterior(COV, pi, 0.05, Views(np.empty((0, 3)), np.empty(0)))
    np.testing.assert_array_equal(posterior.mu, pi)
    np.testing.assert_allclose(posterior.cov, 1.05 * COV, rtol=1e-14)


def test_to_data_adds_the_rate() -> None:
    posterior = BlackLittermanPosterior(np.array([0.01, 0.02, 0.03]), COV)
    data = posterior.to_data(0.02)
    assert data.weight_count == 4
    assert data.portfolio_return([1.0, 0.0, 0.0, 0.0]) == pytest.approx(0.03)
    assert posterior.to_data().weight_count == 3


def test_views_reject_one_dimensional_p() -> None:
    with pytest.raises(ValueError, match="p must be 2-D"):
        Views(np.array([1.0, -1.0]), np.array([0.02]), np.array([[0.0]]))


def test_redundant_views_raise() -> None:
    p = np.array([[1.0, -1.0, 0.0], [2.0, -2.0, 0.0]])
    with pytest.raises(ValueError, match="views are redundant"):
        bl_omega_he_litterman(COV, p, 0.05)


def test_idzorek_rejects_a_short_confidence() -> None:
    p = np.array([[1.0, -1.0, 0.0], [0.0, 1.0, 1.0]])
    with pytest.raises(ValueError, match="confidence has 1 entries, expected 2"):
        bl_omega_idzorek(COV, p, 0.05, [1.0])


def test_p_width_must_match_cov() -> None:
    p = np.array([[1.0, -1.0]])
    with pytest.raises(ValueError, match="p has 2 columns, expected 3"):
        bl_omega_he_litterman(COV, p, 0.05)
    with pytest.raises(ValueError, match="p has 2 columns, expected 3"):
        bl_omega_idzorek(COV, p, 0.05, [0.5])
    with pytest.raises(ValueError, match="p has 2 columns, expected 3"):
        bl_posterior(COV, np.zeros(3), 0.05, Views(p, np.array([0.01]), np.array([[0.001]])))


def test_pi_size_must_match_cov() -> None:
    with pytest.raises(ValueError, match="same asset count"):
        bl_posterior(COV, np.zeros(2), 0.05, Views(np.empty((0, 3)), np.empty(0)))


def test_views_need_omega_when_not_empty() -> None:
    with pytest.raises(ValueError, match="omega is required"):
        Views(np.array([[1.0, -1.0, 0.0]]), np.array([0.02]))
