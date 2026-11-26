import numpy as np
import pytest

from portopt import (
    BlackLittermanPosterior,
    Data,
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


def test_bayes_stein_rejects_a_negative_t() -> None:
    with pytest.raises(ValueError, match="t must be non-negative, got -1"):
        bayes_stein(np.array([0.001, 0.002, 0.0005]), COV / 100, -1)


def test_bayes_stein_rejects_a_float_t() -> None:
    with pytest.raises(TypeError):
        bayes_stein(np.array([0.001, 0.002, 0.0005]), COV / 100, 50.0)  # type: ignore[call-overload]


def test_views_hash_and_compare_by_identity() -> None:
    views = Views(np.empty((0, 3)), np.empty(0))
    assert hash(views) == hash(views)
    assert views == views


def test_bayes_stein_hashes_and_compares_by_identity() -> None:
    first = bayes_stein(np.array([0.001, 0.002, 0.0005]), COV / 100, 50)
    second = bayes_stein(np.array([0.001, 0.002, 0.0005]), COV / 100, 50)
    assert first == first
    assert first != second
    assert not first == second
    assert {first, second} == {first, second}
    mu, grand_mean, weight = first
    assert mu is first[0] and grand_mean == first.grand_mean and weight == first.weight


def test_posterior_hashes_and_compares_by_identity() -> None:
    first = BlackLittermanPosterior(np.array([0.01, 0.02, 0.03]), COV)
    second = BlackLittermanPosterior(np.array([0.01, 0.02, 0.03]), COV)
    assert first == first
    assert first != second
    assert not first == second
    assert {first, second} == {first, second}
    mu, cov = first
    assert mu is first[0] and cov is first.cov


def test_omegas_without_views_are_empty() -> None:
    assert bl_omega_he_litterman(COV, np.empty((0, 3)), 0.05).shape == (0, 0)
    assert bl_omega_idzorek(COV, np.empty((0, 3)), 0.05, np.empty(0)).shape == (0, 0)


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


def test_bayes_stein_of_a_data_is_the_bits_of_the_arrays() -> None:
    mu_hat, cov = np.array([0.001, 0.002, 0.0005]), COV / 100
    expected = bayes_stein(mu_hat, cov, 50)
    data = Data(mu_hat, cov)
    for result in (bayes_stein(data, 50), bayes_stein(data, t=50)):
        assert np.array_equal(result.mu, expected.mu)
        assert result.grand_mean == expected.grand_mean and result.weight == expected.weight


def test_bayes_stein_of_a_data_checks_t() -> None:
    data = Data(np.array([0.001, 0.002, 0.0005]), COV / 100)
    with pytest.raises(ValueError, match="t must be non-negative, got -1"):
        bayes_stein(data, -1)
    with pytest.raises(TypeError):
        bayes_stein(data, 50.0)  # type: ignore[call-overload]
    with pytest.raises(TypeError):
        bayes_stein(data, COV, 50)  # type: ignore[call-overload]


def test_bayes_stein_needs_the_mu_of_the_data() -> None:
    with pytest.raises(ValueError):
        bayes_stein(Data(cov=COV / 100), 50)


def test_one_data_gives_the_bits_of_the_covariance() -> None:
    p, q, tau = np.array([[1.0, -1.0, 0.0], [0.0, 1.0, 1.0]]), np.array([0.02, 0.05]), 0.05
    pi = bl_equilibrium(COV, W_MKT, 2.5)
    # The mu of the Data is not read by the equilibrium and the Omegas, so one Data serves the whole model.
    data = Data(pi, COV)
    assert np.array_equal(bl_equilibrium(data, W_MKT, 2.5), pi)
    assert np.array_equal(bl_equilibrium(Data(cov=COV), W_MKT, 2.5), pi)
    omega = bl_omega_he_litterman(COV, p, tau)
    assert np.array_equal(bl_omega_he_litterman(data, p, tau), omega)
    assert np.array_equal(bl_omega_idzorek(data, p, tau, [0.5, 1.0]), bl_omega_idzorek(COV, p, tau, [0.5, 1.0]))
    views = Views(p, q, omega)
    expected = bl_posterior(COV, pi, tau, views)
    for result in (bl_posterior(data, tau, views), bl_posterior(data, tau=tau, views=views)):
        assert np.array_equal(result.mu, expected.mu)
        assert np.array_equal(result.cov, expected.cov)


def test_a_data_is_checked_like_the_covariance() -> None:
    data = Data(np.zeros(3), COV)
    p = np.array([[1.0, -1.0]])
    with pytest.raises(ValueError, match="p has 2 columns, expected 3"):
        bl_omega_he_litterman(data, p, 0.05)
    with pytest.raises(ValueError, match="p has 2 columns, expected 3"):
        bl_omega_idzorek(data, p, 0.05, [0.5])
    with pytest.raises(ValueError, match="p has 2 columns, expected 3"):
        bl_posterior(data, 0.05, Views(p, np.array([0.01]), np.array([[0.001]])))
    with pytest.raises(ValueError):
        bl_equilibrium(data, [0.5, 0.5], 2.5)


def test_the_posterior_needs_the_mu_of_the_data() -> None:
    with pytest.raises(ValueError):
        bl_posterior(Data(cov=COV), 0.05, Views(np.empty((0, 3)), np.empty(0)))


def test_a_data_with_a_rate_is_refused() -> None:
    data = Data(np.zeros(3), COV, risk_free=0.02)
    with pytest.raises(ValueError):
        bl_equilibrium(data, W_MKT, 2.5)
    with pytest.raises(ValueError):
        bl_posterior(data, 0.05, Views(np.empty((0, 3)), np.empty(0)))
