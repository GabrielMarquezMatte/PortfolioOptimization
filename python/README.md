# pyportopt

Portfolio optimization: minimum variance, maximum return under a volatility cap, maximum Sharpe ratio, mean-variance
utility, importance, Kelly, tracking error, maximum diversification and risk parity, with box bounds, a gross exposure
limit, cash and the multipliers of every constraint. The numerical core is the PortOpt C++23 library, which this
package bundles and calls through its C API.

```
pip install pyportopt
```

```python
import numpy as np
import portopt as po

mu = np.array([0.05, 0.08, 0.11])
cov = np.array([[0.04, 0.006, 0.002], [0.006, 0.09, 0.012], [0.002, 0.012, 0.16]])
data = po.Data(mu, cov)
result = po.solve(po.Spec(po.Objective.MIN_VARIANCE, bounds=(0.0, 1.0)), data)
print(result.status.name, result.weights)      # OPTIMAL [0.6234 0.2356 0.141 ] (rounded)
```

### From a return panel to a portfolio

The estimators are pure functions on arrays. Here: sample moments of daily returns, the Black-Litterman posterior with
one relative view (asset 0 beats asset 1 by 1%), and a maximum Sharpe solve.

```python
import numpy as np
import portopt as po

returns = np.random.default_rng(0).normal(0.0004, 0.01, size=(500, 3))  # daily returns, periods by assets
_, cov = po.sample_moments(returns)
cov = cov * 252                                                         # annualize after the estimation
pi = po.bl_equilibrium(cov, [0.5, 0.3, 0.2], delta=2.5)                 # implied returns, in excess of the rate
p, q = [[1.0, -1.0, 0.0]], [0.01]
views = po.Views(p, q, po.bl_omega_he_litterman(cov, p, tau=0.05))
posterior = po.bl_posterior(cov, pi, 0.05, views)                       # posterior.mu, posterior.cov
data = posterior.to_data(risk_free=0.02)                                # adds the rate back and sets the cash
result = po.solve(po.Spec(po.Objective.MAX_SHARPE, bounds=(0.0, 1.0)), data)
print(result.status.name, result.weights)                               # OPTIMAL [0.49 0.31 0.2 0.] (rounded)
```

`Pi`, `q` and the posterior mean are in excess of the risk-free rate, which is why `to_data` adds it back; without
`risk_free` the mean is left as it is and there is no cash. `bayes_stein(mu_hat, cov, t)` shrinks a mean, and
`bl_omega_idzorek` is the other `Omega`. Compute `bayes_stein` on the original (daily) data and annualize its `mu`
afterwards: its weight depends on the frequency. Each helper also takes a prebuilt `Data` in place of its arrays
(`bl_equilibrium(data, w_mkt, delta)`, `bl_omega_*(data, ...)`, `bayes_stein(data, t)`,
`bl_posterior(data, tau, views)`), so that a covariance used by several calls is validated and factored once; that
`Data` has no `risk_free`.

The distribution is `pyportopt`; the module is `portopt`. Python 3.11 or newer, with `numpy`. Until 1.0 the API may
change between minor versions. Documentation, the C and C++ interfaces and the changelog:
https://github.com/GabrielMarquezMatte/PortfolioOptimization
