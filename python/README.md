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

The distribution is `pyportopt`; the module is `portopt`. Python 3.11 or newer, with `numpy`. Until 1.0 the API may
change between minor versions. Documentation, the C and C++ interfaces and the changelog:
https://github.com/GabrielMarquezMatte/PortfolioOptimization
