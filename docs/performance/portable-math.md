# Portable math

What replacing the platform `std::hypot` (the Givens of the QP, `src/qp/active_set.hpp`) and `std::log` (the Spinu
objective, `src/riskparity/spinu.cpp`) by `po::portable::hypot` and `po::portable::log` (`src/core/portable_math.cpp`)
costs. After the change `src/` calls no libm function besides `sqrt` and `fma`. Tools: [Measuring](measuring.md).

## Machine and runs

AMD Ryzen 7 5700X, Windows 10, on 2026-11-19, in a desktop session that is not isolated. A is commit `9f5c5c9`, B the
same plus the two replacements. Both binaries were kept side by side and run in the same session.

- Windows: `gcc-bench` (g++ 16.2.0, MSYS2 UCRT64, libm of the MinGW runtime).
- Linux: Docker on the same machine, `quay.io/pypa/manylinux_2_28_x86_64` (glibc 2.28), `gcc-toolset-14`
  (`-DCMAKE_BUILD_TYPE=Release -DPO_BUILD_BENCHMARKS=ON -DPO_BUILD_TESTS=OFF`, generic, no `-march=native`). The
  image has no `gcc-toolset-15`, and GCC 14 rejects `std::pair{"c", run_solve_c}` in
  `tests/benchmarks/capi/bm_capi.cpp` (a pair of a function type), so the container built a copy of the tree with
  `&run_solve_c` and `&run_solve_cpp`; no file of the repository changed for that.

Commands, per platform (after one discarded run of each binary): four rounds of `python tests/benchmarks/quick.py run
--exe EXE --out FILE` with A and B alternating, A first in every round. A single A then B pair, the procedure of
Measuring, was taken first and was not enough: on this machine two runs of the same binary moved some medians by 20 to
30% (A against A again), and one B run hit a disturbance that doubled a third of the cases, and on Linux the single
pair put B up to 30% above A on kernels that do not call the replaced functions. The table below is therefore the minimum over the four rounds of each median, plus the ratio of the
medians over the rounds and how many rounds B was not slower than A. A case counts as worse than noise only if B is
slower in both summaries and the gap of the minima exceeds the standard deviation of the rounds.

## Windows

| Case | A (us) | B (us) | B/A of the best rounds | B/A of the median rounds | rounds with B <= A |
|---|---|---|---|---|---|
| `min_variance/long_only/10` | 3.324 | 3.021 | 0.909 | 0.751 | 4/4 |
| `min_variance/long_only/100` | 227.5 | 235.1 | 1.033 | 0.727 | 2/4 |
| `max_return/long_only/10` | 30.92 | 28.84 | 0.933 | 0.952 | 4/4 |
| `max_return/long_only/100` | 2647 | 2517 | 0.951 | 0.960 | 2/4 |
| `max_sharpe/long_only/10` | 3.958 | 3.681 | 0.930 | 0.929 | 3/4 |
| `max_sharpe/long_only/100` | 277.3 | 250.1 | 0.902 | 1.092 | 1/4 |
| `utility/box/10` | 3.939 | 3.787 | 0.961 | 0.963 | 3/4 |
| `utility/box/100` | 456.4 | 444 | 0.973 | 1.036 | 1/4 |
| `tracking_error/long_only/10` | 4.465 | 4.65 | 1.042 | 1.076 | 0/4 |
| `tracking_error/long_only/100` | 980.8 | 906.4 | 0.924 | 0.951 | 3/4 |
| `max_diversification/box/10` | 4.806 | 4.339 | 0.903 | 0.962 | 3/4 |
| `max_diversification/box/100` | 383 | 380.8 | 0.994 | 0.935 | 3/4 |
| `risk_parity/budget/10` | 6.164 | 6.156 | 0.999 | 1.035 | 1/4 |
| `risk_parity/budget/100` | 639.7 | 621.8 | 0.972 | 1.048 | 1/4 |
| `risk_parity/target_volatility/10` | 231.6 | 236.4 | 1.021 | 0.951 | 2/4 |
| `risk_parity/target_volatility/100` | 1.555e+04 | 1.548e+04 | 0.996 | 1.040 | 1/4 |
| `risk_parity/box/10` | 139.5 | 135.4 | 0.971 | 0.961 | 3/4 |
| `risk_parity/box/100` | 3.138e+04 | 2.867e+04 | 0.914 | 0.833 | 4/4 |
| `kernel/cholesky/100` | 34.53 | 30.82 | 0.893 | 0.920 | 4/4 |
| `kernel/solve_lower/100` | 2.247 | 2.256 | 1.004 | 0.943 | 2/4 |
| `kernel/qr/100` | 335.2 | 348 | 1.038 | 1.058 | 1/4 |
| `kernel/qp_iteration/100` | 534.3 | 509.3 | 0.953 | 1.014 | 3/4 |
| `kernel/compensated_dot/100` | 0.3098 | 0.3142 | 1.014 | 1.034 | 2/4 |
| `kernel/covariance_times/100` | 30.59 | 31.44 | 1.028 | 1.034 | 2/4 |
| `kernel/gauss_newton/100` | 65.15 | 69.21 | 1.062 | 1.016 | 2/4 |
| `kernel/empty_active_set/100` | 44.93 | 45.39 | 1.010 | 1.049 | 2/4 |

The one case that tripped the rule, `tracking_error/long_only/10` (B 4% slower in the first four rounds, 0/4), was
repeated with six alternating rounds of 15 repetitions: B was not slower in five of six (medians 4.266 against 4.302
us at the best round); A alone ranges from 4.3 to 9.3 us across those rounds, so the 4% is inside the machine's noise.

## Linux

| Case | A (us) | B (us) | B/A of the best rounds | B/A of the median rounds | rounds with B <= A |
|---|---|---|---|---|---|
| `kernel/cholesky/100` | 29.16 | 30.5 | 1.046 | 1.064 | 2/4 |
| `kernel/solve_lower/100` | 2.266 | 2.706 | 1.194 | 1.174 | 1/4 |
| `kernel/qr/100` | 401.1 | 516 | 1.287 | 1.340 | 1/4 |
| `kernel/qp_iteration/100` | 497 | 480.1 | 0.966 | 1.048 | 2/4 |
| `kernel/compensated_dot/100` | 0.3229 | 0.3035 | 0.940 | 1.020 | 2/4 |
| `kernel/covariance_times/100` | 30.21 | 29.48 | 0.976 | 0.806 | 4/4 |
| `kernel/gauss_newton/100` | 67.44 | 66.7 | 0.989 | 1.004 | 2/4 |
| `kernel/empty_active_set/100` | 44.65 | 46.5 | 1.041 | 1.045 | 1/4 |
| `min_variance/long_only/10` | 2.697 | 2.498 | 0.926 | 0.974 | 3/4 |
| `min_variance/long_only/100` | 418.4 | 387.9 | 0.927 | 0.872 | 3/4 |
| `max_return/long_only/10` | 25.75 | 23.9 | 0.928 | 1.007 | 2/4 |
| `max_return/long_only/100` | 2527 | 2414 | 0.955 | 0.951 | 3/4 |
| `max_sharpe/long_only/10` | 3.281 | 3.042 | 0.927 | 0.955 | 3/4 |
| `max_sharpe/long_only/100` | 288.7 | 261.5 | 0.906 | 0.970 | 3/4 |
| `utility/box/10` | 3.31 | 3.098 | 0.936 | 0.881 | 4/4 |
| `utility/box/100` | 404.9 | 392.5 | 0.969 | 0.901 | 4/4 |
| `tracking_error/long_only/10` | 3.803 | 3.816 | 1.003 | 0.923 | 3/4 |
| `tracking_error/long_only/100` | 1046 | 1066 | 1.019 | 1.040 | 1/4 |
| `max_diversification/box/10` | 3.983 | 3.885 | 0.975 | 0.997 | 3/4 |
| `max_diversification/box/100` | 424.1 | 414.5 | 0.977 | 1.054 | 1/4 |
| `risk_parity/budget/10` | 4.651 | 4.383 | 0.942 | 1.058 | 1/4 |
| `risk_parity/budget/100` | 594.8 | 572.6 | 0.963 | 1.024 | 2/4 |
| `risk_parity/target_volatility/10` | 176.9 | 171.6 | 0.970 | 1.118 | 1/4 |
| `risk_parity/target_volatility/100` | 1.423e+04 | 1.421e+04 | 0.999 | 1.039 | 2/4 |
| `risk_parity/box/10` | 99.39 | 96.84 | 0.974 | 0.867 | 4/4 |
| `risk_parity/box/100` | 3.197e+04 | 2.887e+04 | 0.903 | 0.911 | 4/4 |

The two cases marked in the run, `kernel/solve_lower/100` and `kernel/qr/100`, call neither `hypot` nor `log` (nor
`make_givens`): `solve_lower` is the same machine code in A and B, so the difference is layout or machine state, not
the change. Every `BM_Solve` case is within noise or faster.

## Solutions

`portopt_profile --max-n 100 --json` before and after, `quick.py compare-solutions`:

- Windows: 70 of 70 cases bit for bit identical, same status and iterations. The MinGW `hypot` and `log` already
  rounded as the portable ones on every value of the suite.
- Linux: no status or iteration count changed; 35 of 70 cases differ in the last bits (worst relative weight difference
  3.8e-15, worst objective 1.1e-14, both far inside the KKT bands), because glibc's `hypot` is not the portable one.
  This is the intended effect: the Linux solution now follows the Windows one. The two platforms are still not bit for
  bit with each other on this suite, because the instances (`tests/unit/random_covariance.hpp`) are generated with
  the platform's `std::log`, `std::exp` and `std::cos`; that is the instance exchange of the plan.

## Not measured

macOS: there is no machine, and a shared CI runner has no stable time. The replacement is expected to cost the same
there as on Linux; it has not been measured.
