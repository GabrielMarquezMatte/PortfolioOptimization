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
  `&run_solve_c` and `&run_solve_cpp`; the repository has that fix since (`fix(bench): take the address of the C API
  runners for GCC 14`).

Commands, per platform (after one discarded run of each binary): four rounds of
`python tests/benchmarks/quick.py run --exe EXE --out FILE` with A and B alternating, A first in every round.

A single A then B pair, the procedure of Measuring, was taken first and was not enough. On this machine two runs of
the same binary moved some medians by 20 to 30% (A against A again), and one B run hit a disturbance that doubled a
third of the cases. On Linux the single pair put B up to 30% above A on kernels that do not call the replaced
functions. The tables below are therefore the minimum over the four rounds of each median (A min, B min), the ratio
of those minima, the ratio of the medians over the rounds, how many rounds B was not slower than A, and the
standard deviation of the four round medians of each binary (A sd, B sd).

A case is worse than noise only if B is slower in both ratios and the gap of the minima, B min minus A min, exceeds
the larger of the two standard deviations. The last column applies that rule. A difference under about 5% is not
resolved by this method: four rounds of three repetitions are a screen, not a test.

## Windows

| Case | A min (us) | B min (us) | B/A of the minima | B/A of the medians | rounds with B <= A | A sd | B sd | Worse than noise |
|---|---|---|---|---|---|---|---|---|
| `min_variance/long_only/10` | 3.324 | 3.021 | 0.909 | 0.751 | 4/4 | 1.08 | 0.738 | no |
| `min_variance/long_only/100` | 227.5 | 235.1 | 1.033 | 0.727 | 2/4 | 50.8 | 48.4 | no |
| `max_return/long_only/10` | 30.92 | 28.84 | 0.933 | 0.952 | 4/4 | 5.29 | 3.03 | no |
| `max_return/long_only/100` | 2647 | 2517 | 0.951 | 0.960 | 2/4 | 173 | 633 | no |
| `max_sharpe/long_only/10` | 3.958 | 3.681 | 0.930 | 0.929 | 3/4 | 0.23 | 0.389 | no |
| `max_sharpe/long_only/100` | 277.3 | 250.1 | 0.902 | 1.092 | 1/4 | 16.2 | 37.7 | no |
| `utility/box/10` | 3.939 | 3.787 | 0.961 | 0.963 | 3/4 | 0.047 | 0.202 | no |
| `utility/box/100` | 456.4 | 444 | 0.973 | 1.036 | 1/4 | 10.1 | 41.2 | no |
| `tracking_error/long_only/10` | 4.465 | 4.65 | 1.042 | 1.076 | 0/4 | 0.101 | 0.175 | yes |
| `tracking_error/long_only/100` | 980.8 | 906.4 | 0.924 | 0.951 | 3/4 | 112 | 69.2 | no |
| `max_diversification/box/10` | 4.806 | 4.339 | 0.903 | 0.962 | 3/4 | 0.374 | 0.242 | no |
| `max_diversification/box/100` | 383 | 380.8 | 0.994 | 0.935 | 3/4 | 44.2 | 14.4 | no |
| `risk_parity/budget/10` | 6.164 | 6.156 | 0.999 | 1.035 | 1/4 | 0.0853 | 0.35 | no |
| `risk_parity/budget/100` | 639.7 | 621.8 | 0.972 | 1.048 | 1/4 | 116 | 48.5 | no |
| `risk_parity/target_volatility/10` | 231.6 | 236.4 | 1.021 | 0.951 | 2/4 | 46.7 | 37.5 | no |
| `risk_parity/target_volatility/100` | 1.555e+04 | 1.548e+04 | 0.996 | 1.040 | 1/4 | 3.14e+03 | 958 | no |
| `risk_parity/box/10` | 139.5 | 135.4 | 0.971 | 0.961 | 3/4 | 26.5 | 14.5 | no |
| `risk_parity/box/100` | 3.138e+04 | 2.867e+04 | 0.914 | 0.833 | 4/4 | 3.2e+03 | 750 | no |
| `kernel/cholesky/100` | 34.53 | 30.82 | 0.893 | 0.920 | 4/4 | 2.02 | 1.42 | no |
| `kernel/solve_lower/100` | 2.247 | 2.256 | 1.004 | 0.943 | 2/4 | 0.111 | 0.208 | no |
| `kernel/qr/100` | 335.2 | 348 | 1.038 | 1.058 | 1/4 | 130 | 28 | no |
| `kernel/qp_iteration/100` | 534.3 | 509.3 | 0.953 | 1.014 | 3/4 | 186 | 89 | no |
| `kernel/compensated_dot/100` | 0.3098 | 0.3142 | 1.014 | 1.034 | 2/4 | 0.00921 | 0.134 | no |
| `kernel/covariance_times/100` | 30.59 | 31.44 | 1.028 | 1.034 | 2/4 | 1.19 | 17.2 | no |
| `kernel/gauss_newton/100` | 65.15 | 69.21 | 1.062 | 1.016 | 2/4 | 9.4 | 29 | no |
| `kernel/empty_active_set/100` | 44.93 | 45.39 | 1.010 | 1.049 | 2/4 | 4.6 | 9.11 | no |

The one case that tripped the rule, `tracking_error/long_only/10` (B 4% slower, gap 0.185 us against an sd of
0.175 us, B not slower in 0 of 4 rounds), was repeated with six alternating rounds of 15 repetitions. B was not slower
in five of six rounds (round medians 4.266 against 4.302 us at the best round), and A alone ranges from 4.3 to 9.3 us
across those rounds, so the 4% is inside the noise of the machine.

The `kernel/*` rows where B is slower (`solve_lower`, `qr`, `compensated_dot`, `covariance_times`, `gauss_newton`,
`empty_active_set`, by 0.4 to 6%) call neither `hypot` nor `log`; they only show the spread of the machine. The
repeat below confirms it.

Repeat of the kernels at n = 100, alternating A and B, six rounds of nine repetitions (54 per binary), the median
and the minimum over all repetitions and the standard deviation of the repetitions:

| Case | reps per binary | A median | B median | B/A median | A min | B min | B/A min | A stddev | B stddev |
|---|---|---|---|---|---|---|---|---|---|
| `kernel/cholesky/100` | 54 | 31.54 | 31.74 | 1.006 | 30.16 | 30.37 | 1.007 | 2.38 | 1.3 |
| `kernel/solve_lower/100` | 54 | 2.194 | 2.186 | 0.997 | 2.162 | 2.164 | 1.001 | 0.0443 | 0.157 |
| `kernel/qr/100` | 54 | 306.4 | 312.6 | 1.020 | 296.9 | 300.9 | 1.014 | 12.5 | 33.9 |
| `kernel/qp_iteration/100` | 54 | 515.4 | 485.7 | 0.942 | 473.3 | 458.8 | 0.969 | 29.6 | 47.6 |
| `kernel/compensated_dot/100` | 54 | 0.3148 | 0.3077 | 0.978 | 0.3009 | 0.3009 | 1.000 | 0.0437 | 0.019 |
| `kernel/covariance_times/100` | 54 | 30.85 | 30.66 | 0.994 | 29.95 | 30 | 1.001 | 2.86 | 2.09 |
| `kernel/gauss_newton/100` | 54 | 66.54 | 66.88 | 1.005 | 65.07 | 65.2 | 1.002 | 4.94 | 7.15 |
| `kernel/empty_active_set/100` | 54 | 43.75 | 43.97 | 1.005 | 42.05 | 42.23 | 1.004 | 2.5 | 3.62 |

## Linux

| Case | A min (us) | B min (us) | B/A of the minima | B/A of the medians | rounds with B <= A | A sd | B sd | Worse than noise |
|---|---|---|---|---|---|---|---|---|
| `kernel/cholesky/100` | 29.16 | 30.5 | 1.046 | 1.064 | 2/4 | 7.91 | 3.57 | no |
| `kernel/solve_lower/100` | 2.266 | 2.706 | 1.194 | 1.174 | 1/4 | 0.428 | 0.259 | yes |
| `kernel/qr/100` | 401.1 | 516 | 1.287 | 1.340 | 1/4 | 78.4 | 62.3 | yes |
| `kernel/qp_iteration/100` | 497 | 480.1 | 0.966 | 1.048 | 2/4 | 77 | 47.1 | no |
| `kernel/compensated_dot/100` | 0.3229 | 0.3035 | 0.940 | 1.020 | 2/4 | 0.0198 | 0.0206 | no |
| `kernel/covariance_times/100` | 30.21 | 29.48 | 0.976 | 0.806 | 4/4 | 3.57 | 0.762 | no |
| `kernel/gauss_newton/100` | 67.44 | 66.7 | 0.989 | 1.004 | 2/4 | 4.87 | 2.44 | no |
| `kernel/empty_active_set/100` | 44.65 | 46.5 | 1.041 | 1.045 | 1/4 | 4.44 | 2.66 | no |
| `min_variance/long_only/10` | 2.697 | 2.498 | 0.926 | 0.974 | 3/4 | 0.285 | 0.647 | no |
| `min_variance/long_only/100` | 418.4 | 387.9 | 0.927 | 0.872 | 3/4 | 25.2 | 41.7 | no |
| `max_return/long_only/10` | 25.75 | 23.9 | 0.928 | 1.007 | 2/4 | 1.59 | 2.99 | no |
| `max_return/long_only/100` | 2527 | 2414 | 0.955 | 0.951 | 3/4 | 279 | 243 | no |
| `max_sharpe/long_only/10` | 3.281 | 3.042 | 0.927 | 0.955 | 3/4 | 0.353 | 0.297 | no |
| `max_sharpe/long_only/100` | 288.7 | 261.5 | 0.906 | 0.970 | 3/4 | 25.5 | 30.5 | no |
| `utility/box/10` | 3.31 | 3.098 | 0.936 | 0.881 | 4/4 | 0.196 | 0.198 | no |
| `utility/box/100` | 404.9 | 392.5 | 0.969 | 0.901 | 4/4 | 24.9 | 18.6 | no |
| `tracking_error/long_only/10` | 3.803 | 3.816 | 1.003 | 0.923 | 3/4 | 0.389 | 0.206 | no |
| `tracking_error/long_only/100` | 1046 | 1066 | 1.019 | 1.040 | 1/4 | 42.1 | 265 | no |
| `max_diversification/box/10` | 3.983 | 3.885 | 0.975 | 0.997 | 3/4 | 0.124 | 0.506 | no |
| `max_diversification/box/100` | 424.1 | 414.5 | 0.977 | 1.054 | 1/4 | 23.7 | 53.1 | no |
| `risk_parity/budget/10` | 4.651 | 4.383 | 0.942 | 1.058 | 1/4 | 0.163 | 0.354 | no |
| `risk_parity/budget/100` | 594.8 | 572.6 | 0.963 | 1.024 | 2/4 | 22.9 | 48.1 | no |
| `risk_parity/target_volatility/10` | 176.9 | 171.6 | 0.970 | 1.118 | 1/4 | 4.07 | 23.1 | no |
| `risk_parity/target_volatility/100` | 1.423e+04 | 1.421e+04 | 0.999 | 1.039 | 2/4 | 1.42e+03 | 1.15e+03 | no |
| `risk_parity/box/10` | 99.39 | 96.84 | 0.974 | 0.867 | 4/4 | 7.49 | 8.01 | no |
| `risk_parity/box/100` | 3.197e+04 | 2.887e+04 | 0.903 | 0.911 | 4/4 | 4.89e+03 | 1.98e+03 | no |

`kernel/solve_lower/100` and `kernel/qr/100` tripped the rule in the four rounds (B 19% and 29% above A), although
neither calls `hypot` or `log`, and `solve_lower` is the same machine code in A and B. They were repeated, with
`kernel/cholesky/100` as a control, alternating A and B, eight rounds of nine repetitions (72 per binary):

| Case | reps per binary | A median | B median | B/A median | A min | B min | B/A min | A stddev | B stddev |
|---|---|---|---|---|---|---|---|---|---|
| `kernel/cholesky/100` | 72 | 29.17 | 28.7 | 0.984 | 27.87 | 27.77 | 0.996 | 2.4 | 1.94 |
| `kernel/solve_lower/100` | 72 | 2.227 | 2.272 | 1.020 | 2.163 | 2.157 | 0.997 | 0.132 | 0.188 |
| `kernel/qr/100` | 72 | 406.2 | 413.1 | 1.017 | 388.2 | 381.7 | 0.983 | 41.7 | 48.5 |

B is no slower: the ratio of the minima is 0.997 and 0.983, and that of the medians 1.020 and 1.017, inside the
standard deviation of the repetitions (6 to 8% of the median for `solve_lower`, 10 to 12% for `qr`). The first result was the machine, not the change.
Every `BM_Solve` case is within noise or faster.

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
