# Plan 6 baseline

The time of every case of the benchmark suite before any plan 6 optimization. Where that time goes is in
[Plan 6 profile 1](profile-1.md); how to measure is in [Measuring](measuring.md); the machine and the way the suite is
run are in [Benchmarks](benchmarks.md#machine-and-run). Compare a later run against this page only after measuring the
baseline again in the same session ([Measuring](measuring.md#ab)): on this machine medians move between days by more
than their spread.

## Run

Preset `gcc-bench` (generic x86-64, `-O3 -DNDEBUG -ffp-contract=off`), code at `4cc0b82`, 2026-11-16. Desktop
session, not isolated. Four invocations, one right after the other:

```
build/gcc-bench/portopt_benchmarks --benchmark_filter='^BM_Kernel/' --benchmark_repetitions=5 --benchmark_report_aggregates_only=true --benchmark_out=build/gcc-bench/baseline_kernels.json --benchmark_out_format=json
build/gcc-bench/portopt_benchmarks --benchmark_filter='^BM_Solve/.*/(10|100)$' --benchmark_repetitions=5 --benchmark_report_aggregates_only=true --benchmark_out=build/gcc-bench/baseline_small.json --benchmark_out_format=json
build/gcc-bench/portopt_benchmarks --benchmark_filter='^BM_Solve/risk_parity/.*/1000/' --benchmark_repetitions=1 --benchmark_out=build/gcc-bench/baseline_1000_rp.json --benchmark_out_format=json
build/gcc-bench/portopt_benchmarks --benchmark_filter='^BM_Solve/(min_variance|max_return|max_sharpe|utility|tracking_error|max_diversification)/.*/1000/' --benchmark_repetitions=1 --benchmark_out=build/gcc-bench/baseline_1000_rest.json --benchmark_out_format=json
```

- n = 10 and 100, and every kernel: **5 repetitions**, median and standard deviation of the wall time. The kernels took
  2 min 27 s, the 70 small cases 4 min 17 s.
- n = 1000: **one repetition, one timed solve per case** (`Iterations(1)`), so no spread is measured. Each case also
  pays an untimed probe solve that checks the status first, so the 35 cases took 9 min 42 s of wall time for about
  4 min 50 s of timed solves. Five repetitions would take close to an hour; one is the budget of this baseline. Read a
  difference of less than about 10% at n = 1000 as noise: Task 0's probe of the same cases, four processes at a time,
  gave 46.0 s for `max_return/leverage`, against 32.3 s here alone.
- Every case ended `OPTIMAL` (`status` 0) on every solve.

Spread at n = 10 and 100 and in the kernels: median relative standard deviation 2.5%, 78 of 102 benchmarks under
5%. Over 10%: `max_diversification/free/10` (22%), `tracking_error/long_only/10` (15%),
`min_variance/target_return/10` (13%), `max_return/long_only/10` (11%), `tracking_error/box/10` (11%),
`qp_iteration/1000` (10%) and `gauss_newton/300` (22%). Most of them are cases of a few microseconds.

## Solves

`BM_Solve/<objective>/<variant>/<n>`: median wall time, its standard deviation (sd), `solver_iterations` (it) and
`allocations` per solve (alloc). The instances are in [Measuring](measuring.md#the-instances).

| Case | n = 10 | sd | it | alloc | n = 100 | sd | it | alloc | n = 1000 (1 run) | it | alloc |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `min_variance/free` | 2.51 µs | 0.02 µs | 1 | 20 | 451 µs | 39 µs | 1 | 20 | 387 ms | 1 | 20 |
| `min_variance/long_only` | 4.77 µs | 0.11 µs | 3 | 52 | 756 µs | 7.0 µs | 21 | 158 | 781 ms | 229 | 1073 |
| `min_variance/box` | 5.96 µs | 0.08 µs | 4 | 66 | 869 µs | 5.8 µs | 21 | 261 | 932 ms | 229 | 2076 |
| `min_variance/leverage` | 4.41 µs | 0.05 µs | 3 | 52 | 442 µs | 14 µs | 1 | 25 | 345 ms | 1 | 25 |
| `min_variance/cash` | 3.45 µs | 0.18 µs | 0 | 42 | 456 µs | 3.4 µs | 0 | 141 | 354 ms | 0 | 1050 |
| `min_variance/target_return` | 5.26 µs | 0.69 µs | 3 | 54 | 925 µs | 15 µs | 33 | 162 | 919 ms | 269 | 1078 |
| `max_return/free` | 8.49 µs | 0.47 µs | 3 | 59 | 1.74 ms | 16 µs | 4 | 85 | 1.38 s | 4 | 85 |
| `max_return/long_only` | 44.3 µs | 5.1 µs | 41 | 263 | 10.2 ms | 223 µs | 358 | 457 | 11.4 s | 4174 | 1466 |
| `max_return/box` | 55.5 µs | 1.4 µs | 70 | 303 | 11.3 ms | 61 µs | 358 | 560 | 14.4 s | 4174 | 2469 |
| `max_return/leverage` | 56.2 µs | 0.86 µs | 53 | 709 | 16.6 ms | 137 µs | 646 | 8984 | 32.3 s | 9653 | 183690 |
| `max_return/cash` | 41.0 µs | 0.65 µs | 42 | 282 | 10.2 ms | 71 µs | 365 | 509 | 11.5 s | 3803 | 1724 |
| `max_sharpe/free` | 3.55 µs | 0.07 µs | 1 | 39 | 469 µs | 45 µs | 1 | 45 | 356 ms | 1 | 51 |
| `max_sharpe/long_only` | 6.00 µs | 0.34 µs | 4 | 60 | 888 µs | 9.1 µs | 29 | 168 | 820 ms | 255 | 1086 |
| `max_sharpe/box` | 6.43 µs | 0.31 µs | 4 | 68 | 1.02 ms | 50 µs | 29 | 263 | 1.02 s | 255 | 2078 |
| `max_sharpe/leverage` | 4.98 µs | 0.25 µs | 3 | 59 | 437 µs | 6.8 µs | 1 | 32 | 359 ms | 1 | 32 |
| `max_sharpe/cash` | 5.54 µs | 0.11 µs | 4 | 61 | 852 µs | 58 µs | 27 | 169 | 816 ms | 263 | 1089 |
| `utility/free` | 2.53 µs | 0.09 µs | 1 | 20 | 435 µs | 7.5 µs | 1 | 20 | 354 ms | 1 | 20 |
| `utility/long_only` | 5.06 µs | 0.04 µs | 4 | 53 | 1.22 ms | 27 µs | 51 | 160 | 1.69 s | 766 | 1077 |
| `utility/box` | 6.26 µs | 0.55 µs | 4 | 66 | 1.47 ms | 60 µs | 51 | 263 | 2.38 s | 766 | 2080 |
| `utility/leverage` | 4.40 µs | 0.28 µs | 3 | 52 | 1.32 ms | 55 µs | 71 | 1016 | 2.18 s | 855 | 17568 |
| `utility/cash` | 4.47 µs | 0.14 µs | 3 | 49 | 1.20 ms | 7.5 µs | 51 | 156 | 1.76 s | 766 | 1073 |
| `tracking_error/free` | 4.65 µs | 0.11 µs | 1 | 24 | 1.47 ms | 35 µs | 1 | 24 | 2.02 s | 1 | 24 |
| `tracking_error/long_only` | 6.65 µs | 0.98 µs | 1 | 52 | 2.13 ms | 47 µs | 45 | 164 | 3.25 s | 720 | 1081 |
| `tracking_error/box` | 6.87 µs | 0.76 µs | 1 | 65 | 2.44 ms | 72 µs | 45 | 267 | 4.06 s | 720 | 2084 |
| `tracking_error/leverage` | 4.69 µs | 0.05 µs | 1 | 29 | 1.94 ms | 67 µs | 46 | 598 | 3.25 s | 541 | 11706 |
| `tracking_error/cash` | 5.39 µs | 0.04 µs | 0 | 46 | 2.15 ms | 48 µs | 45 | 160 | 3.60 s | 720 | 1077 |
| `max_diversification/free` | 4.04 µs | 0.88 µs | 1 | 38 | 453 µs | 30 µs | 1 | 44 | 390 ms | 1 | 50 |
| `max_diversification/long_only` | 5.88 µs | 0.29 µs | 4 | 59 | 795 µs | 36 µs | 22 | 167 | 807 ms | 233 | 1085 |
| `max_diversification/box` | 7.02 µs | 0.27 µs | 5 | 68 | 940 µs | 48 µs | 22 | 262 | 1.00 s | 233 | 2077 |
| `max_diversification/leverage` | 5.64 µs | 0.15 µs | 4 | 75 | 442 µs | 7.4 µs | 1 | 31 | 378 ms | 1 | 31 |
| `max_diversification/cash` | 5.48 µs | 0.11 µs | 4 | 60 | 797 µs | 69 µs | 22 | 168 | 803 ms | 233 | 1086 |
| `risk_parity/budget` | 5.39 µs | 0.06 µs | 7 | 49 | 1.09 ms | 20 µs | 13 | 85 | 1.58 s | 14 | 91 |
| `risk_parity/target_return` | 98.5 µs | 2.1 µs | 14 | 706 | 34.3 ms | 843 µs | 37 | 1946 | 68.2 s | 80 | 6986 |
| `risk_parity/target_volatility` | 252 µs | 6.9 µs | 35 | 1566 | 34.1 ms | 832 µs | 39 | 1566 | 26.4 s | 38 | 1541 |
| `risk_parity/box` | 149 µs | 14 µs | 19 | 939 | 73.8 ms | 5.3 ms | 56 | 3051 | 85.3 s | 75 | 9686 |

What the table shows:
- **A QP solve at n = 1000 costs at least about 0.35 s whatever it does.** `min_variance/free`, `utility/free`,
  `max_sharpe/free` and the `leverage` cases that end in one iteration take 345 to 390 ms; `BM_Kernel/empty_active_set`
  alone takes 356 ms. At n = 100 the same cases take 435 to 470 µs against 405 µs for the kernel.
- **Three families take almost all of the n = 1000 time:** `risk_parity` with constraints (180 s for three cases),
  `max_return` with any constraint (70 s for four cases) and `tracking_error` (16 s for five). The other 23 cases
  take 22 s together.
- From n = 100 to 1000 the constrained cases grow by 770 to 2000 times (n³ alone is 1000): the iteration counts grow
  about 10 times with n (`max_return/leverage`: 646 to 9653) and each iteration costs O(n²).
- Allocations: at most 1566 per solve at n = 10 (`risk_parity/target_volatility`), and 183 690 per solve for
  `max_return/leverage` at n = 1000 (19 per path iteration).

## Kernels

`BM_Kernel/<kernel>/<n>`, median of 5 repetitions (see [Measuring](measuring.md#counters) for what each one runs).

| Kernel | n = 10 | n = 100 | n = 300 | n = 1000 | sd at 1000 |
| --- | ---: | ---: | ---: | ---: | ---: |
| `cholesky` | 0.196 µs | 69.6 µs | 2.38 ms | 107 ms | 1.8 ms |
| `solve_lower` | 0.068 µs | 2.25 µs | 25.9 µs | 318 µs | 1.7 µs |
| `qr` (2n × n) | 1.13 µs | 774 µs | 25.0 ms | 1.20 s | 109 ms |
| `qp_iteration` (one long-only QP) | 3.37 µs | 1.56 ms | 45.5 ms | 1.85 s | 187 ms |
| `compensated_dot` | 0.043 µs | 0.401 µs | 1.19 µs | 3.95 µs | 0.11 µs |
| `covariance_times` | 0.441 µs | 39.4 µs | 354 µs | 3.94 ms | 59 µs |
| `gauss_newton` | 0.270 µs | 152 µs | 6.20 ms | 191 ms | 12 ms |
| `empty_active_set` | 1.22 µs | 405 µs | 9.94 ms | 356 ms | 15 ms |

`qp_iteration` took 4, 75, 243 and 895 iterations: 0.84 µs, 20.5 µs, 189 µs and 2.08 ms per iteration, the base
included. Per element: `compensated_dot` 3.95 ns at n = 1000 (a `call fma` per term, [Plan 6 profile
1](profile-1.md#assembly)); `empty_active_set` costs 5.8 Choleskys at n = 100 and 3.3 at n = 1000; `qr` of a 2000 ×
1000 matrix costs 11 Choleskys of the same n.
