# Benchmarks

Results of the C++ Google Benchmark suite (`tests/benchmarks/bm_solve.cpp` and `bm_kernel.cpp`). For quick
measurements, A/B comparisons and profiles, see [Measuring](measuring.md). Every number below comes from the one
machine listed here. Compare only runs from the same machine: a number from another machine says nothing about this
library. No other library is measured, so no comparison with one is made.

## Machine and run

| Item | Value |
| --- | --- |
| CPU | AMD Ryzen 7 5700X, 8 cores, 16 hardware threads, 3.4 GHz base (`Get-CimInstance Win32_Processor`) |
| RAM | 32 GiB (31.9 GiB visible to Windows) |
| OS | Windows 10 Pro 10.0.19045 (64-bit) |
| Compiler | g++ 16.2.0 (MSYS2 UCRT64, Rev4), `-O3 -DNDEBUG -ffp-contract=off -std=c++23`, no `-march` (preset `gcc-bench`, Release) |
| CMake | 4.4.2, Ninja |
| Google Benchmark | v1.9.1 |
| Date | 2026-11-16 |
| Code commit | `4cc0b82` (the commit benchmarked; the results were recorded in a later commit) |
| Repetitions | 5 for n = 10 and 100 and for every kernel (median, with the standard deviation); 1 for n = 1000 |

Reproduce (the exact filters of each slice are in [Plan 6 baseline](baseline.md#run)):

```
export PATH=/c/msys64/ucrt64/bin:$PATH
cmake --preset gcc-bench
cmake --build --preset gcc-bench
build/gcc-bench/portopt_benchmarks --benchmark_filter='/(10|100)$|^BM_Kernel/' --benchmark_repetitions=5 --benchmark_report_aggregates_only=true --benchmark_out=build/gcc-bench/bench.json --benchmark_out_format=json
build/gcc-bench/portopt_benchmarks --benchmark_filter='/1000/' --benchmark_repetitions=1 --benchmark_out=build/gcc-bench/bench_1000.json --benchmark_out_format=json
```

The run was not isolated: a normal desktop session, no pinned cores and no fixed frequency. The standard
deviations in the tables are the measured spread over the 5 repetitions:
- median 2.5% of the median time, under 5% for 78 of the 102 benchmarks;
- 10% to 22% for a few cases of a few microseconds at n = 10 (`max_diversification/free`,
  `tracking_error/long_only`, `min_variance/target_return`, `max_return/long_only`, `tracking_error/box`) and for
  `BM_Kernel/gauss_newton/300` and `qp_iteration/1000`.

Differences smaller than that spread are not results. The n = 1000 cases ran once each, so they have no measured
spread: the same cases run four at a time on the same day were up to 1.4 times slower (`max_return/leverage`, 46.0 s
against 32.3 s). Re-measure before and after a change instead of comparing against these tables.

## Not measured

- More than one repetition at n = 1000. With the untimed probe solve each case pays, 5 repetitions of the 35 cases
  take close to an hour on this machine; the baseline ran one.
- `qp_solves`, `descent_iterations` and `polish_steps`, which need a counter in `src/`
  ([Measuring](measuring.md#counters)). The number of calls of a few functions per solve was counted once with a
  debugger instead ([Plan 6 profile 1](profile-1.md#calls-per-solve)).
- The C API and the Python layer (`BM_CApi`, `bench_overhead.py`): plan 6 Task 5.

## Plan 6 baseline

Measured on 2026-11-16 on the machine above, code at `4cc0b82`, before any plan 6 optimization. The full table by
objective, variant and size, the kernels and the commands are in [Plan 6 baseline](baseline.md); where the time goes
is in [Plan 6 profile 1](profile-1.md). Every solve ended `OPTIMAL`.

### Slowest case of each objective

| Objective | n = 10 | n = 100 | n = 1000 |
| --- | ---: | ---: | ---: |
| `min_variance` | 5.96 µs (`box`) | 925 µs (`target_return`) | 932 ms (`box`) |
| `max_return` | 56.2 µs (`leverage`) | 16.6 ms (`leverage`) | 32.3 s (`leverage`) |
| `max_sharpe` | 6.43 µs (`box`) | 1.02 ms (`box`) | 1.02 s (`box`) |
| `utility` | 6.26 µs (`box`) | 1.47 ms (`box`) | 2.38 s (`box`) |
| `tracking_error` | 6.87 µs (`box`) | 2.44 ms (`box`) | 4.06 s (`box`) |
| `max_diversification` | 7.02 µs (`box`) | 940 µs (`box`) | 1.00 s (`box`) |
| `risk_parity` | 252 µs (`target_volatility`) | 73.8 ms (`box`) | 85.3 s (`box`) |

### Fixed cost of a QP solve

The unconstrained QP cases end in one iteration, so they measure what every `qp::solve` pays before it iterates:

| Measure | n = 10 | n = 100 | n = 1000 |
| --- | ---: | ---: | ---: |
| `BM_Solve/min_variance/free` | 2.51 µs | 451 µs | 387 ms |
| `BM_Solve/utility/free` | 2.53 µs | 435 µs | 354 ms |
| `BM_Kernel/empty_active_set` | 1.22 µs | 405 µs | 356 ms |
| `BM_Kernel/cholesky` | 0.196 µs | 69.6 µs | 107 ms |

At n = 100 and 1000 the initial base of the Goldfarb-Idnani solver is about 90% of such a solve, and it costs 3 to 6
Cholesky factorizations of the same size.

## Plan 6 levers of risk 0

Rewrites that keep, for every output element, the sequence of floating-point operations, so every result is bit for
bit the same. Each one is checked three ways: a unit test that compares by `memcmp` with the old code, kept verbatim
in the test (n in 1, 2, 3, 7, 16, 17, 64, 100; condition 1e2 and 1e8); `quick.py compare-solutions
--require-identical` over the profile dump of every case at n = 10 and 100 (70 cases) and of a slice at n = 1000; and
the full `ctest --preset gcc-release`. Machine of [Machine and run](#machine-and-run), 2026-11-16, `gcc-bench`. A and B
are two copies of the executables, run one right after the other in the same session after a discarded warm-up run of
A; the quick suite is the one of [Measuring](measuring.md#quick-suite-about-35-to-45-s) (3 repetitions, medians), the
kernels 5 repetitions.

### L1: base of the Goldfarb-Idnani solver by rows

`empty_active_set` built J0 = L^-T one column at a time, a solve with L' per unit vector: a sum over all n rows with a
stride of n, half of whose terms are `L_ki * 0`. It now computes J0 by rows, from the last up: row i is
`e_i - sum over k > i of L_ki * (row k of J0)` as axpys over the contiguous rows k, divided by L_ii, then the square is
transposed in place into the column storage of the basis. Each entry (i, c) still subtracts its terms in increasing k
and stops at k = c; the terms of k > c that the column solve also subtracted are `L_ki * (+0)`, a zero that leaves the
sum unchanged for a finite factor (the sum starts at +0 or 1 and subtracting finite products never makes it -0, so
subtracting a zero is exact). Below the diagonal
the basis stays +0, which is what `0 / L_ii` gave for the positive diagonal of a Cholesky factor.

| Measure | A (before) | B (after) | B / A |
| --- | ---: | ---: | ---: |
| `BM_Kernel/empty_active_set/10` | 1.19 µs | 0.542 µs | 0.46 |
| `BM_Kernel/empty_active_set/100` | 410 µs | 45.7 µs | 0.11 |
| `BM_Kernel/empty_active_set/300` | 9.94 ms | 1.66 ms | 0.17 |
| `BM_Kernel/empty_active_set/1000` | 367 ms | 45.8 ms | 0.125 |
| `BM_Kernel/qp_iteration/100` | 1565 µs | 1203 µs | 0.77 |
| `BM_Solve/min_variance/long_only/100` | 830 µs | 425 µs | 0.51 |
| `BM_Solve/max_return/long_only/100` | 9.83 ms | 6.10 ms | 0.62 |
| `BM_Solve/max_sharpe/long_only/100` | 915 µs | 511 µs | 0.56 |
| `BM_Solve/utility/box/100` | 1.50 ms | 1.10 ms | 0.73 |
| `BM_Solve/tracking_error/long_only/100` | 2.24 ms | 1.93 ms | 0.86 |
| `BM_Solve/max_diversification/box/100` | 916 µs | 566 µs | 0.62 |
| `BM_Solve/risk_parity/target_volatility/100` | 36.1 ms | 25.5 ms | 0.71 |
| `BM_Solve/risk_parity/box/100` | 71.4 ms | 54.2 ms | 0.76 |
| `BM_Solve/risk_parity/budget/100` (no QP) | 1.09 ms | 1.09 ms | 1.00, noise |

At n = 10 the quick-suite cases moved by 0.77 to 0.97 (`risk_parity/budget/10`, with no QP, 1.03). At n = 1000, one
solve of each `long_only` case with `portopt_profile --max-n 1000 --repeat 1 --filter long_only/1000` (one run each,
no spread): `min_variance` 745 to 519 ms, `max_sharpe` 780 to 498 ms, `max_diversification` 754 to 451 ms, `utility`
1.63 to 1.35 s, `tracking_error` 3.33 to 2.73 s, `max_return` 10.9 to 8.05 s; all 6 bit for bit identical, as are the
70 cases at n ≤ 100.

The kernel gains 8.0 times at n = 1000 (the exploratory 8.4 times of the research) and 9.0 times at n = 100: a third
of the multiply-subtracts (n³/6 instead of n³/2) and unit stride. In the assembly (`goldfarb_idnani.cpp.obj`, inlined
in `solve_factored`) the inner loop is `movupd`, `movupd`, `mulpd`, `subpd`, `movups`: two entries per instruction
and no chain between entries, with GCC's scalar copy of the loop for overlapping spans; the division by L_ii is a
`divpd` loop. The in-place transpose is scalar and strided (n²/2 swaps); its share of the 45.8 ms was not measured
apart.
