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
