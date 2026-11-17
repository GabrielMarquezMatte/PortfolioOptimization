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

### L2: right-looking Cholesky, two columns at a time

`factor_in_place` was left-looking: entry (i, j) was one dot product of rows i and j, a chain of dependent
subtractions (`mulpd` and then two serial `subsd`). It is now right-looking, two columns at a time: once columns k and
k + 1 are final, every entry (i, j) with j ≥ k + 2 subtracts `L_ik L_jk` and then `L_i,k+1 L_j,k+1`, entries
independent of each other. Each entry still starts from a_ij, subtracts its terms in increasing k and is divided by
L_jj (or square-rooted, after the same pivot test) after its last one. Column k of L is staged in row k above the
diagonal while it updates the block, so that the update reads it with unit stride; the row is cleared after, so no
buffer is allocated and the upper triangle is zero at the end, as before. The factor starts as a copy of the whole
square a (one `memcpy`) instead of zeros and a copy of the lower triangle: the part above the diagonal is overwritten
by the staging before it is read.

| Measure | A (L1) | B (L1 + L2) | B / A |
| --- | ---: | ---: | ---: |
| `BM_Kernel/cholesky/10` | 0.193 µs | 0.227 µs | 1.18 |
| `BM_Kernel/cholesky/100` | 69.4 µs | 30.9 µs | 0.45 |
| `BM_Kernel/cholesky/300` | 2.35 ms | 0.764 ms | 0.33 |
| `BM_Kernel/cholesky/1000` | 108 ms | 27.9 ms | 0.26 |
| `BM_Solve/risk_parity/budget/100` (two runs) | 1059 and 1181 µs | 638 and 655 µs | 0.60 and 0.55 |
| `BM_Solve/risk_parity/budget/1000` (one solve, `portopt_profile`) | 1.54 s | 0.497 s | 0.32 |

The other cases of the quick suite do not resolve it: two A/B pairs of the same two executables gave B/A from 0.66 to
1.35 on cases that do not factor anything inside the solve (`min_variance/long_only/100` 1.03 and then 1.35) and
kernels that did not change (`qp_iteration/100` 1.33 and then 1.03): the desktop was noisy during these runs.
Cholesky was 5 to 9% of `risk_parity` at n = 100 ([Plan 6 profile 1](profile-1.md)), so 2.2 times on it is 3 to 5%
of those cases, under what the quick suite resolves; `risk_parity/target_volatility/100` and `risk_parity/box/100`
gave 1.01 and 0.93, 1.10 and 0.94.

A plain right-looking loop (one column at a time, with a copy of the lower triangle) gave 2.3 times at n = 1000
(106.8 to 46.5 ms) and 1.75 times at n = 100 (68.6 to 39.1 µs); pairing the columns halves the passes over the
trailing block and gave 3.9 and 2.2 times. At n = 10 both are slower than the left-looking form (0.26 µs for the
plain loop, 0.23 µs for the pairs with the `memcpy` start, against 0.19 µs): the loops are a few entries long and the
vectorized ones carry their alias checks and remainders; 34 ns a factorization, under the spread of the n = 10 solves.

All 70 cases at n ≤ 100, the 6 `long_only` cases at n = 1000 (whose covariance factor of `Data::create` is now
computed this way) and `risk_parity/budget/1000` are bit for bit identical. In the assembly (`cholesky.cpp.obj`) the
update loop is `movupd` × 3, `mulpd`, `subpd`, `mulpd`, `subpd`, `movups`: two entries per instruction, the two
subtractions of each entry in order; the square root stays an inline `sqrtsd` (the `call sqrt` is the errno path of a
negative argument, which the pivot test excludes).

### L3: Householder reflections by rows

`householder_qr` reflected each trailing column j on its own (`reflect` with a stride of cols): a projection summed
down the column and an update down the column, both reading one entry per row of a row-major matrix. It now reflects
all trailing columns of reflector k together, a row at a time (`reflect_trailing_columns`): a vector of cols
projections, each summed from zero over the rows i ≥ k in increasing i, then doubled, then `x_ij -= (2 p_j) v_i` row
by row. Column j sees the same operations in the same order as before (`2.0 * projection * v[i]` is `(2.0 *
projection) * v[i]`). One vector of cols doubles is allocated per factorization. `reflect`, still used by
`apply_transposed` on a contiguous vector, lost its stride parameter, which was always 1 there.

`BM_Kernel/qr` factors a 2n × n matrix. The A executable (L1 + L2) was run three times, because its strided loop moved
by up to 40% between runs on the noisy desktop (the baseline table has 1.20 s at n = 1000):

| Measure | A (L1 + L2) | B (L1 + L2 + L3) | B / A |
| --- | ---: | ---: | ---: |
| `BM_Kernel/qr/10` | 1.18 to 1.30 µs | 1.18 to 1.20 µs | 0.92 to 1.0 |
| `BM_Kernel/qr/100` | 874 to 1140 µs | 310 and 353 µs | 0.27 to 0.40 |
| `BM_Kernel/qr/300` | 25.3 to 27.4 ms | 10.8 and 10.9 ms | 0.39 to 0.43 |
| `BM_Kernel/qr/1000` | 1.34 to 1.88 s | 0.381 and 0.403 s | 0.20 to 0.30 |
| `BM_Solve/tracking_error/long_only/100` (two A/B pairs) | 1904 and 1928 µs | 1492 and 1262 µs | 0.78 and 0.65 |
| `BM_Solve/tracking_error/box/1000` (one solve, `portopt_profile`) | 3.65 s | 2.71 s | 0.74 |

The other cases of the quick suite, which do not factor by QR in the solve or do it on a small system (`risk_parity`),
moved between 0.74 and 1.63 in the first pair and between 0.84 and 1.17 in the second, in both directions: noise of
the desktop, not of the change. All 70 cases at n ≤ 100, the 6 `long_only` cases at n = 1000 and
`tracking_error/box/1000` are bit for bit identical.

In the assembly (`qr.cpp.obj`, inlined in `householder_qr`) GCC unrolls the row loop of the projections by two and
jams it: one pass loads the projections once and adds `v_i x_ij` and then `v_i+1 x_i+1,j` (`mulpd`, `addpd`, `mulpd`,
`addpd`), in the order of i, which is legal without reassociation. The doubling is `addpd x, x`, exact like `2.0 *
x`; the update is `movupd`, `movupd`, `mulpd`, `subpd`, `movups`. `column_norm` and the copy of column k in
`make_reflector` still read with a stride; they are O(rows · cols) per factorization against the O(rows · cols²) of the
reflections.

### The three together

`risk_parity/target_volatility/1000`, one solve each with `portopt_profile --max-n 1000 --repeat 1`, the executables of
`1297db0` (before L1) and of L1 + L2 + L3, one right after the other: 25.8 s and 12.7 s (0.49), bit for bit identical
(status, objective, weights, multipliers and residual), with 7 more allocations (the projections of its 7 QR
factorizations).

### L5: unit rows (bounds) read at their one entry

A bound row `x_j >= l` or `-x_j >= -u` has one nonzero entry, but the solver read it as a dense row of n doubles: the
slack of every inequality row in the scan for the most violated one (`dot(row, x)`), and the n projections of the
basis on the candidate's normal in `compute_directions`. `gather_constraints` now records, for every explicit row,
the column of its only nonzero entry (or n when it has more than one). The scan reads such a row as
`r_j x_j + 0.0`, and `compute_unit_directions` projects with `sign * (J_jc r_j + 0.0)`, one entry of each basis
column; the n − q axpys of the primal direction and the back substitution are unchanged. The skipped terms are
`a_i * 0`, zeros that leave `dot()`'s sum unchanged when `a_i` is finite (the sum starts at +0, adding a zero to +0 or
to a nonzero value is exact), and the `+ 0.0` gives the +0 that `dot()` returns when the one product is -0. The scan
checks that x is finite (one pass over x) and falls back to the dense dot otherwise. The projections assume a finite
basis, as L1 assumes a finite factor: J J' = G^-1 and the Givens rotations keep it bounded, but a basis entry that
overflowed would give NaN in the dense projection and not in the unit one. One more allocation per `qp::solve` (the
column of each row).

`ActiveSet.UnitRowDirectionsAreBitIdenticalToTheDenseRow` compares projection, primal and dual directions by
`memcmp` against `compute_directions` on the dense row (with -0 in a third of its zeros), for every column, entries
1, -1 and random, signs ±1, on J0 and on bases rotated by n/2 and n − 1 random rows (n and condition as above);
`ActiveSet.UnitRowDotIsBitIdenticalToTheDenseDot` does the same for the slack against x with +0 and -0 entries.
Dropping the `+ 0.0` makes both fail.

| Measure | A (L1 + L2 + L3) | B (+ L5) | B / A |
| --- | ---: | ---: | ---: |
| `BM_Kernel/qp_iteration/10` (long-only QP, 5 repetitions) | 2.71 µs | 2.84 µs | 1.05, noise |
| `BM_Kernel/qp_iteration/100` | 1219 µs | 559 µs | 0.46 |
| `BM_Kernel/qp_iteration/300` | 38.3 ms | 11.2 ms | 0.29 |
| `BM_Kernel/qp_iteration/1000` | 1.63 s | 0.381 s | 0.23 |
| `BM_Solve/min_variance/long_only/100` | 569 µs | 311 µs | 0.55 |
| `BM_Solve/max_return/long_only/100` | 6.19 ms | 3.04 ms | 0.49 |
| `BM_Solve/max_sharpe/long_only/100` | 533 µs | 315 µs | 0.59 |
| `BM_Solve/utility/box/100` | 1186 µs | 478 µs | 0.40 |
| `BM_Solve/tracking_error/long_only/100` | 1377 µs | 943 µs | 0.68 |
| `BM_Solve/max_diversification/box/100` | 564 µs | 398 µs | 0.71 |
| `BM_Solve/risk_parity/box/100` | 50.8 ms | 37.0 ms | 0.73 |
| `BM_Solve/risk_parity/target_volatility/100` (no active bound) | 21.3 ms | 21.1 ms | 0.99, noise |

The n = 10 cases of the quick suite moved by 0.93 to 1.01 and `risk_parity/budget` (no QP) by 0.96 and 0.98, noise.
At n = 1000, one solve each with `portopt_profile --max-n 1000 --repeat 1` (one run, no spread):

| Case (n = 1000) | A | B | B / A |
| --- | ---: | ---: | ---: |
| `min_variance/long_only` | 447 ms | 153 ms | 0.34 |
| `max_return/long_only` | 7.92 s | 2.26 s | 0.28 |
| `max_sharpe/long_only` | 495 ms | 166 ms | 0.34 |
| `utility/long_only` | 1.34 s | 0.336 s | 0.25 |
| `tracking_error/long_only` | 1.94 s | 0.915 s | 0.47 |
| `max_diversification/long_only` | 464 ms | 161 ms | 0.35 |
| `min_variance/box` | 626 ms | 172 ms | 0.27 |
| `max_return/box` | 11.3 s | 2.59 s | 0.23 |
| `max_sharpe/box` | 680 ms | 377 ms | 0.55 |
| `utility/box` | 1.91 s | 0.383 s | 0.20 |
| `tracking_error/box` | 2.42 s | 0.977 s | 0.40 |
| `max_diversification/box` | 632 ms | 344 ms | 0.54 |
| `risk_parity/box` | 63.3 s | 28.3 s | 0.45 |
| `risk_parity/target_return` | 42.7 s | 27.4 s | 0.64 |

All 70 cases at n ≤ 100 and the 14 cases at n = 1000 above are bit for bit identical, with the same iterations. In
the assembly (`goldfarb_idnani.cpp.obj`, inlined in `solve_factored` and `most_violated_inequality`) the projection
loop is `movsd`/`movhpd` (two basis columns, a stride of n), `mulpd` by the entry, `addpd` of a zero register,
`mulpd` by the sign, `movups`: the `+ 0.0` is kept, as it must be without `-ffast-math`. The scan of a unit row is
`movsd`, `mulsd`, `addsd` of zero, `subsd` of the right-hand side, then its tolerance as before.

### L4: `gauss_newton` by symmetry, in blocks

`gauss_newton` added each row i of J to all n² entries of 2 J'J (`curvature[j][k] += (2 J_ij) J_ik`): it read and
wrote the whole n² curvature n times, 8 MB a pass at n = 1000. It now computes the blocks of 4 × 4 entries on and
above the diagonal and mirrors each one below it. J is first copied into panels of 4 columns (n rows of 4 entries,
zeros past column n), and the 4 columns of the block's rows are doubled once per panel (`2.0 * J_ij`, as before).
Each block keeps its 16 sums in registers over the whole loop in i: entry (j, k) is still summed from +0 in increasing
i, with the same product `(2 J_ij) J_ik`. The entry (k, j) below the diagonal gets the bits of (j, k): `2 J_ij` is
exact, so `(2 J_ij) J_ik` and `(2 J_ik) J_ij` are the same rounding of the same real number, overflow and subnormals
included, and the two sums add the same terms in the same order. This holds while every |J_ij| ≤ max / 2: one pass
over J checks it, and otherwise (an entry above it, inf or NaN) the old loop runs. The old loop also runs below
n = 16, where the packing costs more than the half of the products it saves (a scratch timing of both on a random J,
best of 7: 1.57 times slower at n = 10, 1.25 at n = 12, 0.92 at n = 16, 0.64 at n = 32). Two more allocations per
call at n ≥ 16 (the panels and the doubled panel).

`GaussNewton.SymmetricIsBitIdenticalToTheFullLoop` compares by `memcmp` with the old loop, kept verbatim in the test,
on the Jacobian of `linearize` (n and condition as above); `GaussNewton.SymmetricIsBitIdenticalOnEdgeEntries` does
it on a random J with ±0, subnormals, 1e150, -1e154 (products that overflow) and ±max / 2 in a third of its entries,
and again with one entry of 0.75 max. Reversing the order of i in the blocks makes both fail; removing the check of
max / 2 makes the second fail.

| Measure | A (+ L5) | B (+ L4) | B / A |
| --- | ---: | ---: | ---: |
| `BM_Kernel/gauss_newton/10` (old loop; two interleaved pairs) | 0.266 and 0.269 µs | 0.252 and 0.270 µs | noise |
| `BM_Kernel/gauss_newton/100` | 159 µs | 81.8 µs | 0.51 |
| `BM_Kernel/gauss_newton/300` | 4.63 ms | 2.09 ms | 0.45 |
| `BM_Kernel/gauss_newton/1000` | 178 ms | 59.3 ms | 0.33 |
| `BM_Solve/risk_parity/target_volatility/100` | 20.1 ms | 16.1 ms | 0.80 |
| `BM_Solve/risk_parity/box/100` | 34.3 ms | 30.0 ms | 0.88 |
| `risk_parity/target_volatility/1000` (one solve, `portopt_profile`) | 11.7 s | 6.37 s | 0.55 |
| `risk_parity/box/1000` | 26.0 s | 18.4 s | 0.71 |
| `risk_parity/target_return/1000` | 24.3 s | 16.3 s | 0.67 |

The cases of the quick suite that do not call `gauss_newton` moved by 0.94 to 1.11, and two unchanged kernels by 1.42
and 1.43 (`cholesky/100`, `qr/100`): noise of the desktop. All 70 cases at n ≤ 100 and the three `risk_parity` cases at
n = 1000 are bit for bit identical. In the assembly (`objective.cpp.obj`, everything inlined in `gauss_newton`) the
block loop is four `movupd` loads, the broadcasts and lane swaps of the two pairs (`unpcklpd`, `unpckhpd`,
`shufpd`), then 8 `mulpd` and 8 `addpd` into 8 registers that hold the 16 sums; they are stored once, after the loop.

### L7: J0 kept along the MAX_RETURN path

Every `qp::solve` built J0 = L^-T from its factor (`empty_active_set`, n³/6 multiply-subtracts after L1), and every
point of the MAX_RETURN path solves Q(t) with the same factor, that of `Data`: only the linear term and the rows
change. `qp::solve(problem, cache)` takes a `qp::BasisCache` (the factor it was built from and its J0): it starts
from a copy of the cached J0 when `problem.hessian_factor` has the same bits as the cached factor (one `memcmp` of n²
doubles), and builds J0 and replaces the cache otherwise. The `Path` of a solve holds one cache. The key is the bits
of the factor, so a problem with another factor (even one that compares equal, as -0 and +0 do) never starts from
the wrong J0, and J0 is the output of the same code on the same input whether it is cached or not. `empty_active_set`
is now `active_set_from(initial_basis(factor, n), n)`, the same operations. Two more allocations per path (the cached
factor and J0); each solve still copies J0 into its own basis.

`GoldfarbIdnani.BasisCacheIsBitIdenticalToAFreshBasis` solves, for n and condition as above, six points of a
long-only path (linear term `-t mu`, a row held as an equality and released again), each twice with the cache, and
between them a problem with the factor 2 L, comparing every field of the solution by `memcmp` with `qp::solve`
without the cache, and checking that the cache holds the factor of the last problem. Building J0 only once whatever
the factor makes it fail. `MaxReturn.PathBasisCacheIsBitIdenticalToAFreshBasis` runs `solve_at` on two paths of the
same spec (free, long-only, box, leverage, cash; n = 12) at t = 0 and 1e-4 doubling to 40, emptying the second path's
cache before each point: the points (solution, volatility, shift) are the same bits, and after each point the problem's
factor and the cached factor still have the bits of `Data`'s. On the benchmark instances a temporary count in
`qp::solve` (not committed) gave one rebuild per path and a reuse at every other point: 2, 6, 7, 6, 7 reuses for
`free`, `long_only`, `box`, `leverage`, `cash` at n = 10; 3, 9, 9, 21, 10 at n = 100; 3, 10, 10, 55, 10 at n = 1000.

| Measure | A (+ L4) | B (+ L7) | B / A |
| --- | ---: | ---: | ---: |
| `BM_Solve/max_return/free/100` (9 repetitions, two interleaved pairs) | 283 and 261 µs | 167 and 149 µs | 0.59 and 0.57 |
| `BM_Solve/max_return/long_only/100` | 3.04 and 3.16 ms | 2.79 and 3.04 ms | 0.92 and 0.96 |
| `BM_Solve/max_return/box/100` | 3.63 and 3.83 ms | 3.37 and 3.12 ms | 0.93 and 0.81 |
| `BM_Solve/max_return/leverage/100` | 11.4 and 9.80 ms | 10.5 and 8.54 ms | 0.92 and 0.87 |
| `BM_Solve/max_return/cash/100` | 3.32 and 3.45 ms | 2.82 and 2.91 ms | 0.85 and 0.84 |
| `max_return/free/1000` (one solve, `portopt_profile`) | 189 ms | 91.7 ms | 0.49 |
| `max_return/long_only/1000` | 2.41 s | 2.18 s | 0.90 |
| `max_return/box/1000` | 2.46 s | 1.98 s | 0.80 |
| `max_return/leverage/1000` | 15.9 s | 14.4 s | 0.91 |
| `max_return/cash/1000` | 2.26 s | 2.16 s | 0.96 |

The n = 10 cases moved by 0.84 to 1.10 in both directions (J0 costs 0.5 µs there), and the cases of the quick suite
that do not follow a path by 0.95 to 1.47: noise of the desktop. The saving at n = 1000, 2.4 s over the five cases, is
under the 3.7 s of the profile's estimate (88 reused bases × 45 ms): each point still copies the 8 MB of J0, and the
single solves are not repeated. All 70 cases at n ≤ 100 and the 15 `max_return` cases of the three sizes are bit for
bit identical. In the assembly (`goldfarb_idnani.cpp.obj`) the check is one `call memcmp` and the copy a `memcpy`;
both overloads of `qp::solve` call the same out-of-line `solve_factored`.
