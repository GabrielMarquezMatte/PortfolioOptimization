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
- Anything but one machine, one compiler and generic code: see [What was left out](#what-was-left-out-and-why).

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

### L6a: Dekker's product in the compensated dot

`compensated_dot` took each product's error from `std::fma(a, b, -p)`, which a generic x86-64 build compiles into a
call of the `fma` of `ucrtbase.dll` per term. `two_product_error` now computes it by Dekker's product (Dekker 1971):
Veltkamp's split of each factor into two halves of 26 bits (`c = (2^27 + 1) a`, `a_hi = c - (c - a)`,
`a_lo = a - a_hi`) and `(((a_hi b_hi - p) + a_hi b_lo) + a_lo b_hi) + a_lo b_lo`, every product and every partial sum
exact. It is exact, so it equals what `std::fma` returns, the exact error, when the split cannot overflow and the error
is not below the smallest subnormal (its last bit is ulp(a) ulp(b)): both factors of magnitude in [2^-485, 2^500], or
zero (an error of zero is +0 both ways, since no partial sum is -0). Outside that range (subnormals, |x| < 2^-485,
|x| > 2^500, inf, NaN) it calls `std::fma` as before; the check is per term. A scratch probe of the raw Dekker
product against `std::fma` over random exponents found its first differences at an exponent sum of -998 and at
|a| ≥ 2^997 (where the split overflows), outside the range used. `-ffp-contract=off` keeps every operation rounded as
written.

`CompensatedDot.DekkerMatchesFmaOnEdgeCases` compares `two_product_error(a, b, a b)` with `std::fma(a, b, -(a b))`
by `memcmp`: every pair (both signs) of 26 edge values (±0, the ends of the range and the doubles next to them, 2^996,
2^997, max, 2^-1022, 2^-1000, subnormals, inf, 2^27 + 1 and a few mantissas), then 2000 random pairs near each end of
the range and 2000 over every exponent. Using Dekker's product at every finite magnitude makes it fail (815
mismatches). `CompensatedDot.IsBitIdenticalToTheFmaDot` compares `compensated_dot` with the old loop, kept verbatim, on
covariance rows (n and condition as above) and on vectors of magnitudes 2^-600 to 2^600.

Generic build (`gcc-bench`), three interleaved A/B pairs, 5 repetitions:

| Measure | A (`std::fma`) | B (Dekker) | B / A |
| --- | ---: | ---: | ---: |
| `BM_Kernel/compensated_dot/10` | 0.042 and 0.043 µs | 0.031 and 0.032 µs | 0.74 |
| `BM_Kernel/compensated_dot/100` | 0.403 to 0.443 µs | 0.302 to 0.339 µs | 0.75 |
| `BM_Kernel/compensated_dot/1000` | 3.92 to 4.41 µs | 2.91 to 3.59 µs | 0.74 to 0.81 |
| `BM_Kernel/covariance_times/10` | 0.439 and 0.445 µs | 0.337 and 0.339 µs | 0.77 |
| `BM_Kernel/covariance_times/100` | 40.2 to 41.3 µs | 29.6 to 32.9 µs | 0.74 to 0.80 |
| `BM_Kernel/covariance_times/1000` | 3.92 to 4.39 ms | 2.91 to 2.99 ms | 0.68 to 0.76 |

3.9 ns a term with the call against 2.9 ns, under the exploratory 1.9 times of the research (2.20 ns), whose Dekker
loop had no range check. A version that checked both vectors once before a branch-free Dekker loop was slower than
the per-term check (3.61 µs at n = 1000). On the solves, 9 repetitions, two interleaved pairs:

| Case | B / A | Case | B / A |
| --- | ---: | --- | ---: |
| `risk_parity/target_volatility/10` | 0.91 and 0.93 | `risk_parity/target_volatility/100` | 0.89 and 0.90 |
| `risk_parity/target_return/10` | 0.94 and 0.93 | `risk_parity/target_return/100` | 0.94 and 0.95 |
| `risk_parity/box/10` | 0.94 and 0.96 | `risk_parity/box/100` | 0.96 and 0.96 |
| `risk_parity/budget/10` (no compensated dot) | 0.99 and 1.00 | `risk_parity/budget/100` | 0.98 and 1.00 |

At n = 1000, where `covariance_times` is about 2% of `risk_parity`, one solve each: `target_volatility` 6.79 and
6.66 s, `box` 18.9 and 18.6 s, `target_return` 16.4 and 16.6 s, inside the spread of single runs. All 70 cases at
n ≤ 100 and these three are bit for bit identical. In the assembly (`objective.cpp.obj`, inlined in
`compensated_dot` and `covariance_times`) the term is the range check (`andpd`, `comisd` and a branch for each
factor, taken to the `call fma` only out of range), the two splits (`mulsd`, `subsd`, `subsd`, `subsd`), the four
products and the three additions in the order written, then the TwoSum.

The same A/B in a native build (`native-release` configured with GCC, Ninja and `-march=native`, which turns on
`-mfma` and `-mavx2` on this Zen 3; `std::fma` is then one `vfmsub132sd`), three interleaved pairs:

| Measure | A (`std::fma`, inline) | B (Dekker) | B / A |
| --- | ---: | ---: | ---: |
| `BM_Kernel/compensated_dot/10` | 0.011 to 0.015 µs | 0.029 to 0.037 µs | 2.5 |
| `BM_Kernel/compensated_dot/100` | 0.114 to 0.130 µs | 0.288 to 0.292 µs | 2.3 |
| `BM_Kernel/compensated_dot/1000` | 1.07 to 1.20 µs | 2.77 to 3.25 µs | 2.6 |
| `BM_Kernel/covariance_times/100` | 11.8 to 12.4 µs | 28.9 to 36.9 µs | 2.6 |
| `BM_Kernel/covariance_times/1000` | 1.16 to 1.44 ms | 2.94 to 4.86 ms | 2.5 |

With hardware FMA the call disappears and the fma is 1.1 ns a term, so Dekker's product (2.8 ns) is 2.5 times
slower there. The shipped binaries and the official numbers are generic, where it is 1.35 times faster; a native
build that wants the inline fma back needs the dispatch of L6b (or a compile-time choice of `std::fma` under
`__FMA__`), not done here.

### L8: allocation hygiene, where the counter shows it

The counter ([Measuring](measuring.md#counters)) said where allocations are before anything was changed: 29 054 per
pass of the 70 cases at n ≤ 100 (6369 at n = 10, 22 685 at n = 100). At about 35 ns a `malloc`/`free` pair on this
machine (a scratch loop, 29 to 32 ns for blocks of 80 B to 8 kB), 6369 allocations are about 0.22 ms of a 0.78 ms
pass at n = 10: allocation was a quarter of the time of the small cases, and nothing at n ≥ 100 (about 1% of the
cases there), except in the gross exposure cases, where the generated rows were rebuilt on every iteration. Four
changes, all of data movement, so every solution is bit for bit the same:

- **Bound rows in place** (`problems/constraints.cpp`, `problems/homogeneous.cpp`): each bound was a temporary vector
  of n doubles inserted into a row vector that grew by doubling, and a push into `bounds`: n allocations plus the growth
  per solve. The rows are now written at the end of rows reserved for the number of finite bounds.
- **Family buffers** (`qp/goldfarb_idnani.cpp`): the sign vector and the normal of the most violated row of the gross
  exposure family are two vectors kept in `Constraints` and refilled, and `discard_inactive_generated`, which moved the
  generated rows out and copied the active ones back into three fresh vectors on every loop, copies them into a kept
  `GeneratedRows` and swaps.
- **Reserves** (`qp/active_set.hpp`, `qp/goldfarb_idnani.cpp`): the active list (n entries), the gathered rows and
  right-hand sides (one allocation instead of a copy and a growth) and the list of active inequalities.

| Case | n = 10 before, after | n = 100 before, after | n = 1000 before, after |
| --- | ---: | ---: | ---: |
| `min_variance/long_only` | 53, 26 | 159, 26 | 1073, 26 |
| `min_variance/box` | 67, 26 | 262, 26 | 2076, 26 |
| `max_sharpe/box` | 69, 28 | 264, 28 | 2078, 28 |
| `max_return/long_only` | 276, 203 | 471, 216 | 1466, 233 |
| `max_return/leverage` | 722, 402 | 9010, 1493 | 183 750, 4919 |
| `utility/leverage` | 53, 42 | 1017, 87 | 17 569, 114 |
| `tracking_error/leverage` | 31, 29 | 600, 83 | 11 708, 120 |
| `risk_parity/target_return` | 720, 626 | 2046, 1386 | 2922 after |
| `risk_parity/box` | 959, 786 | 3208, 2006 | 2678 after |
| `risk_parity/target_volatility` | 1606, 1606 | 1700, 1700 | 1669 after |
| all 70 cases at n ≤ 100 | 6369, 4998 | 22 685, 8350 | - |

Same session, A (before) and B (after) two copies of the executables, quick suite twice interleaved, and a filtered
A/B of 4 rounds of 5 repetitions at n = 10:

| Measure | B / A |
| --- | ---: |
| `min_variance/long_only/10`, `box/10`, `cash/10` | 0.72, 0.67, 0.66 |
| `utility/long_only/10`, `box/10`, `cash/10` | 0.74, 0.68, 0.72 |
| `max_sharpe/long_only/10`, `box/10`, `cash/10` (after `homogeneous.cpp`) | 0.73, 0.71, 0.74 |
| `max_diversification/long_only/10`, `box/10` | 0.74, 0.69 |
| `tracking_error/long_only/10`, `box/10` | 0.83, 0.73 |
| `max_return/long_only/10`, `box/10`, `leverage/10` | 0.60 to 0.87, 0.83, 0.77 |
| `risk_parity/{box, target_volatility}/10` (only the `qp::solve` reserves) | 0.88 to 0.93, 0.91 to 0.94 |
| n = 100 cases with bounds (quick suite, two pairs) | 0.87 to 1.16, within the noise of the desktop |
| `utility/leverage/100`, `max_return/leverage/100` (4 rounds) | 0.85, 0.94 |
| `utility/leverage/1000` (one solve, two pairs) | 2.40 and 2.21 s, 1.06 and 1.22 s |
| `tracking_error/leverage/1000` | 2.05 and 2.28 s, 1.36 and 1.45 s |
| `max_return/leverage/1000` | 13.7 and 13.8 s, 12.7 and 11.6 s |

All 70 cases at n ≤ 100 and 12 cases at n = 1000 with bounds, cash or leverage (`min_variance/{box, cash,
target_return}`, `max_return/{leverage, cash}`, `utility/leverage`, `tracking_error/leverage`,
`max_sharpe/{box, cash, long_only}`, `max_diversification/box`, `risk_parity/target_volatility`) are bit for bit identical (`compare-solutions
--require-identical`); the fuzz batteries of the gross exposure family and `ctest --preset gcc-release` pass.

Not done, with the number that decides it: `risk_parity` still allocates 1600 times per solve at n = 10 (1627 before the
`qp::solve` reserves) in about 40 call sites of `riskparity/` and `problems/` that return vectors (`linearize`,
`moments_of`, `value_change`, `step_at`, `record`, the KKT and the polish), which is at most 25% of the 0.22 ms of
`target_volatility/10` and 0.4% of its n = 100 time. Removing them needs a workspace passed through the descent, a
change of its interfaces for a gain only at n ≈ 10 that the quick suite cannot resolve above its noise at that size
(14% of standard deviation in the baseline): a proposal in [Plan 6 profile 2](profile-2.md#proposals-for-the-author).

### Exception safety of the L7 cache

`qp::solve(problem, cache)` assigned `cache.factor` before `initial_basis` ran: a `bad_alloc` in the basis left the new
key over the old basis, and the next solve with that factor would have started from the wrong J0. Both are now built
into locals and moved in. `GoldfarbIdnani.BasisCacheKeepsItsPairWhenAnAllocationFails` fails the n-th allocation of a
rebuild for every n (the test executable replaces `operator new` with a countdown that is off outside that test) and
checks that the cache is the old pair or the new one, and that the next solve equals a fresh one; it fails on the old
code at the first allocation.

### No accidental FMA (guard 3)

On the native build above (`-march=native`: `-mfma` enabled, `g++ -Q --help=target`), `objdump -d` of the 20
objects of `portopt_core` has three `vfm*` instructions, each a `vfmsub132sd`: one in `two_product_error`, the only
function that calls `std::fma`, and one in each of `compensated_dot` and `covariance_times`, where it is inlined. No
other function has one: `-ffp-contract=off` holds on every product of the core. `ctest` on that build passed (337
tests, the benchmark self-tests included).

The C API. `portopt_capi` was compiled without `-ffp-contract=off` until the end of plan 6: it does no arithmetic of
its own, but it includes core headers, and with LTO (or an inline function) core code compiled in its units would be
contracted. The flag is now on it (`build: -ffp-contract=off on the C API`). Before the change the 7 objects of
`portopt_capi` and `libportopt.dll` of the native build had the same three `vfm*`, in `two_product_error`,
`compensated_dot` and `covariance_times` of the core, and the C interface gave the same bits as the C++ path on 70 of
70 cases at n ≤ 100 and on the 18 `min_variance` cases of the three sizes; after it, the same three and the same bits
(`portopt_profile --via capi` against `--via cpp`, `compare-solutions --require-identical`). The guard for the library:

```
objdump -d -C -M intel --no-show-raw-insn build/native-release/libportopt.dll | grep -cE 'vf(n?m)(add|sub)'
```

must print 3, and the three must be those functions (`objdump` with the function headers shows which). On a build with
`-flto` and the flag only on the core it was 25, in `po_solve` and the metrics ([profile 2](profile-2.md#l8-l9-and-l10)).

## Plan 6 final: before and after

The baseline against the end of plan 6, both measured in one session on 2026-11-17, one right after the other on the
machine of [Machine and run](#machine-and-run). The comparison with the table of [Plan 6 baseline](#plan-6-baseline)
(2026-11-16) is not made: runs of different days move by more than their spread
([Measuring](measuring.md#ab)).

| | Before (A) | After (B) |
| --- | --- | --- |
| Code | `1297db0`, the last commit before any change of `src/` (its `src/` is that of `4cc0b82`, the baseline; the commits between them are documents), built in a separate worktree | `6419150`, HEAD; its `src/` is that of the last lever (`66ca389`), the commits after it change build files, tests, benchmarks and documents only |
| Build | `gcc-bench`, generic x86-64, `-O3 -DNDEBUG -ffp-contract=off` | the same |
| Harness | `portopt_benchmarks` and `portopt_profile` of `1297db0` | those of HEAD; the cases, instances and counters are the same (the diff of `tests/benchmarks` between the two is the C API benchmark and `--via`) |

- n = 10 and 100: 5 repetitions, the median of the wall time of `BM_Solve`; the 7 objectives were run as separate
  slices, A and B alternating and the order swapped on every slice (about 4 minutes of wall time for each side).
- n = 1000: one repetition, one timed solve of each of the 35 cases (`portopt_profile --max-n 1000 --repeat 1 --filter
  /1000`, B and then A), so no spread is measured; the same cases run four at a time were up to 1.4 times slower than
  alone ([Machine and run](#machine-and-run)). The sum of the 35 solves is 290.1 s before and 78.8 s after (3.7
  times).
- Spread at n = 10 and 100: the median relative standard deviation over the 140 runs of 5 repetitions is 4.7%, 65
  of the 140 are above 5% and the largest is 33%. Ratios between 0.9 and 1.1 are not results.

| Case | n = 10 before | after | ratio | n = 100 before | after | ratio | n = 1000 before | after | ratio |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `min_variance/box` | 6.15 µs | 3.69 µs | 0.60 | 936.57 µs | 231.29 µs | 0.25 | 1.03 s | 193.01 ms | 0.19 |
| `min_variance/cash` | 3.59 µs | 1.83 µs | 0.51 | 502.72 µs | 74.54 µs | 0.15 | 563.35 ms | 62.12 ms | 0.11 |
| `min_variance/free` | 2.53 µs | 2.00 µs | 0.79 | 458.82 µs | 79.67 µs | 0.17 | 386.06 ms | 58.35 ms | 0.15 |
| `min_variance/leverage` | 4.38 µs | 3.68 µs | 0.84 | 445.62 µs | 70.30 µs | 0.16 | 424.71 ms | 57.45 ms | 0.14 |
| `min_variance/long_only` | 4.90 µs | 3.34 µs | 0.68 | 797.62 µs | 231.13 µs | 0.29 | 854.46 ms | 237.21 ms | 0.28 |
| `min_variance/target_return` | 5.02 µs | 3.21 µs | 0.64 | 939.43 µs | 280.05 µs | 0.30 | 848.98 ms | 191.06 ms | 0.23 |
| `max_return/box` | 54.74 µs | 42.31 µs | 0.77 | 11.35 ms | 2.97 ms | 0.26 | 16.23 s | 2.45 s | 0.15 |
| `max_return/cash` | 47.03 µs | 32.46 µs | 0.69 | 13.51 ms | 2.80 ms | 0.21 | 10.85 s | 2.17 s | 0.20 |
| `max_return/free` | 7.87 µs | 5.73 µs | 0.73 | 1.80 ms | 149.33 µs | 0.08 | 1.57 s | 72.96 ms | 0.05 |
| `max_return/leverage` | 57.12 µs | 43.92 µs | 0.77 | 18.79 ms | 7.65 ms | 0.41 | 32.88 s | 13.22 s | 0.40 |
| `max_return/long_only` | 38.42 µs | 29.98 µs | 0.78 | 9.79 ms | 2.66 ms | 0.27 | 12.10 s | 2.42 s | 0.20 |
| `max_sharpe/box` | 6.60 µs | 4.16 µs | 0.63 | 1.06 ms | 447.33 µs | 0.42 | 1.04 s | 488.53 ms | 0.47 |
| `max_sharpe/cash` | 6.48 µs | 3.88 µs | 0.60 | 872.50 µs | 266.55 µs | 0.31 | 904.69 ms | 209.91 ms | 0.23 |
| `max_sharpe/free` | 5.89 µs | 2.98 µs | 0.51 | 545.38 µs | 109.87 µs | 0.20 | 351.75 ms | 60.94 ms | 0.17 |
| `max_sharpe/leverage` | 5.68 µs | 4.18 µs | 0.74 | 447.83 µs | 78.22 µs | 0.17 | 378.47 ms | 59.69 ms | 0.16 |
| `max_sharpe/long_only` | 6.79 µs | 3.77 µs | 0.56 | 913.00 µs | 289.73 µs | 0.32 | 795.83 ms | 183.63 ms | 0.23 |
| `utility/box` | 6.07 µs | 3.77 µs | 0.62 | 1.61 ms | 430.65 µs | 0.27 | 2.44 s | 427.86 ms | 0.18 |
| `utility/cash` | 4.97 µs | 2.88 µs | 0.58 | 1.25 ms | 413.19 µs | 0.33 | 1.66 s | 423.95 ms | 0.26 |
| `utility/free` | 2.63 µs | 2.04 µs | 0.77 | 463.96 µs | 81.55 µs | 0.18 | 366.72 ms | 56.15 ms | 0.15 |
| `utility/leverage` | 5.44 µs | 3.73 µs | 0.69 | 1.43 ms | 859.48 µs | 0.60 | 1.96 s | 1.37 s | 0.70 |
| `utility/long_only` | 5.80 µs | 3.56 µs | 0.61 | 1.20 ms | 420.52 µs | 0.35 | 1.88 s | 401.00 ms | 0.21 |
| `tracking_error/box` | 6.73 µs | 4.58 µs | 0.68 | 2.61 ms | 931.32 µs | 0.36 | 3.60 s | 1.16 s | 0.32 |
| `tracking_error/cash` | 5.46 µs | 4.08 µs | 0.75 | 2.25 ms | 974.81 µs | 0.43 | 3.06 s | 1.14 s | 0.37 |
| `tracking_error/free` | 4.81 µs | 3.91 µs | 0.81 | 1.62 ms | 644.82 µs | 0.40 | 1.97 s | 742.97 ms | 0.38 |
| `tracking_error/leverage` | 4.71 µs | 4.06 µs | 0.86 | 2.03 ms | 1.07 ms | 0.53 | 2.77 s | 1.74 s | 0.63 |
| `tracking_error/long_only` | 5.90 µs | 4.47 µs | 0.76 | 2.19 ms | 1.00 ms | 0.46 | 3.40 s | 1.17 s | 0.35 |
| `max_diversification/box` | 7.93 µs | 5.11 µs | 0.64 | 977.13 µs | 435.46 µs | 0.45 | 973.72 ms | 391.46 ms | 0.40 |
| `max_diversification/cash` | 6.07 µs | 4.17 µs | 0.69 | 832.81 µs | 263.07 µs | 0.32 | 792.86 ms | 183.43 ms | 0.23 |
| `max_diversification/free` | 4.35 µs | 3.34 µs | 0.77 | 473.99 µs | 115.93 µs | 0.24 | 357.11 ms | 59.65 ms | 0.17 |
| `max_diversification/leverage` | 5.76 µs | 4.91 µs | 0.85 | 475.37 µs | 80.52 µs | 0.17 | 372.33 ms | 49.80 ms | 0.13 |
| `max_diversification/long_only` | 5.58 µs | 3.94 µs | 0.71 | 850.79 µs | 241.64 µs | 0.28 | 810.90 ms | 197.79 ms | 0.24 |
| `risk_parity/box` | 165.66 µs | 132.87 µs | 0.80 | 72.71 ms | 28.86 ms | 0.40 | 86.81 s | 21.52 s | 0.25 |
| `risk_parity/budget` | 5.55 µs | 6.19 µs | 1.12 | 1.09 ms | 630.90 µs | 0.58 | 1.60 s | 490.70 ms | 0.31 |
| `risk_parity/target_return` | 99.21 µs | 91.32 µs | 0.92 | 36.24 ms | 15.69 ms | 0.43 | 68.25 s | 18.11 s | 0.27 |
| `risk_parity/target_volatility` | 259.70 µs | 228.72 µs | 0.88 | 35.87 ms | 14.95 ms | 0.42 | 25.82 s | 6.98 s | 0.27 |

Geometric mean of the ratio B / A by objective over its cases (the range in parentheses):

| Objective | n = 10 | n = 100 | n = 1000 |
| --- | ---: | ---: | ---: |
| `min_variance` | 0.67 (0.51 to 0.84) | 0.21 (0.15 to 0.30) | 0.17 (0.11 to 0.28) |
| `max_return` | 0.75 (0.69 to 0.78) | 0.22 (0.08 to 0.41) | 0.16 (0.05 to 0.40) |
| `max_sharpe` | 0.60 (0.51 to 0.74) | 0.27 (0.17 to 0.42) | 0.23 (0.16 to 0.47) |
| `utility` | 0.65 (0.58 to 0.77) | 0.32 (0.18 to 0.60) | 0.25 (0.15 to 0.70) |
| `tracking_error` | 0.77 (0.68 to 0.86) | 0.43 (0.36 to 0.53) | 0.40 (0.32 to 0.63) |
| `max_diversification` | 0.73 (0.64 to 0.85) | 0.28 (0.17 to 0.45) | 0.22 (0.13 to 0.40) |
| `risk_parity` | 0.92 (0.80 to 1.12) | 0.45 (0.40 to 0.58) | 0.27 (0.25 to 0.31) |
| all cases | 0.71 (0.51 to 1.12) | 0.29 (0.08 to 0.60) | 0.23 (0.05 to 0.70) |

The slowest case of each objective, before and after:

| Objective | n = 10 | n = 100 | n = 1000 |
| --- | ---: | ---: | ---: |
| `min_variance` | 6.15 µs to 3.69 µs (`box`) | 939 µs to 280 µs (`target_return`) | 1.03 s (`box`) to 237 ms (`long_only`) |
| `max_return` | 57.1 µs to 43.9 µs (`leverage`) | 18.8 ms to 7.65 ms (`leverage`) | 32.9 s to 13.2 s (`leverage`) |
| `max_sharpe` | 6.79 µs (`long_only`) to 4.18 µs (`leverage`) | 1.06 ms to 447 µs (`box`) | 1.04 s to 489 ms (`box`) |
| `utility` | 6.07 µs to 3.77 µs (`box`) | 1.61 ms (`box`) to 859 µs (`leverage`) | 2.44 s (`box`) to 1.37 s (`leverage`) |
| `tracking_error` | 6.73 µs to 4.58 µs (`box`) | 2.61 ms (`box`) to 1.07 ms (`leverage`) | 3.60 s (`box`) to 1.74 s (`leverage`) |
| `max_diversification` | 7.93 µs to 5.11 µs (`box`) | 977 µs to 435 µs (`box`) | 974 ms to 391 ms (`box`) |
| `risk_parity` | 260 µs to 229 µs (`target_volatility`) | 72.7 ms to 28.9 ms (`box`) | 86.8 s to 21.5 s (`box`) |

What the table says:

- The gain grows with n: all cases 0.71 at n = 10, 0.29 at n = 100, 0.23 at n = 1000. The levers are on the n³ kernels
  and on the iteration of the active set, which are the time at n ≥ 100; at n = 10 what is left is call overhead and
  allocation. The best ratios are the cases that were mostly the fixed cost of the base: `max_return/free` 0.05 at
  n = 1000 (J0 cached along the path), `min_variance/cash` 0.11, `min_variance/{free, leverage}` 0.14 to 0.15.
- The smallest gains are the cases whose time is not in the kernels the levers touched: `utility/leverage` and
  `tracking_error/leverage` (0.70 and 0.63 at n = 1000: the gross exposure family, whose dense generated rows are the
  cost, [profile 2](profile-2.md#hot-functions)), `max_return/leverage` (0.40), `max_sharpe/box` (0.47) and
  `tracking_error` in general (0.32 to 0.63: the QR of the series, 0.365 s of its 0.98 s at n = 1000).
- `risk_parity` ends at 0.25 to 0.31 at n = 1000 and 0.40 to 0.58 at n = 100; at n = 10 it is 0.80 to 0.92, the
  allocator and `compensated_dot` of its small solves.
- **One case is slower: `risk_parity/budget/10`, 5.55 to 6.19 µs in the table.** Three more pairs of 9 repetitions
  gave 5.64, 5.47 and 5.41 µs before and 5.93, 5.90 and 5.88 after (1.05 to 1.09): a real loss of about 0.4 µs, the
  one case of the 105 that the levers did not help. It has no QP and no compensated dot; it factors by Cholesky (L2:
  the right-looking form with the `memcpy` start is 34 ns slower than the left-looking one at n = 10,
  [L2](#l2-right-looking-cholesky-two-columns-at-a-time)) and runs the Spinu Newton steps. Not isolated further: 7% of
  6 µs is under what a caller would see, and the same case is 0.58 at n = 100 and 0.31 at n = 1000.

### What was left out, and why

Everything below has its number in [profile 2](profile-2.md) or in a lever above; none is implemented.

| Left out | Why |
| --- | --- |
| L9, a `Moments` cache in the Armijo halvings | At most 2.8% of the samples of `risk_parity` at n ≤ 100, under what the quick suite resolves (3 to 5% between two runs of one binary); not visible at n = 1000 |
| L10, `-flto` on the core | 0.87 to 1.00 at n = 10 (about 0.93 geometric mean), nothing at n ≥ 100, a serial link of 25 s; it needs `-ffp-contract=off` on every target that includes inline core code (done for `portopt_capi`, below) and a guard. A small separate task for callers of many small problems |
| ISA multiversioning (AVX2 and FMA chosen at run time) | 0.75 to 0.87 at n = 100 with no numerical risk (70 of 70 bit-identical to the generic build), but it needs a dispatch by hand (MinGW has no `ifunc`; MSVC another mechanism) and keeps the generic binary the shipped one; the native run below repeats the gain |
| FMA-aware Dekker (the inline `std::fma` under `__FMA__`) | In a native build Dekker's product is 2.5 times slower than the inline `vfmsub132sd`; the shipped binaries are generic, where Dekker is 1.35 times faster than the `call fma` of `ucrtbase`. It needs the dispatch above |
| Class-B kernels (reassociated reductions: four or eight accumulators) | They change the rounding; 2.8 to 6.4 times on a scratch dot product, about 0.5 ms of the 1.2 ms of an iteration of `max_return/leverage` at n = 1000 by arithmetic, not by experiment. They need the A/B of the weights inside the bands of `expect_user_kkt` and `expect_risk_parity_kkt`, which this plan did not loosen |
| A warm start of the QP of `risk_parity` | It changes the path of the QP and `iterations`; the iterations are not wasted (one per active row), so it would have to build the active set at once; the ceiling is 43 to 60% of the slow `risk_parity` cases if the iterations were free |
| A workspace for the descent of `risk_parity` | 25% of 0.22 ms at n = 10, 0.4% at n = 100, about 40 call sites; only for callers of many small `RISK_PARITY` problems |
| The counters `qp_solves`, `descent_iterations`, `polish_steps` | They need a counter in `src/` ([Measuring](measuring.md#counters)); counted once with a debugger instead |
| More than one repetition at n = 1000 | Five repetitions of the 35 cases took close to an hour before the levers and would take about 15 minutes now; the single solve is a size of effect, not an interval |
| Isolated, pinned or fixed-frequency runs; other machines, Linux, macOS, Clang | Not available here; the numbers are of one desktop |
| A comparison with another library | None is measured, so none is made |

### The builds other than `gcc-bench`

**`msvc-release`** (Visual Studio 18 2026, MSVC 19.51, `cmake --preset msvc-release`; the build takes 1 min 40 s): the
test battery did not build at `1297db0` and did not at HEAD before the fix, for reasons older than plan 6: the `getenv`
and `fopen` of the two measurement tests (C4996) and an `int64_t` range converted to `vector<double>` in
`max_sharpe_test.cpp` (C4244), both errors under `/W4 /WX`; and, once it builds, five self-tests of the reference
(`RiskParityReference.Kkt*` and `HessianByDifferencesMatchesTheObjective`) that need an extended `long double` and
failed where `long double` is `double`, the same five at `1297db0`. Fixed in the test code only (`3570570`):
`_CRT_SECURE_NO_WARNINGS` on `portopt_tests` for MSVC, a loop instead of the range constructor, and the
`kLongDoubleIsExtended` skip that the sibling tests already use. Nothing in `src/` and nothing of the plan 6 tests
failed on MSVC: the replaced `operator new` of `goldfarb_idnani_test.cpp` builds and
`BasisCacheKeepsItsPairWhenAnAllocationFails` passes. Result: **`ctest --preset msvc-release`, 100% of 336 tests
passed**, 16 of them skipped (the 11 that already skipped for the `long double` and the 5 above), 16 disabled (the
measurement tests). The speed of the MSVC build was not measured beyond the samples of [profile 2](profile-2.md)
(1.3 to 1.7 times slower than `gcc-bench` on the slow cases).

**`native-release`.** The preset alone selects no compiler and sets no `-march`, so it is not native. It was configured as
`cmake --preset native-release -G Ninja -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_FLAGS=-march=native`
(`znver3`; `-mfma` and `-mavx2` enabled, `__FMA__` defined). `ctest` there: 338 tests passed (the 336 and the two
benchmark self-tests), and its 70 solutions at n ≤ 100 are bit for bit identical to those of `gcc-bench`
(`compare-solutions --require-identical`). The difference from the generic build, quick suite run twice interleaved
(the shipped binaries and every official number are generic; this only records the gap):

| Measure | native / generic (two runs) |
| --- | ---: |
| solves at n = 100 | 0.71 to 0.94: `max_return`, `max_sharpe`, `utility`, `max_diversification` and `risk_parity/budget` 0.76 to 0.81, `tracking_error` 0.70 to 0.75, `risk_parity/target_volatility` 0.89 to 0.91, `risk_parity/box` 0.86 to 0.94 |
| solves at n = 10 | 0.90 to 1.02 (`risk_parity/budget` 1.01 and 1.02, `max_return/long_only` 0.92 and 1.00) |
| `cholesky/100`, `qr/100`, `qp_iteration/100`, `empty_active_set/100` | 0.71 to 0.79, 0.73 to 0.90, 0.80 to 0.82, 0.86 to 0.94 |
| `gauss_newton/100`, `compensated_dot/100`, `covariance_times/100` | 1.03 to 1.04, 0.94 to 1.03, 0.96 |

That is [proposal 1 of profile 2](profile-2.md#proposals-for-the-author) again: AVX2 and FMA buy 6 to 30% at n = 100 and
nothing at n = 10; `gauss_newton`, already register-blocked, and the compensated dot, which Dekker's product keeps
off the hardware fma, do not gain.

### The C interface and Python

The cost of the layers above `po::solve` is in [Python and C layers against the solve](python-vs-arch.md): under 0.3 µs
per cycle in the C layer; in Python 30 µs per `Spec`, `Data` and `solve` (17.6 µs of it in `solve()`), which is 640% of
a 4.7 µs solve at n = 10, 16% of a 193 µs one at n = 100 and 0.02% of a 166 ms one at n = 1000.

### Verification at the end of plan 6

- `compare-solutions --require-identical`, `1297db0` against HEAD, both `gcc-bench`: **70 of 70 bit-identical at
  n ≤ 100 and 35 of 35 at n = 1000**, which is every case of the suite (status, objective, weights, multipliers and
  residual as hex-float), worst relative differences 0, iteration ratio 1.0000. Every lever of plan 6 is of risk 0.
- The C interface against the C++ path, `portopt_profile --via capi` against `--via cpp`: 70 of 70 identical at n ≤ 100
  in the generic build and in the native one, before and after `-ffp-contract=off` on `portopt_capi`, and the 18
  `min_variance` cases of the three sizes in the native one.
- `ctest --preset gcc-release`: 336 tests passed; `pytest`: 59 passed; `mypy --strict` of `python/src`,
  `python/tests` and `python/benchmarks`: clean; `python tests/style/check_style.py .`: clean.
