# MAX_RETURN tangent crossing

What the exact tangent crossing of `MAX_RETURN` (the face solve and its one-point KKT certificate in
`src/problems/parametric_path.cpp`) costs against the crossing it replaced. That crossing interpolated two QP
solutions and confirmed the root with one more QP solve. Tools and commands: [Measuring](measuring.md).

## Machine and run

AMD Ryzen 7 5700X, Windows 10, g++ 16.2.0 (MSYS2 UCRT64), `gcc-bench`, on 2026-11-19, in a desktop session that is not
isolated. No compiler or test process was running during the runs. A is commit `7287d84` (the interpolated crossing)
and B is the tangent crossing. Both binaries were kept side by side and run in the same session, each twice, with the
first run of each discarded:

```
portopt_benchmarks --benchmark_filter='^BM_Solve/max_return/.*/(10|100)$' --benchmark_min_time=0.3s \
    --benchmark_repetitions=3 --benchmark_report_aggregates_only=true --benchmark_out=mr_a.json --benchmark_out_format=json
python tests/benchmarks/quick.py compare mr_a.json mr_b.json
```

The flags are those of `quick.py run`. The filter takes every `max_return` variant, where the quick suite has only
`long_only`.

## Medians (microseconds)

| Case | A | B | B/A |
| --- | ---: | ---: | ---: |
| `max_return/free/10` | 5.67 | 7.49 | 1.321 |
| `max_return/free/100` | 187.4 | 356.8 | 1.904 |
| `max_return/long_only/10` | 30.12 | 29.02 | 0.964 |
| `max_return/long_only/100` | 2755 | 2629 | 0.954 |
| `max_return/box/10` | 42.81 | 39.68 | 0.927 |
| `max_return/box/100` | 2888 | 2797 | 0.969 |
| `max_return/leverage/10` | 42.57 | 40.94 | 0.962 (noise) |
| `max_return/leverage/100` | 7786 | 4010 | 0.515 |
| `max_return/cash/10` | 32.08 | 30.94 | 0.965 |
| `max_return/cash/100` | 2790 | 2747 | 0.985 (noise) |

The quick suite (`quick.py run`, the same pairing) has `max_return/long_only/10` at 0.985 and `/100` at 0.946. Every
other row of the suite is noise or within 0.956 to 1.042.

## Why `free` is slower

The face solve costs O(n p m + n r^2), with:

- p the free weights of the face,
- m its remaining rows,
- r = p - m the dimension of the face.

The weights held by an active bound row are fixed exactly and leave the solve, so `long_only`, `box` and `cash` reach
the cap on a face of small r. `free` has the budget as its only row: p = n, m = 1, r = n - 1. The QR of `L'Z`
(n x (n - 1)) is about 1.3e6 flops at n = 100.

The step it replaces is cheap in that case. It was one more Goldfarb-Idnani solve of `Q(t*)`, of 1 to 3 iterations
from the cached basis `L^-T`. In the bounded variants that solve took tens of iterations, and the tangent also saves
bisection steps (QP iterations 358 to 320 at `long_only/100`, 646 to 310 at `leverage/100`).

A first version held the active bound rows as columns of the null-space QR (p = n, m up to 90 at `long_only/100`). It
measured 1.107 and 1.109 on `max_return/long_only` in the quick suite. Fixing the bound weights is what brought those
cases below 1.

## Verdict

The `free` regression (1.32x at n = 10, 1.90x at n = 100) is accepted by the author. The other `max_return` variants
are faster or noise.

## The range space on a face with every weight free

The table above stopped at n = 100. At n = 1000 `free` was 2.9 times slower than before the tangent crossing (66 ms
to 192 ms): the QR of `L'Z` is n x (n - 1), O(n^3). The face solve now has a second branch, in the range space of
`Sigma = L L'`: with `G = L^-1 N = Q R`, `y = L'w = Q [R'^-1 b; (Q' t L^-1 v)_2]`. It needs no factorization beyond the
`L` of `Data`, so it applies only when no weight is fixed, and it is taken there when its `m + 4` triangular solves,
`(m + 4) n^2` flops, cost less than the null space's `4 n^2 m + 2 n (n - m)^2` (the rule in `range_space_cheaper`).
A crossing it does not certify is computed again by the null space. A face with fixed weights stays in the null space:
its range space would need the factor of `Sigma_FF`, `n p^2 / 3` to `n p^2` flops to form and `p^3 / 3` to factor,
the order of the null space itself. A first version did that, and its refusals at `leverage/1000` (p = 838, six of
seven tangents refused, then recomputed by the null space) made that case 1.12 times slower, so it was dropped.

### Machine and run

The same machine, on 2026-11-19, quiet (no compiler or test process). A is commit `9d66181`, B the range-space branch,
both `gcc-bench` binaries side by side; H is the head of the previous plan (`build/gcc-bench/pb_head.exe`, the
interpolated crossing, built with the instances of that time, whose covariances were generated with the platform
`std::exp`, `std::log` and `std::cos`, so H is a reference for the curve, not a matched A/B). Five rounds, A then B
(then H on `free`), the first discarded; per round:

```
portopt_benchmarks --benchmark_filter='^BM_Solve/max_return/.*/(10|100|1000)(/iterations:1)?$' \
    --benchmark_min_time=0.3s --benchmark_repetitions=3 --benchmark_report_aggregates_only=true
python tests/benchmarks/quick.py run --exe EXE --out FILE
```

The columns and the rule are those of [Portable math](portable-math.md): minimum and median over the four rounds of
each round's median, rounds with B not slower, the standard deviation of the round medians; worse than noise only if
B is slower in both ratios by more than the larger deviation.

### All `max_return` cases (microseconds)

| Case | A min | B min | B/A of the minima | B/A of the medians | rounds with B <= A | A sd | B sd | Worse than noise |
|---|---|---|---|---|---|---|---|---|
| `max_return/free/10` | 11.39 | 10.34 | 0.907 | 0.908 | 4/4 | 0.0747 | 0.069 | no |
| `max_return/free/100` | 368.3 | 179.3 | 0.487 | 0.504 | 4/4 | 7.08 | 12.9 | no |
| `max_return/free/1000` | 1.916e+05 | 6.287e+04 | 0.328 | 0.327 | 4/4 | 2.69e+03 | 620 | no |
| `max_return/long_only/10` | 14.61 | 14.51 | 0.993 | 1.006 | 1/4 | 0.197 | 0.213 | no |
| `max_return/long_only/100` | 2388 | 2330 | 0.976 | 0.998 | 2/4 | 62.9 | 84.6 | no |
| `max_return/long_only/1000` | 1.764e+06 | 1.771e+06 | 1.004 | 1.007 | 0/4 | 6.9e+03 | 9.57e+03 | no |
| `max_return/box/10` | 15.95 | 15.79 | 0.990 | 1.003 | 1/4 | 0.0339 | 0.146 | no |
| `max_return/box/100` | 2583 | 2538 | 0.982 | 0.995 | 2/4 | 44.5 | 120 | no |
| `max_return/box/1000` | 1.821e+06 | 1.806e+06 | 0.992 | 0.994 | 4/4 | 7.48e+03 | 7.55e+03 | no |
| `max_return/leverage/10` | 12.18 | 11.13 | 0.914 | 0.912 | 4/4 | 0.183 | 0.0254 | no |
| `max_return/leverage/100` | 3613 | 3574 | 0.989 | 0.976 | 3/4 | 41.3 | 64.9 | no |
| `max_return/leverage/1000` | 6.294e+06 | 6.272e+06 | 0.997 | 1.000 | 2/4 | 3.95e+04 | 2.35e+04 | no |
| `max_return/cash/10` | 16.21 | 16.26 | 1.004 | 1.005 | 1/4 | 0.108 | 0.0799 | no |
| `max_return/cash/100` | 2383 | 2395 | 1.005 | 0.999 | 1/4 | 17 | 16.2 | no |
| `max_return/cash/1000` | 1.782e+06 | 1.79e+06 | 1.005 | 1.005 | 0/4 | 5.79e+03 | 1.09e+04 | no |

`long_only`, `box` and `cash` reach the cap on a face with fixed weights at every n, and `leverage` at 100 and 1000,
so they run the null space as in A: their solutions are bit for bit those of A, and their rows are noise.

### The curve of `free` (minimum of the round medians, microseconds)

| n | H (interpolated crossing) | A (null space) | B (range space) | B/H |
|---|---:|---:|---:|---:|
| 10 | 5.48 | 11.39 | 10.34 | 1.89 |
| 100 | 146.5 | 368.3 | 179.3 | 1.22 |
| 1000 | 65350 | 191600 | 62870 | 0.96 |

At n = 1000 the regression is gone. At 10 and 100 assets `free` is still slower than with the interpolated crossing,
which needed one Goldfarb-Idnani solve of 1 to 3 iterations from the cached basis there: what remains is the fixed
cost of the tangent and its certificate (the face rows, the least squares of the multipliers, the products with `L`),
not the O(n^3) factorization.

The quick suite in the same rounds: every row noise except `min_variance/long_only/10` (1.007, sd 0.009 us) and
`tracking_error/long_only/10` (1.006); repeated in six rounds of nine repetitions, 1.005 and 1.006, the first again
over its deviation. Neither reaches `parametric_path.cpp` (only `max_return.cpp` includes it). In the two binaries
(`objdump -d`), 60 of the 106 functions of `po::` outside `po::problems` and the harness, `po::qp::solve`, the
`po::linalg` kernels and `po::portable` among them, sit at the same addresses with the same instructions up to call
targets; the others (destructors, copies, `po::riskparity`) and most of `po::problems` were linked after
`parametric_path.cpp`, which grew, and moved. The 0.5% is the placement of unchanged code.
