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
