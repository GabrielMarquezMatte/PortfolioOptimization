# Python and C layers against the solve

What the layers above `po::solve` cost: the C interface (`include/portopt/portopt.h`, one guard and one handle per
call) and the Python package (`ctypes` and `numpy` over it), against the time of the solve they wrap. No other
library is measured, so no comparison with one is made; the file keeps the name of its model in the author's
VolEngine (`python-vs-arch.md`), whose table compared a package with `arch`. Every number comes from the one machine
below; compare only runs from the same machine.

Tools:

- `tests/benchmarks/capi/bm_capi.cpp`: `BM_CApi/*` in `portopt_benchmarks`, the C calls against the C++ calls under
  them, on the instances of [Measuring](measuring.md#the-instances).
- `python/benchmarks/bench_overhead.py`: the same steps through `ctypes`, one row per step.
- `portopt_profile --via capi`: every benchmark case solved through the C interface, for `compare-solutions`
  ([Measuring](measuring.md#profile-and-comparison-of-solutions)).

## Machine and run

| Item | Value |
| --- | --- |
| Machine | The one of [Benchmarks](benchmarks.md#machine-and-run): AMD Ryzen 7 5700X, Windows 10, g++ 16.2.0 |
| Build | `gcc-bench` (generic x86-64, `-O3 -DNDEBUG -ffp-contract=off`), code at `6419150` |
| Python | 3.14.4, numpy 2.5.4, `portopt` 0.1.0 over that build's `libportopt.dll` (`PORTOPT_LIB`) |
| Date | 2026-11-17 |
| C layer | `BM_CApi`: median of 21 repetitions interleaved at n = 10 and 100, of 9 at n = 1000, of 5 for the calls and the data handle at n = 1000 |
| Python layer | `bench_overhead.py --repeats 15`, run twice; the second run is tabled, the range of the two in the text |

The run was not isolated: a normal desktop session, no pinned cores and no fixed frequency. The first `BM_CApi` run
of the same day moved the full-cycle rows by 2 to 3 times (a noisy minute of the desktop) and was repeated with the
repetitions interleaved; the spread inside the tabled runs is 5 to 20% for the sized rows and up to 60% for the
solves at n = 1000, whose few iterations catch every burst of the desktop. A difference smaller than that is not a
result.

## The C layer

The case is `MIN_VARIANCE` with the weights bounded below by 0 (`long_only`) and, for the fixed cost of a solve that
ends at once, `free`, the two cheapest solves of their size, so the ratio of the C layer to the solve is the largest
one it can have. `c` is the C entry point (`po_solve` and the free of its result), `C++` the call under it
(`po::solve` and the destruction of its result).

| Measure | n = 10 | n = 100 | n = 1000 |
| --- | ---: | ---: | ---: |
| `po_abi_version()` (a call, nothing else) | 0.98 ns | - | - |
| `po_data_asset_count()` (the cheapest call through the guard) | 3.5 ns | - | - |
| Spec: `po_spec_create`, `po_spec_set_bounds`, `po_spec_free` | 42.8 ns | - | - |
| `po_data_create` and its free, `c` | 0.727 µs | 65.7 µs | 35.7 ms |
| `po::Data::create` and its destruction, C++ | 0.656 µs | 66.8 µs | 33.0 ms |
| `po_solve`, `free` case, `c` | 2.15 µs | 79.4 µs | 53.7 ms |
| `po::solve`, `free` case, C++ | 2.24 µs | 79.8 µs | 53.3 ms |
| `po_solve`, `long_only` case, `c` | 3.37 µs | 240 µs | 198 ms |
| `po::solve`, `long_only` case, C++ | 3.31 µs | 249 µs | 165 ms (noisy) |
| The 11 readers of a result into a `po::Result` (status, weights, iterations, 4 multipliers, bound multipliers, residual) | 0.128 µs | 0.234 µs | 1.24 µs |
| Whole cycle (spec, data, `po_solve`, readers, frees), `c` | 4.35 µs | 339 µs | 204 ms |
| Whole cycle, the C++ calls under it | 4.10 µs | 334 µs | 212 ms (noisy) |

The C interface adds nothing that these runs resolve: `po_solve` and `po::solve` differ by less than their spread
at every size (`free` at n = 10: 2.15 against 2.24 µs, spread 13 to 20%), and so do the data handles from n = 100 up
and the whole cycle. What it does add is measurable only at n = 10, in the tens of nanoseconds: the guard (3.5 ns), one
`new` of the handle and the move of the `Result` into it (`po_data_create` 0.07 µs over `Data::create`, 10%), and
the readers, 0.13 µs at n = 10 and 1.2 µs at n = 1000 (the copy of two vectors of n doubles, 2 allocations of the
`std::vector` of the harness). The fixed cost of the C layer is under 0.3 µs per cycle, at most 6% of a 4 µs solve
and nothing at n ≥ 100.

`po_data_create` is the one call that is not fixed cost: it copies the covariance and factors it
(Cholesky), 66 µs at n = 100 and 33 ms at n = 1000, as `Data::create` does. A caller that solves many problems on
the same data builds the handle once.

## The Python layer

`MIN_VARIANCE`, long only, on a covariance of condition 1e2 and means in [0.02, 0.12] drawn with numpy (the size and
conditioning of the C++ instance, not its bits: the raw solve is 4.7 µs here against 3.4 µs there). Each row is the
median of 15 samples; a sample of a fast row is the mean of 5000 calls, of the solves 1000 / n calls.

| Measure (microseconds, median) | n = 10 | n = 100 | n = 1000 |
| --- | ---: | ---: | ---: |
| no-op Python call | 0.027 | 0.027 | 0.029 |
| `ctypes` call, no arguments (`po_abi_version`) | 0.135 | 0.134 | 0.140 |
| `ctypes` call through the guard (`po_data_asset_count`, prebuilt `byref`) | 0.466 | 0.444 | 0.453 |
| the same, through `read_scalar` (a new `c_size_t`, `byref`, `check`, `.value`) | 0.810 | 0.770 | 0.910 |
| `po_data_create` and its free through `ctypes` | 2.66 | 60.8 | 35 553 |
| `po_solve` and `po_result_free` through `ctypes` | 4.71 | 192.7 | 165 849 |
| `Spec(...)` and its free | 5.81 | 5.76 | 7.48 |
| `Data(mu, cov)` and its free | 9.13 | 69.1 | 32 572 |
| `read_result` (11 `ctypes` calls, 2 numpy arrays, the dataclass) | 15.8 | 15.8 | 16.4 |
| `solve(spec, data)` | 22.3 | 217.9 | 162 286 |

The two runs agree to 5% on every row but `po_data_create` and `po_solve` at n = 1000 (33.0 to 43.9 ms and 166 to
181 ms: the native work, in the noise of the desktop) and the raw solve at n = 100 (193 and 221 µs).

What the Python layer costs, besides the native work:

- A `ctypes` call is 0.13 µs with no arguments and 0.45 µs through the guard with two arguments (the C call is 3.5 ns of
  it); through `read_scalar`, which builds the output cell on every call, 0.8 µs. The calls are the cost, not the
  data: there is no copy of the arrays (`np.ascontiguousarray` returns a contiguous float64 array as it is, and its
  pointer is passed).
- `Spec(...)` costs 5.8 µs: the dataclass, the `canonical_bounds` normalization of the bounds and four calls.
- `Data(mu, cov)` costs 6.5 to 8 µs over the native `po_data_create` at n = 10 and 100 (the validations of shape and the
  four calls of the constructor); at n = 1000 the difference is inside the noise of the 33 ms.
- `solve()` costs 17.6 µs over `po_solve`: `read_result` is 15.8 to 16.4 µs of it (11 calls through `read_scalar`, 8.6
  µs, and the rest the two `np.empty` and the read-only flags, the `SolveStatus`, the frozen dataclass), the handle
  of the result the rest. It does not depend on n: the vectors are copied once, 8 KB at n = 1000.

| Python overhead per `Spec`, `Data` and `solve` (30 µs) | n = 10 | n = 100 | n = 1000 |
| --- | ---: | ---: | ---: |
| `po_solve` of the case above | 4.7 µs | 193 µs | 166 ms |
| Overhead as a share of that solve | 640% | 16% | 0.02% |
| Overhead of `solve()` alone, with the `Spec` and `Data` reused (17.6 µs) | 370% | 9% | 0.01% |

The cases above are the shortest solve of each size that has a bound. At n = 10, `RISK_PARITY` with a target
volatility takes 229 µs (the 17.6 µs of `solve()` is 8% of it) and `MAX_RETURN` with a leverage 44 µs. At n = 100 the
cases that end at once (`free`, `leverage`, `cash`) take 70 to 110 µs, where the 30 µs is 27 to 43%; the cases that
iterate take 0.23 to 29 ms, 13% of the cheapest and under 7% of the others. At n = 1000 it is negligible for every
case, the cheapest `po_solve` being 49.8 ms.

## Reading

- The cost of a call into the library is in the Python layer, not in the C layer: 0.45 µs of `ctypes` against 3.5 ns
  of guard, 17.6 µs of `solve()` against 0.13 µs of readers. The C layer is not worth optimizing.
- The layer matters for one use, many solves of problems of about 10 assets, where the Python overhead is more than
  the solve itself. For one solve at n ≥ 100 it is not visible; for n = 1000 it is not measurable.
- What would reduce the fixed cost, none of it done: building the `ctypes` output cells once instead of per call in
  `read_scalar` (the 11 reads are 8.6 of the 15.8 µs), filling the `Result` in a single pass without the
  frozen-dataclass `__init__`, and reusing the `Spec` and `Data` between solves (already possible: they are
  immutable handles). A reader that returns the whole result in one call would need a new C function, a change of the
  C interface that the numbers above do not justify.

Reproduce, in the `python/.venv` environment, with the library of the build to measure:

```
cmake --preset gcc-bench && cmake --build --preset gcc-bench
build/gcc-bench/portopt_benchmarks --benchmark_filter='^BM_CApi' --benchmark_repetitions=21 --benchmark_enable_random_interleaving=true --benchmark_report_aggregates_only=true
PORTOPT_LIB=$PWD/build/gcc-bench/libportopt.dll python python/benchmarks/bench_overhead.py --repeats 15
```
