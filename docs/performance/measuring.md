# Measuring

How to answer, in seconds, where the time goes and whether a change helped. The official numbers come from the
generic `gcc-bench` build (no `-march=native`), because the shipped binaries are generic. The method is the one of
the author's VolEngine, with PortOpt's own content: instances of portfolio problems, the counters the library
exposes, and a comparison of solutions instead of a log-likelihood.

Run every executable with `/c/msys64/ucrt64/bin` first on `PATH` (Git's mingw64 `libstdc++` lacks symbols and the
executable dies on start with a heap-corruption status):

```
export PATH=/c/msys64/ucrt64/bin:$PATH
cmake --preset gcc-bench
cmake --build --preset gcc-bench
```

The benchmarks are opt-in (`PO_BUILD_BENCHMARKS`, off by default and on in the `gcc-bench` and `msvc-profile`
presets), so the default builds and `ctest --preset gcc-release` do not see them. They fetch Google Benchmark
v1.9.1. Two executables come out:

| Target | What it is |
|---|---|
| `portopt_benchmarks` | Google Benchmark suite: `BM_Solve/*` and `BM_Kernel/*` |
| `portopt_profile` | One pass over every case, with the per-case table and the dump of solutions |

Both link `portopt_core` directly (no C API) and are compiled with `-ffp-contract=off`, like the core: the
compensated products need every `a * b + c` rounded as written. They also replace the global `operator new` to count
allocations; no other target does. `ctest --preset gcc-bench` runs the two self-tests: `portopt_benchmarks
--selftest` (the allocation counter against allocations whose number is known) and `quick.py selftest` (the
comparison code on synthetic files).

## The instances

`tests/benchmarks/instances.hpp` and `cases.hpp`. An instance is fixed by its size and a few constants, so two runs and
two builds see the same numbers:

- covariance `random_covariance(n, seed 1, condition 1e2)` of `tests/unit/random_covariance.hpp`; mu uniform in
  [0.02, 0.12]; risk-free rate 0.01 in the `cash` cases;
- the generation is O(n³) and takes seconds at n = 1000, so each covariance and each `Data` is built once per process,
  outside the timed loop (the first benchmark that needs a size pays for it);
- `tracking_error`: T = 2n observations of scale 2% and a benchmark near their cross-sectional mean.

Cases are `BM_Solve/<objective>/<variant>/<n>`, `po::solve` end to end (the C++ entry point under the C API):

| Objective | Variants | n |
|---|---|---|
| `min_variance` | `free`, `long_only`, `box`, `leverage`, `cash`, `target_return` | 10, 100; 1000 for `long_only` |
| `max_return` | `free`, `long_only`, `box`, `leverage`, `cash` (the volatility cap is 1.5 times the minimum volatility of the feasible set) | idem |
| `max_sharpe`, `utility`, `tracking_error`, `max_diversification` | `free`, `long_only`, `box`, `leverage`, `cash` | idem |
| `risk_parity` | `budget` (the Spinu point is the answer), `target_return`, `target_volatility` (1.5 times the minimum volatility), `box` (a cap of 1.5 / n that cuts the Spinu point) | 10, 100; 1000 for `budget` and `box` |

`box` is [0, max(0.3, 2/n)] on every weight, `leverage` a gross exposure limit of 1.5 with free weights, `long_only`
a lower bound of 0. `IMPORTANCE` and `KELLY` share the builder of `UTILITY` and are not benchmarked. The n = 1000
cases are registered with `Iterations(1)` (their name ends in `/iterations:1`); a variant other than the ones in the
table at n = 1000 takes minutes per solve and has no case.

Every case asserts `status == OPTIMAL` and finite weights, before the loop (on a first solve) and on every solve
inside it, and ends with `SkipWithError` otherwise: a benchmark that measures a solve that failed would publish the
time of the failure.

## Counters

| Counter | Where | Source | State |
|---|---|---|---|
| `iterations` | `BM_Solve`, profile | `Result::iterations`, counted in a method-specific way (QP iterations; for `max_return` those of every QP of the path; for `risk_parity` the Spinu Newton steps, the steps of every descent and the Newton steps on the KKT system) | real |
| `allocations` | `BM_Solve`, profile | calls of the global `operator new` per solve | real |
| `status`, `residual` | `BM_Solve`, profile | `Result::status` (0 is OPTIMAL); `Result::residual` of `risk_parity` | real |
| `qp_iterations`, `seconds_per_qp_iteration` | `BM_Kernel/qp_iteration` | `qp::Solution::iterations` | real |
| `qp_solves` | - | no result field or accessor counts the calls of `qp::solve`; observing it needs a counter in `src/` | **unavailable** |
| `descent_iterations`, `polish_steps` | - | folded into `Result::iterations` for `risk_parity`, not separable without touching `src/` | **unavailable** |

The `iterations` counter shares its JSON key with Google Benchmark's own iteration count; in the JSON the counter is
written last and wins. Aligned (over-aligned type) allocations keep the library's allocation function and are not
counted; nothing measured here uses them.

`BM_Kernel/<kernel>/<n>` for n in 10, 100, 300, 1000: `cholesky`, `solve_lower`, `qr` (a 2n x n matrix, the shape of the
`tracking_error` series), `qp_iteration` (a long-only QP over the covariance, solved by `qp::solve`: the time of one
solve and of one iteration), `compensated_dot`, `covariance_times`, `gauss_newton`, `empty_active_set`. `solve_lower`
copies its right-hand side into the work vector on every iteration, because it solves in place.

## Quick suite (about 35 to 45 s)

```
python tests/benchmarks/quick.py run --out build/gcc-bench/quick.json
```

Three repetitions of a fixed filter (`QUICK_FILTER` in `tests/benchmarks/quick.py`), with `--benchmark_min_time=0.3s`:
one case of each objective at n = 10 and 100 (`min_variance/long_only`, `utility/box`, `max_return/long_only`,
`max_sharpe/long_only`, `tracking_error/long_only`, `max_diversification/box` and the three `risk_parity` variants
`budget`, `target_volatility` and `box`) and the eight kernels at n = 100. It prints its own wall time at the end.
Measured on 2026-11-16 on the development desktop (AMD Ryzen 7 5700X, Windows 10, MinGW-GCC 16.2.0, `gcc-bench`,
a desktop session that is not isolated): 44.9 s with `--benchmark_min_time=0.4s`, then 44.0 s and 36.1 s with 0.3 s,
the last two one right after the other. Two consecutive runs of the same binary moved some medians by up to 3x
(`tracking_error/long_only/10`: 19.5 and 5.76 microseconds; `max_diversification/box/10` 18.1 and 6.6), well beyond the
`noise` mark of `compare`, which only sees spread inside one run: the first run after a build or a pause is the slow
one, so discard it or run the suite twice before an A/B.

## A/B

Measure before and after in the same session, one right after the other. On this machine, runs on different days
move some medians by more than their spread, so never compare against an old table.

```
python tests/benchmarks/quick.py run --out build/gcc-bench/quick_a.json
# change, rebuild
python tests/benchmarks/quick.py run --out build/gcc-bench/quick_b.json
python tests/benchmarks/quick.py compare build/gcc-bench/quick_a.json build/gcc-bench/quick_b.json
```

`compare` prints the medians, the ratio B/A and `noise` when |B - A| is within the larger standard deviation. For a
gain below 5%, repeat the run interleaved or use Google Benchmark's `compare.py` with at least 9 repetitions: three
repetitions and a deviation rule are a screen, not a test.

## Profile and comparison of solutions

```
build/gcc-bench/portopt_profile [--filter TEXT] [--max-n N] [--repeat K] [--json PATH]
```

Runs every case once through `po::solve` (`--max-n` limits the size, 100 by default; `--filter` is a substring of the
case name; `--repeat` takes the best of K timings, 3 by default; on Git Bash a filter must not start with `/`, which
the shell turns into a path) and prints, per case, the status, the solver iterations, the time and the allocations.
With `--json` it writes the solutions: status, iterations, objective value (computed by `metrics::`; the residual
for `risk_parity`), weights, the four named multipliers, bound multipliers and residual, all as hex-float (`%a`, exact).

Compare two dumps, for example before and after a change to a kernel:

```
build/gcc-bench/portopt_profile --max-n 100 --json build/gcc-bench/solutions_a.json
# change, rebuild
build/gcc-bench/portopt_profile --max-n 100 --json build/gcc-bench/solutions_b.json
python tests/benchmarks/quick.py compare-solutions build/gcc-bench/solutions_a.json build/gcc-bench/solutions_b.json
```

It prints the number of cases that are bit for bit identical (status, objective, weights, multipliers and residual
equal as hex-float text), the worst relative weight difference (`max|dw| / max|w|`) and relative objective
difference with the case that has each, the cases whose status or iteration count changed, the ratio of iterations
and the cases present on one side only. Cases are matched by name, so the order does not matter. With
`--require-identical` the exit status is 1 unless every case is present on both sides and identical.

The acceptance of a change follows from its risk class: a change that keeps the sequence of operations of every
element (risk 0) must give all cases identical; a change that rounds differently must keep weights and objectives
inside the band of the instance (`expect_user_kkt`, `expect_risk_parity_kkt`) with status and iterations unchanged;
a change of the algorithm's path needs its own study.

## Sampling profile (AMD uProf)

AMD uProf on Windows does not symbolize the MinGW-gcc build: `AMDuProfCLI report` crashes on a session of a DWARF
binary, and GCC's `-gcodeview` fails in the assembler. The `msvc-profile` preset builds the same two executables with
PDB symbols (`/O2 /Ob2 /DNDEBUG /Zi`, linker `/DEBUG /OPT:REF /OPT:ICF`), and serves only to locate hot spots. MSVC
code is slower than GCC's, so read the shares as indicative; how much a change gains always comes from `gcc-bench`.

```
cmake --preset msvc-profile
cmake --build --preset msvc-profile --target portopt_profile
"C:/Program Files/AMD/AMDuProf/bin/AMDuProfCLI.exe" collect --config tbp -o build/uprof build/msvc-profile/Release/portopt_profile.exe --filter risk_parity --max-n 100
"C:/Program Files/AMD/AMDuProf/bin/AMDuProfCLI.exe" report -i build/uprof/<session> --bin-path <abs path to build/msvc-profile/Release> --symbol-path <same>
```

The report is `report.csv` in the session directory; without `--bin-path` and `--symbol-path` the functions show as
bare addresses. The `msvc-profile` preset has not been built on the machine of this document: only the GCC
toolchain was used to write the suite.

## Assembly

`objdump -d -C -M intel --no-show-raw-insn` of an object of the `gcc-bench` build (or `-S` with `-g`) and
`-fopt-info-vec-missed`. Look for calls that should be inlined, allocations in a loop, copies, divisions and
transcendental functions that can be hoisted, and, in a generic build, `call fma`: `std::fma` is a call into `libm`
without `-mfma`. Findings go to `docs/performance/`, never into a comment.

## Guards that stay green during any optimization

1. `ctest --preset gcc-release` in full (this author runs `gcc-release` only, not `gcc-debug`): the fuzz batteries and
   the `style` test included.
2. `quick.py compare-solutions ... --require-identical` over the profile dump for a change of risk 0.
3. No accidental FMA: in a `native-release` build, no `vfm*` instruction in the objects of `portopt_core` outside the
   functions that call `std::fma` on purpose.
4. No band of `expect_user_kkt` or `expect_risk_parity_kkt` loosened to accommodate an optimization.
