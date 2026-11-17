# Plan 6 profile 2

Where the time of the slow cases goes after the risk-0 levers of plan 6 (L1 to L7, L6a) and the allocation hygiene
(L8), measured with the routine of [Plan 6 profile 1](profile-1.md), and what is decided for the phase after them.
Tools and commands: [Measuring](measuring.md). The levers and their numbers: [Benchmarks](benchmarks.md).

## Machine and run

The machine of [Benchmarks](benchmarks.md#machine-and-run): AMD Ryzen 7 5700X, Windows 10, g++ 16.2.0 (MSYS2
UCRT64), MSVC 19.51 for the samples. Date 2026-11-17. Desktop session, not isolated. The code is the one of commit
`62f9499` (allocation hygiene before the change to `homogeneous.cpp`, which no sampled case exercises beyond the
`max_sharpe` and `max_diversification` cases of the two "all" sessions).

The three measurements of profile 1, the same cases and commands (the uProf sessions with `--call-graph-mode fpo
--call-graph-type user`, the breakpoint counts of `gdb -batch` with `ignore`, the `BM_Kernel` medians of 5
repetitions), plus the shares are computed from the flat report of each session (self time of the function, shares of
all the samples of the process, the generation of the instance included).

| Session | Filter, `--max-n`, `--repeat` | MSVC solve time | `gcc-bench` time (profile 1) | Samples (s) |
| --- | --- | ---: | ---: | ---: |
| RP target vol 1000 | `target_volatility/1000`, 1000, 1 | 9.49 s | 7.02 s (26.4 s) | 11.1 |
| RP box 1000 | `risk_parity/box/1000`, 1000, 1 | 30.98 s | 19.98 s (85.3 s) | 32.0 |
| RP target return 1000 | `risk_parity/target_return/1000`, 1000, 1 | 25.70 s | 17.35 s (68.2 s) | 27.0 |
| MR leverage 1000 | `max_return/leverage/1000`, 1000, 1 | 15.92 s | 11.57 s (32.3 s) | 17.4 |
| TE box 1000 | `tracking_error/box/1000`, 1000, 1 | 1.66 s | 0.98 s (4.06 s) | 3.3 |
| RP ≤ 100 | `risk_parity`, 100, 30 | 0.081 s per pass (0.168) | - | 2.5 |
| all ≤ 100 | every case, 100, 20 | 0.128 s per pass (0.269) | - | 2.7 |
| all n = 10 | every case, 10, 3000 | - | - | 2.5 |

The `gcc-bench` times of the n = 1000 cases are single solves of one pass of the 35 cases
(`portopt_profile --max-n 1000 --repeat 1`, 74 s of wall time for the 105 cases of the three sizes, 72.1 s of solves
at n = 1000, against the 287.5 s of timed solves that profile 1 estimated for the same 35 cases). Profile 1's times
are from another day: read the ratios as the size of the plan's effect, not as an A/B. The MSVC build is 1.3 to 1.7
times slower than `gcc-bench` on these cases (1.07 to 1.31 in profile 1), so the shares are indicative.

The samples of RP target vol 1000 include `minimum_volatility`, the harness's own solve of the volatility target
(15%), and every n = 1000 session the generation of the instance (`random_covariance`, `orthogonal_factor`: 4.6 to
14.4%, 46% in TE box, whose solve is a second long).

## Hot functions

Share of all the samples (self time). `-`: under 0.4%. `base` is `inverse_transposed_by_rows` (the J0 of
`empty_active_set`); `gauss_newton` is its block loop `add_block`; `compute_directions` is the one of the dense row of
`add_constraint`; `most_violated` is the scan of the explicit rows and the family.

| Function | RP tv 1000 | RP box 1000 | RP tr 1000 | MR lev 1000 | TE box 1000 | RP ≤ 100 | all ≤ 100 | n = 10 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `gauss_newton` (`add_block`) | 43.5% | 24.0% | 29.7% | - | - | 16.7% | 11.4% | - |
| `append_active` (Givens on the basis) | - | 25.8% | 20.3% | 25.0% | 8.1% | 10.9% | 17.4% | 2.6% |
| `directions_from_projection` | - | 18.4% | 14.5% | 21.9% | 6.1% | 8.5% | 12.3% | 1.9% |
| `compute_directions` (dense row) | - | - | - | 36.6% | - | 1.1% | 5.0% | 0.8% |
| `base` | 10.5% | 7.7% | 9.6% | 0.4% | 1.2% | 5.0% | 4.5% | 1.6% |
| Cholesky (`update_trailing_block`) | 8.9% | 5.7% | 6.9% | - | 0.7% | 4.0% | 2.5% | 0.6% |
| `compensated_dot` | 3.9% | 1.8% | 2.3% | - | - | 13.7% | 8.4% | 10.7% |
| `reflect_trailing_columns` (QR) | 8.9% | 1.4% | 2.3% | - | 19.8% | 2.2% | 3.8% | 1.4% |
| `normal_matrix` | - | - | - | - | 12.0% | - | 1.5% | 0.4% |
| `most_violated` | - | 0.8% | - | 0.2% | 0.9% | 0.9% | 1.1% | 2.5% |
| `append_generated` | - | - | - | 1.8% | - | - | - | - |
| operator new (`_Allocate`, deep) | - | - | - | - | - | - | - | 14.1% |
| the harness (instance, `minimum_volatility`) | 14.4% | 4.6% | 5.5% | 8.6% | 46.0% | - | - | - |

What changed against [profile 1](profile-1.md#hot-functions):
- The strided solve `solve_lower_transposed` (20 to 51% of the slow cases) is gone; the base is 5 to 10% of RISK_PARITY
  and 0.4% of MR leverage (J0 is cached along the path).
- `gauss_newton` is the largest function in the MSVC profile of RISK_PARITY, more than in profile 1's (18.8 to 41.0%
  before L4): the MSVC build does not turn the block of L4 into the register-blocked code of GCC. In `gcc-bench` the call
  is 54 ms at n = 1000 (191 ms before), see the next section: read the MSVC share as an upper bound for it.
- The GI iteration (`append_active`, `directions_from_projection`, `compute_directions`) is now the first cost of every
  case with bounds or a family: 44% of RP box, 35% of RP target return, 83% of MR leverage, 14% of TE box, and 35% of
  all the cases at n ≤ 100. `slack` (20.9% of RP box in profile
  1) is below 0.4% in every session: L5.
- `reflect_trailing_columns` is 19.8% of TE box (`reflect` was 20.8%): the QR of the 2000 × 1000 series is 0.365 s in
  `gcc-bench` (1.20 s before), the same share of a solve that is now a quarter of its time.
- At n = 10 `compensated_dot` is 10.7% (12.1% in profile 1: the MSVC build keeps `fma` as a call into the
  `ucrtbase.dll`) and the allocator 14.1% of the samples (15.7%).
- At n ≤ 100 `compensated_dot` is 13.7% of RP, 8.4% of all (10.3% and 6.3% in profile 1).

Callers and the rest: `gauss_newton` is called by `set_model` (the model curvature of every descent step), `linearize`
with the Hessian (the KKT checks and the polish) and `initial_proximal` as in profile 1; every sampled function of the
GI is called by `add_constraint` (100%); `value_change` is 6.0% of RP ≤ 100 and one of its two `moments_of` 47% of
that: see [L9](#l8-l9-and-l10).

## Calls per solve

Breakpoint hits in one solve of the `gcc-bench` build, minus the harness (one `cholesky` in `Data::create`, one
`qp::solve` in `minimum_volatility` for `max_return` and `target_volatility`; the path of `max_return` calls
`qp::solve(problem, cache)`, counted on its own breakpoint). They are the counts of profile 1, to the unit: the levers
change the cost of a call, not the number of calls.

| Case | `qp::solve` | `gauss_newton` | `cholesky` | `householder_qr` | `covariance_times` |
| --- | ---: | ---: | ---: | ---: | ---: |
| `risk_parity/target_volatility/1000` | 29 | 46 | 41 | 7 | 156 |
| `risk_parity/box/1000` | 62 | 71 | 74 | 7 | 208 |
| `risk_parity/target_return/1000` | 67 | 76 | 79 | 7 | 223 |
| `max_return/leverage/1000` | 56 (one basis built) | 0 | 0 | 2 | 0 |
| `tracking_error/box/1000` | 1 | 0 | 0 | 1 | 0 |
| `risk_parity/target_volatility/100` | 31 | 48 | 42 | 7 | 158 |
| `risk_parity/box/100` | 44 | 53 | 55 | 7 | 154 |
| `risk_parity/target_return/100` | 25 | 34 | 36 | 7 | 97 |
| `max_return/leverage/100` | 22 | 0 | 0 | 2 | 0 |

### Time of each kernel inside a case (`gcc-bench`, calls × seconds per call)

Seconds per call, `BM_Kernel` medians of 5 repetitions on the same day (profile 1's in parentheses). At n = 1000:
`empty_active_set` 46.4 ms (356), `gauss_newton` 54.0 ms (191), `cholesky` 26.8 ms (107), `covariance_times` 2.90 ms
(3.94), `qr` 2000 × 1000 0.365 s (1.20), one `qp_iteration` solve 0.375 s (1.85). At n = 100: 42.1 µs (405), 62.7 µs
(152), 31.05 µs (69.6), 29.9 µs (39.4), `qr` 305 µs (774), `qp_iteration` 496 µs (1560).

| Case | Time | GI base | `gauss_newton` | `cholesky` | `covariance_times` | Not priced here |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| RP target vol 1000 | 7.02 s | 1.34 s (19%) | 2.48 s (35%) | 1.10 s (16%) | 0.45 s (6.5%) | 1.64 s (23%) |
| RP box 1000 | 19.98 s | 2.87 s (14%) | 3.83 s (19%) | 1.98 s (9.9%) | 0.60 s (3.0%) | 10.7 s (53%) |
| RP target return 1000 | 17.35 s | 3.11 s (18%) | 4.10 s (24%) | 2.12 s (12%) | 0.65 s (3.7%) | 7.4 s (43%) |
| MR leverage 1000 | 11.57 s | 0.05 s (0.4%) | - | - | - | 11.5 s (99%) |
| TE box 1000 | 0.98 s | 0.05 s (4.7%) | - | - | - | QR 0.365 s (37%), the rest 0.57 s |
| RP target vol 100 | 15.0 ms | 1.31 ms (8.7%) | 3.01 ms (20%) | 1.30 ms (8.7%) | 4.72 ms (31%) | 4.7 ms (31%) |
| RP box 100 | 29.0 ms | 1.85 ms (6.4%) | 3.32 ms (11%) | 1.71 ms (5.9%) | 4.60 ms (16%) | 17.5 ms (60%) |
| RP target return 100 | 16.8 ms | 1.05 ms (6.3%) | 2.13 ms (13%) | 1.12 ms (6.7%) | 2.90 ms (17%) | 9.6 ms (57%) |

The times at n = 1000 are single solves of the pass above; at n = 100 the `target_volatility` and `box` times are the
medians of the quick suite and `target_return` is the best of three of the profile. The unpriced remainder of RP box
and RP target return at n = 1000 (53% and 43%) matches the samples of the GI iteration (`append_active` +
`directions_from_projection`: 44% and 35% in MSVC); in MR leverage 1000 the 9653 iterations cost 1.2 ms each, almost
all of it `compute_directions` on the dense generated rows, the Givens rotations and the axpys of the directions.
`covariance_times` at n = 100 went from 8 to 18% to 16 to 31% of RP: it dropped less (1.3 times) than the rest.

## Where the time went, profile 1 to profile 2

| Case | Profile 1 | Profile 2 | Ratio | What moved |
| --- | ---: | ---: | ---: | --- |
| RP target vol 1000 | 26.4 s | 7.02 s | 0.27 | base 10.3 to 1.34 s, `gauss_newton` 8.8 to 2.48 s, `cholesky` 4.4 to 1.10 s |
| RP box 1000 | 85.3 s | 19.98 s | 0.23 | base 22.1 to 2.87 s, `gauss_newton` 13.6 to 3.83 s, `slack` and `compute_directions` (L5) |
| RP target return 1000 | 68.2 s | 17.35 s | 0.25 | same |
| MR leverage 1000 | 32.3 s | 11.57 s | 0.36 | base (19.9 s, L1 and L7), allocation (1.3 s, L8) |
| TE box 1000 | 4.06 s | 0.98 s | 0.24 | QR 1.20 to 0.365 s, `slack` (L5) |
| RP target vol 100 | 34.1 ms | 15.0 ms | 0.44 | base, `gauss_newton`, Cholesky; `covariance_times` now the largest term (31%) |
| RP box 100 | 73.8 ms | 29.0 ms | 0.39 | the GI iteration now 60% |
| MR leverage 100 | 16.6 ms | 8.1 ms | 0.49 | base (54% of it) |

All 35 cases at n = 1000 together: 287.5 s estimated for profile 1's baseline against 72.1 s of solves now (one
solve of each; the baseline numbers are single runs of another day, so 4 times is the size of the effect, not an
interval).

## L8, L9 and L10

**L8, allocation hygiene.** Done where the counter showed it, see [Benchmarks](benchmarks.md#l8-allocation-hygiene-where-the-counter-shows-it):
29 054 to 13 348 allocations per pass of the 70 cases at n ≤ 100, 0.63 to 0.87 on the small cases with bounds, 0.5 to 0.9
on the gross exposure cases at n = 1000, bit for bit identical. The allocator is
still 14.1% of the MSVC samples at n = 10 (15.7% before), mostly RISK_PARITY, whose 1600 allocations per solve are
spread over about 40 sites: a proposal below.

**L9, `Moments` cache in the Armijo halvings: stays dropped.** `value_change` is 6.0% of the RP ≤ 100 samples (MSVC,
4.5% in profile 1) and the cache would save at most the recomputation of one of its two `moments_of`, 47% of that: 2.8%
of the samples: under what the quick suite resolves (3 to 5% between two runs of the same binary). Not visible at
n = 1000. The number is, as in profile 1, below the noise, so the lever stays out.

**L10, `-flto` on the core: not adopted, decided with these numbers.** Verification of `-ffp-contract=off` at the link,
on a `-march=native` build (where `-mfma` makes a contraction visible as `vfm*`), `objdump -d` of the executables:

| Build | `vfm*` instructions |
| --- | --- |
| `-march=native`, no LTO, objects of `portopt_core` ([Benchmarks](benchmarks.md#no-accidental-fma-guard-3)) | 3 (`two_product_error` and its two inlinings) |
| `-flto`, `portopt_profile` (core and benchmark TUs all with the flag) | 1, in `compensated_dot`, where the `std::fma` is inlined |
| `-flto`, `libportopt.dll`, `portopt_capi` without the flag (as in `CMakeLists.txt`) | 25: `po_solve` 8, `po_sharpe_ratio` 4, `po_utility` 4, `po_tracking_error` 3, `po_portfolio_return` 3, `po_diversification_ratio` 2, and 1 in `compensated_dot` |
| `-flto`, `libportopt.dll`, the flag on every TU | 1 |
| control, toy `cholesky.cpp` with the flag and `main.cpp` without, `-flto` | 26, all in `main` (the inlined Cholesky); both with the flag: 0 |

The flag of each translation unit is kept at the link (no `vfm*` in `gauss_newton`, `cholesky` or the QP code of the
LTO executable, with no flag at the link), but a function that is inlined into a unit compiled without it is compiled
under the caller's setting: with LTO, `portopt_capi` (which today does not carry the flag, because it does no
arithmetic) would contract the core code that GCC inlines into the C entry points. Adopting LTO therefore needs
`-ffp-contract=off` on `portopt_capi` and on every target that includes inline core code, and a guard for it (no
`vfm*` outside `std::fma` in the `native-release` library).

The gain, generic builds, `gcc-bench` against the same tree with `CMAKE_INTERPROCEDURAL_OPTIMIZATION=ON`, quick suite
twice interleaved: n = 10 cases 0.87 to 1.00 (min_variance/long_only 0.90 and 0.89, max_return/long_only 0.89 and 0.91,
utility/box 0.89 and 0.96, max_diversification/box 0.87 and 0.99, risk_parity/box 0.94 and 0.90, target_volatility 0.98
and 0.93), the geometric mean about 0.93; n = 100 cases 0.99 to 1.15 (no gain: the kernels do not cross a unit); the
unchanged kernels 0.87 to 1.06. The 70 cases at n ≤ 100 are bit for bit identical to the non-LTO build. The cost: a
link of about 25 s that is serial in GCC 16.2 ("using serial compilation of 9 LTRANS jobs"), a flag for the library
target, and a guard to keep. Decision: not adopted in the default build in this task; it buys 5 to 10% of solves of 3
to 5 µs and nothing at n ≥ 100, where all the time is in kernels that do not cross a unit. It is a small separate task
if the author wants it for callers of many small problems. MSVC's `/GL` and `/LTCG` were not looked at.

## Proposals for the author

None of these is implemented. Each has the measured potential on this machine; an ISA or reassociation change needs its
own plan (the risk class of each is in [Measuring](measuring.md#profile-and-comparison-of-solutions)).

**1. ISA multiversioning (class 0 for elementwise kernels).** Same tree, `-march=native` (AVX2 and FMA on this Zen 3,
`-ffp-contract=off`, so no contraction) against generic, quick suite twice interleaved: n = 100 solves 0.75 to 0.87
(min_variance/long_only 0.83, max_return/long_only 0.79, max_sharpe/long_only 0.80, utility/box 0.83,
tracking_error/long_only 0.75, max_diversification/box 0.80, risk_parity/budget 0.78, target_volatility 0.87, box
0.87), kernels `cholesky/100` 0.73, `qr/100` 0.74, `qp_iteration/100` 0.81, `empty_active_set/100` 0.90,
`gauss_newton/100` 0.97; n = 10 cases 0.92 to 1.03. All 70 cases at n ≤ 100 are bit for bit identical to the generic
build, so the gain is free of numerical risk; at n = 1000 it was not measured. It needs a runtime dispatch (MinGW has
no `ifunc`, so `target("avx2")` and `__builtin_cpu_supports` by hand; MSVC another mechanism) and keeps the generic
binary as the shipped default.

**2. FMA-aware Dekker (class 0, needs 1 or a compile-time choice).** In a native build the Dekker product of L6a
is 2.5 times slower than the inline `std::fma` ([L6a](benchmarks.md#l6a-dekkers-product-in-the-compensated-dot)).
Scratch build (not committed) with `two_product_error` returning `std::fma(a, b, -product)` under `__FMA__`:
`compensated_dot/100` 0.38 of the generic one (and 0.40 of the native Dekker build), `covariance_times/100` 0.39,
risk_parity/target_volatility/100 0.70 against 0.87 for native Dekker and 1.00 generic (second run 0.67), risk_parity/box/100
0.76 against 0.87, target_volatility/10 0.80 against 0.92; the 70 cases are bit for bit identical to the generic build
(the fma is the exact error in the Dekker range, and outside it L6a already called it). Together with proposal 1:
about 0.70 of RP target vol 100 and 0.76 of RP box 100, from `covariance_times`, which is 16 to 31% of RISK_PARITY at
n = 100.

**3. Class-B kernels: reassociated reductions.** Not measured in the solver (they change the rounding; hazard 1 of
the [research](../research/2026-11-16-performance-method.md#5-hazards): the band of RISK_PARITY with cond 1e8 has
margin 6). A scratch dot product (`dotvar.cpp`, random vectors, best of 5): the ordered sum of the core costs 0.41 ns
a term at n = 100 and 0.64 at n = 1000 (a chain of dependent additions); four accumulators 0.147 and 0.174, eight
0.154 and 0.137 (generic), 0.108 and 0.166 against 0.093 and 0.100 with AVX2: 2.8 to 6.4 times. The first use is
`compute_directions` of a dense row (n dots of n per iteration: 0.64 ms of the 1.2 ms of an MR leverage iteration at
n = 1000, so about 0.5 ms saved per iteration, 11.6 s toward 6 to 7 s, by arithmetic, not by experiment), then
`slack` of dense rows and the `dot` of `append`/`refine`. The compensated dot (Dot2) is not a candidate: its
accuracy depends on the order. It needs the A/B of weights inside the band of `expect_user_kkt` and
`expect_risk_parity_kkt`, with status and iterations unchanged.

**4. RISK_PARITY: a warm start of the QP (class A).** The trace of the active sets of the consecutive QPs of a descent
(a scratch counter in `step_at`, not committed) on RP box 100, target return 100 and box 1000: the iterations of a
solve equal the number of its active rows (ratio 1.00, 1.01 and 1.02: the dual method adds a row per iteration
and almost never drops one), and consecutive solves share 98% of their active inequality rows. So the iterations
are not wasted and a warm start cannot save them one by one: it would have to build the guessed active set at once (one
solve of the KKT system) and check it. The ceiling is the GI iteration of RISK_PARITY: 53% of RP box 1000 and 43% of RP
target return 1000 (the unpriced remainder above; none in `target_volatility`, whose QP has one equality), 60 and
57% at n = 100, if the iterations were free, which they will not be. It changes the path of the QP (rounding of
the solution) and `iterations`, so it needs its own study and the 576 forms of `RiskParitySolveSweep`.

**5. RISK_PARITY: a workspace for the descent (class 0).** The 1600 allocations per solve of `risk_parity/target_volatility`
at n = 10 are 25% of its 0.22 ms (35 ns each), 0.4% at n = 100: a workspace passed through `linearize`, the KKT
and the polish (about 40 sites) removes at most that, only at n ≈ 10. Worth it only for callers that solve many small
RISK_PARITY problems.

**6. Not proposed: more `qp::solve` work for `max_return/leverage`.** After L1, L5 and L7 the 9653 iterations of MR
leverage 1000 take 1.2 ms each; the dense family rows are the cost (proposal 3), not their count.

## What is left

The order of the remaining cost, from the shares and the call counts above, at n = 1000 (RP = RISK_PARITY with
constraints, the slow cases):

1. **The GI iteration** (`append_active`, `directions_from_projection`, `compute_directions`): 53% of RP box, 43% of
   RP target return, 99% of MR leverage, and 40 to 60% of RP at n = 100. The Givens rotation and the axpy over the
   basis are elementwise (proposal 1); the dense dots are reductions (proposal 3); the number of iterations is
   the structure of the active set (proposal 4).
2. **`gauss_newton`**: 35%, 19% and 24% of RP target vol, box and target return, 54 ms a call at n = 1000.
3. **The GI base**: 14 to 19% of RP (46 ms a call: n³/6 multiply-subtracts, plus the transposition and the copy of
   J0).
4. **`cholesky`**: 10 to 16% of RP (26.8 ms a call).
5. **`covariance_times`**: 16 to 31% of RP at n = 100, 3 to 6% at n = 1000 (proposal 2).
6. **QR of the TRACKING_ERROR series**: 37% of TE box 1000 (0.365 s of 0.98 s).

At n = 10 the allocator (about 25% of RISK_PARITY's time) and `compensated_dot` (10% of the samples).
