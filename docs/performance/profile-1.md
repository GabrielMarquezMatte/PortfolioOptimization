# Plan 6 profile 1

Where the time of the slow cases goes, measured before any plan 6 optimization, and which of the hot-spot hypotheses
of the performance research (H1 to H10) hold. Tools and commands: [Measuring](measuring.md). Times of every case:
[Plan 6 baseline](baseline.md).

## Machine and run

The machine of [Benchmarks](benchmarks.md#machine-and-run): AMD Ryzen 7 5700X (8 cores, 16 hardware threads), Windows
10, g++ 16.2.0 (MSYS2 UCRT64). Code at `4cc0b82`. Date 2026-11-16. Desktop session, not isolated.

Three kinds of measurement, each answering a different question:
- **Shares of time**: AMD uProf 5.0 time-based sampling (1 ms, call stacks by unwinding) of `portopt_profile` built
  by the `msvc-profile` preset (MSVC 19.51, `/O2 /Ob2 /Zi`), because uProf does not symbolize the GCC build. The
  MSVC build is 1.07 to 1.31 times slower than `gcc-bench` on the cases sampled (below), so the shares are
  indicative. uProf reports the session as running under Hyper-V (Windows virtualization-based security); sampling
  worked.
- **Calls per solve**: a debugger breakpoint count on the `gcc-bench` build of `portopt_profile` (`gdb -batch`,
  `ignore` on each breakpoint, one solve per case). The library has no counter of `qp::solve` calls (`qp_solves` is
  unavailable, [Measuring](measuring.md#counters)).
- **Seconds per call**: the `BM_Kernel` medians of [Plan 6 baseline](baseline.md#kernels) (`gcc-bench`). Calls ×
  seconds per call gives the `gcc-bench` time of a kernel inside a case. It is an estimate (isolated kernel median ×
  call count): the kernels run on the covariance of the same n (`qr` on a 2n × n matrix), but in a loop with warm
  caches, and the Jacobian of `gauss_newton` changes from call to call (not its size).

```
"C:/Program Files/AMD/AMDuProf/bin/AMDuProfCLI.exe" collect --config tbp --call-graph-mode fpo --call-graph-type user -o build/uprof/<name> build/msvc-profile/Release/portopt_profile.exe --filter <case> --max-n <n> --repeat <k>
"C:/Program Files/AMD/AMDuProf/bin/AMDuProfCLI.exe" report -i build/uprof/<name>/<session> --bin-path <abs path to build/msvc-profile/Release> --symbol-path <same> --cutoff 40 --detail -g --report-output build/uprof/<name>/callgraph.csv
```

`--src-path <abs path to src>` with `--detail` adds the samples per source line; `--cutoff 40` lists 40 functions
instead of 10.

## The cases sampled

| Session | Filter, `--max-n`, `--repeat` | MSVC solve time | `gcc-bench` time | Samples (s) |
| --- | --- | ---: | ---: | ---: |
| RP target vol 1000 | `target_volatility/1000`, 1000, 1 | 34.5 s | 26.4 s | 36.4 |
| RP box 1000 | `risk_parity/box/1000`, 1000, 1 | 106.2 s | 85.3 s | 106.3 |
| RP target return 1000 | `risk_parity/target_return/1000`, 1000, 1 | 78.3 s | 68.2 s | 79.2 |
| MR leverage 1000 | `max_return/leverage/1000`, 1000, 1 | 38.1 s | 32.3 s | 39.6 |
| TE box 1000 | `tracking_error/box/1000`, 1000, 1 | 4.32 s | 4.06 s | 6.0 |
| RP ≤ 100 | `risk_parity`, 100, 30 | 0.168 s per pass | 0.144 s per pass | 5.2 |
| all ≤ 100 | every case, 100, 20 | 0.269 s per pass | - | 5.8 |
| all n = 10 | every case, 10, 3000 | - | - | 3.3 |

The samples of an n = 1000 session include the generation of the instance (`random_covariance`,
`orthogonal_factor`, 2 to 7% of the samples; 28% for TE box); the shares below are of all the samples of the
process.

## Hot functions

Share of the samples, self time (`deep` where marked: the function and everything it calls). `initial_state` is
`empty_active_set`, inlined; its time is all in `solve_lower_transposed` (99.6 to 99.99% of whose samples come from
`initial_state`). `add_constraint` is, by source line, `compute_directions` (inlined), apart from 0.1%.

| Function | RP tv 1000 | RP box 1000 | RP tr 1000 | MR lev 1000 | TE box 1000 | RP ≤ 100 | all ≤ 100 | n = 10 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `solve_lower_transposed` (base of the GI) | 30.9% | 20.3% | 28.6% | 50.6% | 5.6% | 22.1% | 26.5% | 6.8% |
| `gauss_newton` | 41.0% | 18.8% | 25.2% | - | - | 21.3% | 13.8% | - |
| `factor_in_place` (Cholesky) | 13.3% | 7.5% | 10.5% | 0.3% | 1.8% | 6.0% | 3.4% | 0.8% |
| `slack` (violation scan) | - | 20.9% | 8.4% | - | 17.1% | 9.2% | 9.3% | - |
| `add_constraint` (`compute_directions`) | 0.2% | 18.3% | 14.0% | 27.4% | 13.2% | 9.6% | 15.0% | 5.7% |
| `append_active` (Givens rotations) | - | 8.4% | 6.8% | 11.5% | 5.0% | 6.0% | 7.6% | - |
| `compensated_dot` (deep, with `fma`) | - | - | - | - | - | 10.3% | 6.3% | 12.1% |
| `reflect` (Householder) | - | 0.6% | - | - | 20.8% | 1.3% | 2.5% | 1.1% |
| `discard_inactive_generated` (deep) | - | - | - | 3.3% | - | - | - | - |
| operator new (`_Allocate`, deep) | - | - | - | 0.2% | - | - | - | 15.7% |

`-`: not among the 40 hottest functions of the session. `discard_inactive_generated` rebuilds the list of generated
rows through `append_generated` (98% of the samples of the latter); 72% of it is the copy of the rows
(`vector::insert`).

By source line (`--detail`): 98% of `solve_lower_transposed` is the strided `sum -= l[k * n + i] * b[k]`
(`cholesky.cpp:72`); 99% of `gauss_newton` the update `curvature[j * n + k] += scaled * jacobian[i * n + k]`
(`objective.cpp:197`); all of `slack` its `dot(row, x)` (`goldfarb_idnani.cpp:124`); 99% of `append_active` the
`rotate_basis` call (`active_set.hpp:183`), and 0.7% of it (0.06% of the samples) the `make_givens` line with its
`std::hypot`.

Callers of `gauss_newton` in RP tv 1000: `set_model` (the model curvature of every descent step) 62%, `linearize`
with the Hessian (`evaluate_kkt`, KKT checks and the polish) 27%, `descend` (`initial_proximal`) 11%. Callers of
`compensated_dot`: `covariance_times`, 100%; in the MSVC build 57% of its time (RP ≤ 100) is spent in the
`ucrtbase.dll` code it calls, its `fma`.

## Calls per solve

Breakpoint hits in one solve of the `gcc-bench` build, minus the calls of the harness (one `cholesky` in `Data::create`;
one `qp::solve` in `minimum_volatility` for `max_return` and `target_volatility`):

| Case | `qp::solve` | `gauss_newton` | `cholesky` | `householder_qr` | `covariance_times` |
| --- | ---: | ---: | ---: | ---: | ---: |
| `risk_parity/target_volatility/1000` | 29 | 46 | 41 | 7 | 156 |
| `risk_parity/box/1000` | 62 | 71 | 74 | 7 | 208 |
| `risk_parity/target_return/1000` | 67 | 76 | 79 | 7 | 223 |
| `max_return/leverage/1000` | 56 | 0 | 0 | 2 | 0 |
| `max_return/free/1000` | 4 | 0 | 0 | 1 | 0 |
| `max_return/long_only/1000`, `box`, `cash` | 11 | 0 | 0 (`cash`: 1) | 2 | 0 |
| `risk_parity/budget/1000` | 0 | 0 | 14 | 0 | 0 |
| `tracking_error/box/1000`, `leverage` | 1 | 0 | 0 | 1 | 0 |
| `max_sharpe/box`, `max_diversification/box`, `utility/leverage`, `min_variance/target_return` (1000) | 1 | 0 | 0 | 0 | 0 |
| `risk_parity/target_volatility/100` | 31 | 48 | 42 | 7 | 158 |
| `risk_parity/box/100` | 44 | 53 | 55 | 7 | 154 |
| `risk_parity/target_return/100` | 25 | 34 | 36 | 7 | 97 |
| `max_return/leverage/100` | 22 | 0 | 0 | 2 | 0 |

The `householder_qr` of `risk_parity` (the polish and the multipliers) are of the free-variable system, not 2n × n,
and are not priced below.

### Time of each kernel inside a case (`gcc-bench`, calls × seconds per call)

Seconds per call at n = 1000: `empty_active_set` 0.356 s, `gauss_newton` 0.191 s, `cholesky` 0.107 s,
`covariance_times` 3.94 ms, `qr` (2000 × 1000) 1.20 s; at n = 100: 405 µs, 152 µs, 69.6 µs, 39.4 µs.

| Case | Time | GI base (`empty_active_set`) | `gauss_newton` | `cholesky` | `covariance_times` | `qr` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| RP target vol 1000 | 26.4 s | 10.3 s (39%) | 8.8 s (33%) | 4.4 s (17%) | 0.61 s (2.3%) | - |
| RP box 1000 | 85.3 s | 22.1 s (26%) | 13.6 s (16%) | 7.9 s (9.3%) | 0.82 s (1.0%) | - |
| RP target return 1000 | 68.2 s | 23.9 s (35%) | 14.5 s (21%) | 8.5 s (12%) | 0.88 s (1.3%) | - |
| MR leverage 1000 | 32.3 s | 19.9 s (62%) | - | - | - | - |
| TE box 1000 | 4.06 s | 0.36 s (8.8%) | - | - | - | 1.20 s (30%) |
| RP target vol 100 | 34.1 ms | 12.6 ms (37%) | 7.3 ms (21%) | 2.9 ms (8.6%) | 6.2 ms (18%) | - |
| RP box 100 | 73.8 ms | 17.8 ms (24%) | 8.1 ms (11%) | 3.8 ms (5.2%) | 6.1 ms (8.2%) | - |
| RP target return 100 | 34.3 ms | 10.1 ms (30%) | 5.2 ms (15%) | 2.5 ms (7.3%) | 3.8 ms (11%) | - |
| MR leverage 100 | 16.6 ms | 8.9 ms (54%) | - | - | - | - |

What is left in RP box 1000 and RP target return 1000 (48% and 30%) matches the samples of the GI iterations in the
MSVC profile (`slack` + `compute_directions` + `append_active`: 48% and 29%). The GI base is a larger share in
`gcc-bench` than in the MSVC samples: the strided solve costs about the same in both builds (0.33 s a call in the MSVC
TE box session, one call), while the rest of the GCC code is faster.

## Assembly

`objdump -d -C -M intel --no-show-raw-insn` of the `portopt_core` objects of `gcc-bench`:
- `compensated_dot` (`objective.cpp.obj`): `call fma` once per term (relocation `IMAGE_REL_AMD64_REL32 fma`); in the
  executable `fma` is a `jmp [__imp_fma]`, the `fma` of `ucrtbase.dll`. The same holds for `std::fmax` (15 call sites
  in `polyhedral_descent.cpp.obj`) and `std::hypot` (2 in `goldfarb_idnani.cpp.obj`), which the profile does not show
  as hot.
- `solve_lower_transposed`: the inner loop loads `l[k * n + i]` and `l[(k + 1) * n + i]` with `movsd` / `movhpd`
  (stride n), one `mulpd`, then two `subsd` in order: a gather on a strided column, and a serial chain.
- Cholesky (`factor_in_place`, inlined in `cholesky`): `movupd` × 2, `mulpd`, then two dependent `subsd` per pair of
  terms. GCC vectorizes the multiplications and keeps the sum in order: the loop is bound by the latency of the
  subtraction chain.
- `gauss_newton`: the inner loop is `movupd`, `mulpd`, `addpd`, `movupd` (store): vectorized, but it reads and writes
  the n² curvature once per row i of the Jacobian, 8 MB n times at n = 1000.
- `slack`: the `dot` of a dense row of n doubles (`mulpd`, two serial `addsd`), also for a bound row that has one
  non-zero.

## Hypotheses

| Hypothesis | Verdict | Number |
| --- | --- | --- |
| H1. Base of the GI, n³/2 strided, in every `qp::solve` | **confirmed**, the largest | 26% to 62% of the slow n = 1000 cases (calls × 0.356 s); 20 to 51% of the MSVC samples; 90% of a one-iteration QP at n = 100 and 1000 (`min_variance/free` 387 ms against 356 ms) |
| H2. RISK_PARITY descent: `gauss_newton`, Cholesky and a QP per iteration | **confirmed** | `gauss_newton` 16 to 33% (46 to 76 calls, 0.191 s each), Cholesky 9 to 17% (41 to 79 calls), base of the QP as in H1; 62% of the `gauss_newton` samples come from the model curvature, 27% from the KKT checks |
| H3. Dot2 calls `fma` in `libm` in the generic build | **confirmed at n ≤ 100**, not at 1000 | `covariance_times` 8 to 18% of RP at n = 100 (97 to 158 calls, 39.4 µs each), 12% of the samples at n = 10; 1 to 2.3% at n = 1000; 3.95 ns per term, `call fma` in the object |
| H4. Dense bound rows scanned and projected per candidate | **confirmed where there are bounds** | `slack` 20.9% of RP box 1000, 8.4% of RP target return 1000, 17.1% of TE box 1000 (24% of its solve); plus part of `compute_directions` (18.3%, 14.0%, 13.2%); 0% with leverage only |
| H5. GI iteration linear algebra (`compute_directions`, Givens) | **confirmed** | `compute_directions` + `append_active` 39% of MR leverage 1000, 27% of RP box 1000, 21% of RP target return 1000; `std::hypot` under 1% |
| H6. Householder reflection with stride | **confirmed for TRACKING_ERROR only** | `reflect` 20.8% of the TE box 1000 samples (29% of its solve; one `qr` of 2000 × 1000, 1.20 s); under 2% in RISK_PARITY and MAX_RETURN |
| H7. Cholesky in order, a latency chain | **confirmed** (assembly), moderate share | 9 to 17% of RP at n = 1000, 5 to 9% at n = 100; `mulpd` with two serial `subsd` |
| H8. Allocations in the GI loop and the problems | **confirmed at n = 10 only** | operator new 15.7% of the samples at n = 10 (about 1 point of it the harness); under 2% at n ≥ 100; the rebuild of the generated rows (`discard_inactive_generated`, copies) 3.3% of MR leverage 1000; up to 1566 allocations per solve at n = 10, 183 690 in MR leverage 1000 |
| H9. The implicit family of 2^n rows is not a hot spot | **confirmed** (not hot) | `most_violated_inequality` 0.2% (deep, its `slack` aside) of MR leverage 1000 |
| H10. Minor: divisions in `model_curvature`, `jacobian_of`, `gradient_noise`, `volatility_of`, `std::hypot` | **refuted** (none hot) | `set_model` self (with the `k / n` loop) 0.25% of RP tv 1000; `jacobian_of` + `second_order_term` 1.1% of RP ≤ 100; `gradient_noise` 0.2%; `make_givens` under 1% |

## Levers for Tasks 2 to 4

Ordered by the estimated seconds each one saves on one pass of the n = 1000 cases (287.5 s of timed solves, which
dominate the suite). The factors are the exploratory experiments of the research (L1 8.4 times, L2 2.0, L3 4.0, Dekker
1.9 times) and the halving of `gauss_newton` by symmetry: **hypotheses to measure, not results**. Each lever is
measured A/B when it is made, and the profile is taken again after the first ones (the bottleneck moves).

The estimate is `time × savable fraction × (1 − 1 / factor)`, with the time from one of two bases:
- **(a)** gcc calls × kernel time ([Calls per solve](#calls-per-solve)), for the base, `gauss_newton`, `cholesky`,
  `qr` and `covariance_times`. Calls were counted on every case where they are not obvious: 10 more n = 1000 cases
  (`max_return` `free`, `long_only`, `box`, `cash`: 4, 11, 11 and 11 `qp::solve` calls; `risk_parity/budget`: 14
  `cholesky`, no `qp::solve` and no `gauss_newton`; one `qp::solve` in `max_sharpe/box`, `max_diversification/box`,
  `utility/leverage`, `min_variance/target_return` and `tracking_error/leverage`). The other QP objectives are taken as
  one `qp::solve` per case, like the five counted.
- **(b)** MSVC seconds of the function × (gcc solve time / MSVC solve time) of the case, for the GI iteration
  (`slack`, `compute_directions`), which has no kernel. **Rough.**

| Lever | Basis | Time it acts on (n = 1000) | Savable, factor | Saved (s) | n ≤ 100 |
| --- | --- | --- | --- | ---: | --- |
| L1 | (a) | base: 277 `qp::solve` × 0.356 s = 98.6 s (RP 158, MR 93, 26 QP cases × 1) | all, 8.4 | **86.9** | base 24 to 54% of the slow cases at n = 100 |
| L5 | (b), rough | RP box: `slack` 17.9 s + `compute_directions` 15.6 s; RP target return: 5.8 + 9.7 s; TE box: 0.96 + 0.74 s | `slack` 1, `compute_directions` 0.5 (its dot half), factor ≈ n | **37.6** profiled; **+24.4** extrapolated, very rough | `slack` 9.2% of RP ≤ 100 |
| L4 | (a) | 193 `gauss_newton` × 0.191 s = 36.9 s (RP only) | all, 2 | **18.4** | `gauss_newton` 11 to 21% of RP at n = 100 |
| L2 | (a) | 208 `cholesky` × 0.107 s = 22.3 s (RP, budget included) | all, 2.0 | **11.1** | Cholesky 5 to 9% of RP at n = 100 |
| L3 | (a) | 5 `qr` (2000 × 1000) × 1.20 s = 6.0 s (TE only) | all, 4.0 | **4.5** | `reflect` 2.5% of all ≤ 100 |
| L7 | (a), after L1, conditional | 88 repeated bases of the MR path × 0.356 s / 8.4 = 3.7 s | all, ∞ (not recomputed) | **3.7** (31.3 before L1) | MR leverage 100: 21 repeated bases |
| L6a | (a) | 587 `covariance_times` × 3.94 ms = 2.3 s (RP only) | all, 1.9 | **1.1** | `covariance_times` 8 to 18% of RP at n = 100; 12% of n = 10 |
| L8 | (b), rough | rebuild of the generated rows of MR leverage: 3.3% × 32.3 s × 0.85 = 0.9 s | at most all | **≤ 0.9** | operator new 15.7% of the n = 10 samples (0.84 ms of solves per pass) |

Notes on the rough entries:
- L5 profiled: `slack` over a box or long-only problem reads bound rows almost only, so all of it is taken as savable;
  `compute_directions` computes n dots and n − q axpys over the basis, and only the dots shrink for a unit row, so
  half of it is taken (between half and all, depending on q). The extrapolation applies the savable share of the GI
  iteration of TE box (64% of `slack` + `compute_directions` + `append_active` + the rest of
  `most_violated_inequality`) to the GI-iteration time of the bounded cases that were not profiled: `max_return`
  `long_only`, `box`, `cash` (case time − base: 25.5 s, so 16.3 s) and the `long_only`, `box`, `cash` cases of the
  five QP objectives (case time − base − `qr`: 12.6 s, so 8.1 s). Even without the extrapolation L5 ranks second.
- L7 assumes that every `qp::solve` of the MAX_RETURN path has the same Hessian factor, so the same J₀. Checked by
  reading the code, not numerically: `make_path` builds the one `Path` of a solve (`max_return.cpp:290`) with the
  factor of `mean_variance_problem`, and `solve_with`, `hold_face` and `release_face`, the only writers of the path's
  `qp::Problem`, change its linear term and its rows, never `hessian` or `hessian_factor`; in the MR leverage 1000
  samples 99.1% of the `qp::solve` time comes from `solve_with` (the rest is the harness). The gain is conditional on
  that, and on what L1 leaves.
- L8: the n = 10 gain is a share of a fraction of a millisecond per pass, nothing in seconds; it matters only for
  callers that solve many small problems.

Order: **L1, L5, L4, L2, L3, L7, L6a, L8**. L5 is the least certain of the large ones (basis (b)); L4 comes first for
`risk_parity/target_volatility`, which has no active bound rows (no `slack` in its samples).

Dropped:
- **L9, cache of `Moments` in the Armijo halvings.** `value_change` is 4.5% of the RP ≤ 100 samples and the cache
  would save at most the recomputation of one of its two `moments_of`, about 2%: below what the quick A/B resolves.
  Not visible at n = 1000.

Deferred, not measured here (Task 4 decides with its own profile):
- **L10, `-flto` on the core.** Needs `-ffp-contract=off` to hold at link time, not yet verified; the profile cannot
  price it (cross-TU calls are not separable in the samples). Decided in [profile 2](profile-2.md#l8-l9-and-l10): not adopted.
- **L6b, ISA dispatch (`target("avx2,fma")`).** Hardware FMA in the compensated dot after L6a; it acts on the same
  `covariance_times` time as L6a (2.3 s at n = 1000, 8 to 18% of RP at n = 100) and needs a manual dispatch on MinGW.

Not a lever of Tasks 2 to 4 but shown by the profile: after L1 and L5, the rest of the GI iteration (H5,
`compute_directions` on general rows and the Givens rotations of `append_active`: 39% of MR leverage 1000) is the next
hot spot of the QP; Task 4's profile decides what follows.
