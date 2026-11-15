# PortOpt Style Guide

Code style rules. Every rule is literal and checkable. If a rule blocks a clearly better solution, stop
and ask instead of breaking it silently.

Code, identifiers, error messages and repository documentation (README, STYLEGUIDE,
`docs/guide`) are in English. Working documents under `docs/superpowers/` are in Portuguese.

## Priorities

1. **Correctness** — the code must be right. No shortcuts that sacrifice safety or accuracy.
2. **Performance and readability, together** — code is read far more than it is written, and this library exists
   to be fast. The two are not traded against each other by default:
   - A change that is faster and reads as well (an elided copy, an allocation removed, a move instead of a copy) is
     always taken, anywhere in the code.
   - A change that is faster but reads worse is taken only on a measured hot path, with the measurement recorded in
     `docs/` (see [Performance Notes](#performance-notes)).
   - A change that reads better but is slower on a measured hot path is not taken.

## Automatic checks

`tests/style/check_style.py` (standard library only) runs in `ctest` (test `style`) and in CI:

```
python tests/style/check_style.py [ROOT]
```

It prints one `path:line: rule: message` line per violation and exits 1 if there is any. It scans
`src/`, `include/`, `cpp/`, `tests/`, `python/src/`, `python/tests/` and `python/benchmarks/`, and skips `build/`
and `python/.venv/`.

| Rule             | C / C++                | Python                              |
|------------------|------------------------|-------------------------------------|
| `else`           | token outside strings and comments | `else`/`elif`/`for…else`/`try…else` by AST |
| `ternary`        | review only            | by AST: a ternary whose branches are not a name, literal, attribute or negated literal, that is nested, or that spans more than one line |
| `folder-modules` | more than 8 modules    | more than 8 `.py` files (not counting `__init__.py`) |
| `function-length`| clang-tidy             | more than 60 lines, by AST          |
| `function-params`| clang-tidy             | more than 5 parameters, by AST      |
| `nesting`        | clang-tidy             | more than 3 nested blocks, by AST   |

In C++, function length, parameter count and nesting are enforced by `clang-tidy`
(`readability-function-size`: `LineThreshold: 60`, `ParameterThreshold: 5`, `NestingThreshold: 4`,
configured in `.clang-tidy`) in CI and locally when available. clang-tidy counts the function body as level 1, so
its threshold of 4 is the limit of 3 below. Formatting comes from `.clang-format`.

Copy elision on return ([Return values](#return-values--rvo-and-nrvo)) is enforced by the compiler: GCC 14+ builds
with `-Wnrvo` as an error, and GCC and Clang reject pessimizing and redundant moves (`-Wpessimizing-move`,
`-Wredundant-move`, part of `-Wall`/`-Wextra`).

---

## Never `else`

Forbidden: `else`, `else if`; in Python also `elif`, `for…else`, `while…else`, `try…else`. Use guard
clauses and early returns. To choose among several cases, use a sequence of `if` with `return`, a
`switch`, or a table. (The rule does not apply to CMake files.)

```cpp
// NOT OK
double sign_of(double x)
{
    if (x > 0.0)
    {
        return 1.0;
    }
    else if (x < 0.0)
    {
        return -1.0;
    }
    else
    {
        return 0.0;
    }
}

// OK
double sign_of(double x)
{
    if (x > 0.0)
    {
        return 1.0;
    }
    if (x < 0.0)
    {
        return -1.0;
    }
    return 0.0;
}
```

---

## Control Flow — Prefer Early Return

Avoid wrapping the happy path in a condition. Guard and exit early; keep the main logic at the lowest
indentation level.

```cpp
// NOT OK: the body is nested under the guard
void append_decoded(std::span<const std::byte> source)
{
    if (!source.empty())
    {
        reserve(length_ + source.size());
        length_ += decode(source, buffer_.data() + length_);
    }
}

// OK: the guard exits early
void append_decoded(std::span<const std::byte> source)
{
    if (source.empty())
    {
        return;
    }
    reserve(length_ + source.size());
    length_ += decode(source, buffer_.data() + length_);
}
```

### `if` always has braces, body on its own line

Never `if (cond) return;` on one line, never a brace-less body.

```cpp
// NOT OK
if (series.empty()) return 0.0;
if (series.empty())
    return 0.0;

// OK
if (series.empty())
{
    return 0.0;
}
```

---

## Ternaries

Use a ternary only with one condition and operands that are names, literals or members. Never nested,
never with a function call in a branch, never broken across lines.

`check_style.py` enforces this for Python (rule `ternary`: each branch must be a name, a literal, an
attribute or a negated literal such as `-1.0`; no nested ternary; one source line). In C++ ternaries are
checked by review only.

```cpp
// OK: two literals
double sign = inverted ? -1.0 : 1.0;

// OK: one member on each side
std::size_t start = options.has_start ? options.start : 0;

// NOT OK: nested
double scale = a ? (b ? 1.0 : 2.0) : 3.0;

// NOT OK: a call in a branch
double value = positive ? std::sqrt(x) : 0.0;

// Use if instead:
double value = 0.0;
if (positive)
{
    value = std::sqrt(x);
}
```

When in doubt, use `if`.

---

## Scope Nesting — Maximum 3 Levels

Count every block that introduces indentation: `if`, `for`, `while`, `switch`, `try`, `with`. A
function body is level 0. Nesting beyond 3 is a signal to extract a function.

```
level 0: function body
level 1: for (asset : universe)
level 2:   if (asset.valid())
level 3:     for (value : asset)       <- maximum allowed
level 4:       if (std::isnan(value))  <- extract this into a function
```

When you reach level 4, ask what the inner block is doing. Name it and move it.

---

## Function Length — Maximum 60 Lines

Count every line from the signature to the closing brace, excluding blank lines and documentation
comments. If a function exceeds 60 lines, split it. The target is 30.

Do not game the limit by packing logic onto one line. The limit exists to keep functions focused, not
to reward line-golf. A function that packs three operations on one line to stay under 60 is worse than
one that spreads them across 61 readable lines and is then split.

```cpp
// NOT OK: three statements on one line to save two lines
double a = f(x); double b = g(a); return h(a, b);

// OK: split by responsibility
double variance = portfolio_variance(weights, covariance);
double risk_penalty = risk_aversion * variance;
return expected_return - risk_penalty;
```

Common splits: validate vs. compute; compute one term vs. combine terms; find a boundary vs. process
what is inside it.

---

## Parameters — Maximum 5 Per Function

This includes the C interface. More than 5 means a handle with setters (C) or an aggregate type
(C++/Python).

```cpp
// NOT OK
po_status po_spec_configure(po_spec* spec, int32_t objective, int32_t long_only, double risk_aversion,
                            double target_return, int32_t max_assets);

// OK (C interface): one setter per option on an opaque handle
po_status po_spec_set_objective(po_spec* spec, int32_t objective);
po_status po_spec_set_long_only(po_spec* spec, int32_t long_only);

// OK (C++): an aggregate for internal code
struct SolveConfig
{
    bool long_only;
    double risk_aversion;
    std::size_t max_iterations;
};
Result solve(const Data& data, const Spec& spec, const SolveConfig& config);
```

---

## Modules and Folders — Maximum 8 Per Folder

A module is a `.hpp`/`.cpp` pair with the same name, or a lone file. In Python, a `.py` file not
counting `__init__.py`. At the limit, create a subfolder per context. A file above ~300 lines is a
sign to split it.

```
// NOT OK: src/problems/ with 9 modules
src/problems/{spec,translate,solve,constraints,objective,bounds,report,export,cache}.{hpp,cpp}

// OK: split by context
src/problems/{spec,translate,solve}.{hpp,cpp}
src/problems/constraints/{budget,bounds,target}.{hpp,cpp}
```

No `utils`, `helpers`, `misc`, `common`: the module name states its context.

```cpp
// NOT OK
#include "core/utils.hpp"

// OK
#include "core/data.hpp"
```

No speculative abstraction: no base class with one implementation, no factory for one product, no
configuration for a value that never changes.

---

## Layers

Dependencies point one way: `core` (spec, data, result, failure) <- `problems` (translates spec and data to a QP and
solves it) <- `capi`. `linalg` (Cholesky) depends on nothing; `qp` (Goldfarb-Idnani) depends only on `linalg` and is
used by `problems`. `metrics` (portfolio return and volatility) depends on `core` and is used by `capi`. Nothing in
`src/` includes from `src/capi/` except `capi` itself.

```cpp
// NOT OK: src/problems/solve.cpp
#include "capi/error.hpp"

// OK: src/capi/solve_api.cpp
#include "problems/solve.hpp"
```

---

## Naming

- Types: `PascalCase` (`SolveResult`, `DataView`)
- Functions and variables: `snake_case` (`solve_problem`, `variance`)
- Private members: trailing underscore (`length_`)
- Constants: `kPascalCase` (`kMaxLags`)
- C interface: prefix `po_` for functions and types, `PO_` for macros and constants
- Whole words (`index`, not `idx`; `source`, not `src`) unless the abbreviation is idiomatic
  (`len`, `buf`, `pos`) or mathematical notation (`mu`, `lambda`, `sigma`)
- Braces in Allman style, 4 spaces of indentation (`.clang-format`)

```cpp
// NOT OK
struct solve_result { double* Weights; int nIdx; };

// OK
struct SolveResult
{
    std::vector<double> weights;
    std::size_t asset_count;
};
```

### Public API naming must be symmetric across problem kinds

What exists for one objective has the same shape for the others. A caller who learns one entry point should be able
to guess the next: `po_spec_set_risk_aversion` and `po_spec_set_target_return`, not one of those plus a bare
`po_spec_return_target`. When you add a member for one kind, find its siblings and match them, or add the missing
siblings in the same change.

---

## Comments

**The default is no comment.** Code explains itself through naming and flow. Reaching for a comment is
the signal to look at the code again: if a reader cannot follow the logic, that is a defect in the
code, not a gap in its documentation. Rename the variable, extract the function, flatten the branch.

A comment is justified only when the reason lives outside the code and no naming can bring it in:

- A constraint imposed by an external convention (a paper's definition), a platform or a format.
- An invariant a reader could silently break, where breaking it still compiles and passes review.

```cpp
// OK: an external convention the code cannot express
// Goldfarb and Idnani (1983) require a positive definite Hessian.

// NOT OK: describes what the code does
// Loop over each asset and solve it.

// NOT OK: the code should have said this itself
// p is the position of the first non-finite value
std::size_t p = first_non_finite(returns);   // name it first_non_finite_position and delete the comment
```

Never reference a plan, spec, task or ticket in a comment.

### Findings do not belong in code

Measurements, benchmark results and the reasoning that picked one approach over another are **not
comments**. They go in `docs/`, with the numbers, the machine and the data, so the next person can
re-run them and disagree. What may stay in the code is the constraint itself, without its evidence:

```cpp
// OK: the ceiling, not the benchmark that found it
// At most 4096 assets per problem; larger universes are rejected.

// NOT OK: a finding, and it rots the moment anything around it changes
// Measured 11 ms vs 36 ms for the scalar loop on 5000 assets, so the blocked version stays.
```

### A comment must stay true, or it is worse than no comment

- Do not describe a state the code is no longer in. When you refactor, reread the comments on what you
  touched.
- Do not cite a file that can disappear (plans, specs, section markers). A committed document under
  `docs/` is different in kind: it is permanent, and pointing at one is how a finding stays out of the
  code. The test is whether the target is expected to still exist in a year.

---

## Numerics

- Only `double`.
- Never form an explicit inverse to solve a system: use a Cholesky factorization or a solve.
- Log-determinants come from the Cholesky factor, never `log(det(...))`.
- Keep covariance matrices symmetrized (`0.5 * (P + P^T)`).
- Break long formulas into named terms that match the documented equation.
- Parameter constraints are imposed by a documented reparametrization (`exp` for variances, logistic
  for values in (0, 1)), never by clamping inside the objective.
- A missing value (NaN) inside the returns is an argument error, never silently filled.

```cpp
// NOT OK: one line, no names
double u = w.dot(m) - 0.5 * a * w.dot(S * w);

// OK: terms named after the documented equation
double expected_return = weights.dot(mean_returns);
double variance = weights.dot(covariance * weights);
double utility = expected_return - 0.5 * risk_aversion * variance;
```

---

## Error Handling

Report errors at trust boundaries: every function of the C interface, constructor arguments, and
malformed external data. Do not add defensive checks for conditions that are impossible given the
internal invariants.

- Every C function that can fail returns `po_status`; the detailed message is in `po_last_error()`,
  which is thread-local.
- No exception crosses the C interface. Every `extern "C"` function catches everything:
  `std::bad_alloc` becomes `PO_OUT_OF_MEMORY`, anything else becomes `PO_INTERNAL_ERROR`.
- Match the status to whose fault it is. A bad argument is `PO_INVALID_ARGUMENT`; an infeasible
  problem is a status on the result, never an error of the call.
- Internal code (below `capi`) may throw; the `capi` guard translates at the boundary.

```cpp
// NOT OK: an exception escapes through the C interface
extern "C" po_status po_solve(po_problem* problem, po_result** result)
{
    return run(problem, result);   // may throw
}

// OK: the boundary guards everything
extern "C" po_status po_solve(po_problem* problem, po_result** result)
{
    return guard([&] { return run(problem, result); });
}
```

---

## Untrusted Input

Every argument that crosses the C interface is untrusted: callers are Python, C, or anything that can
load the library. These rules are not optional in `capi` code.

- **Validate every length, index, count and stride before it drives an allocation, an index or a loop
  bound.** Check against what the caller's buffer can actually hold (the declared `T`, `N`, strides and
  ranges) and reject with `PO_INVALID_ARGUMENT` before allocating, not after.
- **When you add a guard for one argument, apply it to every sibling argument.** The expensive bugs are
  inconsistency, not ignorance: a validated asset count next to an unvalidated weight offset. If one
  field needs the guard, assume its neighbours do too.
- **Check size arithmetic for overflow.** `count * stride * sizeof(double)` computed in `size_t` can wrap
  into a small number that sails past a bounds check. Use a checked multiply/add helper and reject on
  overflow.
- **Never let a malformed value silently change semantics.** A negative stride, a null pointer with a
  non-zero length, or a length that wraps must be reported, not reinterpreted.
- **Validate before you allocate, not before you use.** Allocating first and discovering the
  inconsistency during the walk still hands the caller the allocation.
- New or changed argument handling needs a matching test that forges the malformed call. Assert the
  status code and the message, not merely that something failed.

```cpp
// NOT OK: the product can overflow before it is compared
if (count * stride > capacity)
{
    return PO_INVALID_ARGUMENT;
}

// OK: checked arithmetic, rejected before allocating
std::size_t required = 0;
if (!checked_mul(count, stride, required) || required > capacity)
{
    return set_error(PO_INVALID_ARGUMENT, "asset range exceeds the data");
}
```

---

## Return values — RVO and NRVO

Return results by value and write the function so the compiler builds the result directly in the caller's storage:
no copy and no move. This applies everywhere, not only on hot paths: it costs no readability.

- **Return by value**, not through an output parameter. An output parameter (`std::span<double>`, `T&`) is for a
  buffer the caller reuses across calls on a hot path, or for filling an object the caller itself returns by name.
- **Every `return` in a function returns the same named local, or every `return` is a prvalue.** GCC elides the
  named return (NRVO) only then: a single guard clause that returns some other object disables it for the whole
  function, even a guard placed before the local is declared.
- **Guards that return something else live in a caller.** Split the function: the outer one keeps its guard
  clauses and returns prvalues (`return unavailable(k);`, `return from_factor(factor, middle, k);`); the inner one
  builds the named result and returns only it.
- **To return "nothing" from a named `std::optional`, reset it and return it** (`hessian.reset(); return
  hessian;`), never `return std::nullopt;` beside `return hessian;`.
- **The local has the exact return type.** A `Matrix` returned from a function that returns `std::optional<Matrix>`
  is moved into the optional, not elided: declare `std::optional<Matrix> hessian(std::in_place, m * m);` instead.
- **Never `return std::move(local);`**: it forbids the elision it was meant to help. The same for
  `return std::move(parameter);`, which is already an implicit move.
- **Do not build a fallback result up front** (`SolveResult failed = ...;` returned by early exits): it allocates on
  the paths that never use it and blocks NRVO on the one that returns the real result.

```cpp
// NOT OK: the guard returns another object, so `errors` is moved, not elided
StandardErrors from_hessian(const Matrix& hessian, const Matrix& middle, std::size_t k)
{
    Matrix factor = negated(hessian);
    if (!cholesky_factor(factor, k))
    {
        return unavailable(k);
    }
    StandardErrors errors = ...;
    return errors;
}

// OK: the guard stays in the outer function; the inner one returns only its named local
StandardErrors from_hessian(const Matrix& hessian, const Matrix& middle, std::size_t k)
{
    Matrix factor = negated(hessian);
    if (!cholesky_factor(factor, k))
    {
        return unavailable(k);
    }
    return from_factor(factor, middle, k);
}
```

---

## Ownership and RAII

Resources are owned by objects, released by destructors, on every exit path including exceptions.

- No naked `new`/`delete` and no manual `free`: `std::vector`, `std::unique_ptr` or a small RAII type.
- A buffer that lives beyond one stack frame has one owner, stated by its type.
- Bound reads by the logical length you tracked, never by the capacity of the container.
- Handles that cross the C interface are opaque; freeing one never invalidates another. Results are
  immutable once produced, so concurrent readers need no lock.

```cpp
// NOT OK: leaks if solve() throws
double* work = new double[length];
solve(work, length);
delete[] work;

// OK
std::vector<double> work(length);
solve(work);
```

---

## Tests and Benchmarks

Harness code is held to the same standard as library code, because a broken harness reports success.

- **Fail loudly, never silently.** A fixture that was never built, a setup that never ran, an input
  that parsed to nothing: assert it. The dangerous harness is not the one that crashes, it is the one
  that measures zero work and publishes a number.
- Assert real behavior: exact values, status codes, messages. A test that calls a function and asserts
  "did not throw" documents nothing.
- Name tests for what they cover (`Abi.RuntimeVersionMatchesHeader`). A hand-written return series is not a
  "real world" fixture, and calling it one hides the gap it was meant to close.
- Python tests are `test_*` functions with plain `assert`; use `pytest.mark.parametrize` instead of a
  loop and no `unittest` classes.
- Benchmarks against another library compare matched work, or say plainly that they do not. State the
  machine for every published number, and never derive a ratio across different hardware.
- Tests never touch the network.

```cpp
// NOT OK
TEST(Solve, Works)
{
    EXPECT_NO_THROW(solve(problem));
}

// OK
TEST(Solve, InfeasibleBoundsReportInfeasible)
{
    EXPECT_EQ(solve_status(impossible_bounds), PO_INFEASIBLE);
}
```

---

## Performance Notes

The goal is the fastest code the problem allows, proved by evidence and not by intuition. How a hot-path change is
measured (benchmark, profile, assembly) is recorded in `docs/performance/`.

These rules apply only when a function is on a measured hot path (the QP active-set loop, the Cholesky
factorization):

- Prefer `std::span<const double>` over `const std::vector<double>&` or a pointer and length pair for
  data the function only reads.
- Allocate once, outside the loop: workspaces are sized up front and reused across observations and
  iterations.
- Keep results bit-identical for any number of threads: no reduction whose order depends on scheduling.
- Split cold paths (error construction, growth) out of the hot function so they do not bloat it.
- Optimize against a measurement, and record what it showed in `docs/`, never in a comment. An
  optimization with no number behind it is a guess that costs readability.

```cpp
// NOT OK: allocates on every observation
for (std::size_t t = 0; t < length; ++t)
{
    std::vector<double> scratch(parameter_count);
    accumulate(scratch, t);
}

// OK
std::vector<double> scratch(parameter_count);
for (std::size_t t = 0; t < length; ++t)
{
    accumulate(scratch, t);
}
```

Outside the hot path (setup, argument validation, result assembly), ordinary allocations are fine.

---

## Checklist before finishing

- [ ] `python tests/style/check_style.py` exits 0 (also run by `ctest`).
- [ ] No `else` in C/C++/Python; every `if` has braces and its body on its own line.
- [ ] No function above 60 lines, with more than 5 parameters, or nested deeper than 3 blocks.
- [ ] No folder with more than 8 modules; no `utils`/`helpers`/`misc`/`common`.
- [ ] No comment that describes what the code does, cites a plan or task, or records a measurement.
- [ ] Results are returned by value with NRVO or RVO: the GCC build (`-Wnrvo`) has no warning, and no
      `return std::move(local)`.
- [ ] A hot-path change comes with before/after benchmark numbers, and its assembly was read; findings are in
      `docs/performance/`.
- [ ] Every new C function validates its arguments, returns `po_status` and lets no exception escape.
- [ ] Tests assert exact values and pass.
