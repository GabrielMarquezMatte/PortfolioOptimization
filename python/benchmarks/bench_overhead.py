"""Fixed cost of the Python layer against the native solve: ctypes calls, handles, numpy, the readers.

    PORTOPT_LIB=<path of libportopt.dll> python python/benchmarks/bench_overhead.py [--sizes 10,100,1000]
                                                                                    [--repeats 9]

The instance is the one of tests/benchmarks (a covariance of condition 1e2, means in [0.02, 0.12]) drawn with numpy,
so its numbers are not the bits of the C++ one; the case is MIN_VARIANCE long only. Every row is the median of
`repeats` samples in microseconds; a sample of a fast call is the mean of `inner` calls. See
docs/performance/python-vs-arch.md.
"""

from __future__ import annotations

import argparse
import ctypes
import math
import statistics
import time
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

import portopt as po
from portopt import _native, _result
from portopt._native import lib
from portopt._spec import doubles_pointer

CONDITION = 100.0
SEED = 1
SLOW_SECONDS = 0.05


@dataclass(frozen=True)
class Row:
    name: str
    call: Callable[[], object]
    inner: int = 1


def instance(n: int) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    generator = np.random.default_rng(SEED)
    q, _ = np.linalg.qr(generator.standard_normal((n, n)))
    eigenvalues = 0.04 * np.exp(np.linspace(0.0, 1.0, n) * math.log(CONDITION))
    covariance = (q * eigenvalues) @ q.T
    mu = generator.uniform(0.02, 0.12, n)
    return mu, 0.5 * (covariance + covariance.T)


def sample_ns(call: Callable[[], object], inner: int) -> float:
    started = time.perf_counter_ns()
    for _ in range(inner):
        call()
    return (time.perf_counter_ns() - started) / inner


def median_us(row: Row, repeats: int) -> float:
    row.call()
    samples = [sample_ns(row.call, row.inner) for _ in range(repeats)]
    return statistics.median(samples) / 1e3


def batch_of(n: int) -> int:
    """Calls per sample of a row whose time grows with n: enough that the 100 ns clock of Windows resolves it."""
    return max(1, 1000 // n)


def native_rows(spec: po.Spec, data: po.Data, mu: NDArray[np.float64], cov: NDArray[np.float64]) -> list[Row]:
    batch = batch_of(mu.size)
    count_out = ctypes.c_size_t()
    count_ref = ctypes.byref(count_out)
    spec_pointer, data_pointer = spec.pointer, data.pointer
    mu_pointer, cov_pointer = doubles_pointer(mu), doubles_pointer(cov)

    def solve_raw() -> None:
        out = ctypes.c_void_p()
        _native.check(lib.po_solve(spec_pointer, data_pointer, ctypes.byref(out)))
        lib.po_result_free(out)

    def data_create_raw() -> None:
        _native.create_handle(lib.po_data_create, lib.po_data_free, mu.size, mu_pointer, cov_pointer).close()

    return [
        Row("no-op Python call", lambda: None, 5000),
        Row("ctypes call, no arguments (po_abi_version)", lib.po_abi_version, 5000),
        Row("ctypes call through the guard (po_data_asset_count)",
            lambda: lib.po_data_asset_count(data_pointer, count_ref), 5000),
        Row("same, through read_scalar", lambda: _native.read_scalar(ctypes.c_size_t, lib.po_data_asset_count,
                                                                     data_pointer), 2000),
        Row("po_data_create + free through ctypes", data_create_raw, batch),
        Row("po_solve + po_result_free through ctypes", solve_raw, batch),
    ]


def python_rows(spec: po.Spec, data: po.Data, mu: NDArray[np.float64], cov: NDArray[np.float64]) -> list[Row]:
    result = _native.create_handle(lib.po_solve, lib.po_result_free, spec.pointer, data.pointer)
    batch = batch_of(mu.size)
    count = data.weight_count
    return [
        Row("Spec(...) and its free", lambda: po.Spec(po.Objective.MIN_VARIANCE, bounds=(0.0, math.inf)), 500),
        Row("Data(mu, cov) and its free", lambda: po.Data(mu, cov), batch),
        Row("read_result (10 ctypes calls, 1 buffer)", lambda: _result.read_result(result.pointer, count), 200),
        Row("solve(spec, data)", lambda: po.solve(spec, data), batch),
    ]


def estimation_rows(data: po.Data, cov: NDArray[np.float64]) -> list[Row]:
    """The equilibrium from the covariance (a Data built, validated and factored per call) and from a prebuilt Data."""
    batch = batch_of(cov.shape[0])
    w_mkt = np.full(cov.shape[0], 1.0 / cov.shape[0])
    return [
        Row("bl_equilibrium(cov, w_mkt, delta)", lambda: po.bl_equilibrium(cov, w_mkt, 2.5), batch),
        Row("bl_equilibrium(data, w_mkt, delta)", lambda: po.bl_equilibrium(data, w_mkt, 2.5), batch),
    ]


def check_case(spec: po.Spec, data: po.Data) -> None:
    result = po.solve(spec, data)
    if result.status != po.SolveStatus.OPTIMAL:
        raise RuntimeError(f"the case ended {result.status.name}, not OPTIMAL")


def measure(n: int, repeats: int) -> dict[str, float]:
    mu, cov = instance(n)
    spec = po.Spec(po.Objective.MIN_VARIANCE, bounds=(0.0, math.inf))
    data = po.Data(mu, cov)
    check_case(spec, data)
    rows = native_rows(spec, data, mu, cov) + python_rows(spec, data, mu, cov) + estimation_rows(data, cov)
    return {row.name: median_us(row, repeats) for row in rows}


def table(sizes: list[int], results: dict[int, dict[str, float]]) -> str:
    header = "| Measure (microseconds, median) | " + " | ".join(f"n = {n}" for n in sizes) + " |"
    divider = "| --- | " + " | ".join("---:" for _ in sizes) + " |"
    lines = [header, divider]
    for name in results[sizes[0]]:
        cells = " | ".join(f"{results[n][name]:.3f}" for n in sizes)
        lines.append(f"| {name} | {cells} |")
    return "\n".join(lines)


def parse() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sizes", default="10,100,1000")
    parser.add_argument("--repeats", type=int, default=9)
    return parser.parse_args()


def main() -> None:
    arguments = parse()
    sizes = [int(text) for text in arguments.sizes.split(",")]
    results = {n: measure(n, arguments.repeats) for n in sizes}
    print(table(sizes, results))


if __name__ == "__main__":
    main()
