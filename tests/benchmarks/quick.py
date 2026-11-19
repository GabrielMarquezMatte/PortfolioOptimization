"""Quick benchmark suite and A/B comparisons. Standard library only. See docs/performance/measuring.md.

    python tests/benchmarks/quick.py run --out FILE [--exe PATH]
    python tests/benchmarks/quick.py compare A.json B.json
    python tests/benchmarks/quick.py compare-solutions A.json B.json [--require-identical]
    python tests/benchmarks/quick.py selftest

`compare` reads two reports of `run`; `compare-solutions` reads two dumps of `portopt_profile --json`.
"""

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

DEFAULT_EXE = Path("build/gcc-bench/portopt_benchmarks" + (".exe" if sys.platform == "win32" else ""))
QUICK_FILTER = (
    "^BM_Solve/(min_variance/long_only|utility/box|max_return/long_only|max_sharpe/long_only"
    "|tracking_error/long_only|max_diversification/box|risk_parity/(budget|target_volatility|box))/(10|100)$"
    "|^BM_Kernel/(cholesky|solve_lower|qr|qp_iteration|compensated_dot|covariance_times|gauss_newton"
    "|empty_active_set)/100$"
)
SOLUTIONS_FORMAT = "portopt-solutions-1"
IDENTITY_KEYS = ("status", "objective", "weights", "multipliers", "bound_multipliers", "residual")


def medians(report: dict[str, Any]) -> dict[str, tuple[float, float]]:
    """run_name -> (median, standard deviation) of real_time, from a Google Benchmark JSON with aggregates."""
    values: dict[str, dict[str, float]] = {}
    for entry in report["benchmarks"]:
        if entry.get("run_type") != "aggregate":
            continue
        values.setdefault(entry["run_name"], {})[entry["aggregate_name"]] = entry["real_time"]
    return {name: (kinds["median"], kinds.get("stddev", 0.0)) for name, kinds in values.items() if "median" in kinds}


def compare_rows(a: dict[str, Any], b: dict[str, Any]) -> list[tuple[str, float, float, float, bool]]:
    """(name, median A, median B, B / A, within noise) for every benchmark in both runs."""
    before = medians(a)
    after = medians(b)
    rows = []
    for name, (median_a, sd_a) in before.items():
        if name not in after:
            continue
        median_b, sd_b = after[name]
        rows.append((name, median_a, median_b, median_b / median_a, abs(median_b - median_a) <= max(sd_a, sd_b)))
    return rows


def solutions_by_name(dumped: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """name -> solution record of a dump of `portopt_profile --json`; a repeated name is an error."""
    if dumped.get("format") != SOLUTIONS_FORMAT:
        raise ValueError(f"not a solutions dump: format is {dumped.get('format')!r}, expected {SOLUTIONS_FORMAT!r}")
    indexed: dict[str, dict[str, Any]] = {}
    for record in dumped["solutions"]:
        if record["name"] in indexed:
            raise ValueError(f"solution {record['name']} appears twice")
        indexed[record["name"]] = record
    return indexed


def hex_difference(a: str, b: str) -> float:
    """|a - b| of two hex-float texts: 0 for equal values (two NaNs included), infinity for NaN against a number."""
    x = float.fromhex(a)
    y = float.fromhex(b)
    if x == y or (math.isnan(x) and math.isnan(y)):
        return 0.0
    if math.isnan(x) or math.isnan(y):
        return math.inf
    return abs(x - y)


def weight_difference(a: list[str], b: list[str]) -> float:
    """max|a - b| / max|a| of two weight vectors given as hex-float texts; infinity when the sizes differ."""
    if len(a) != len(b):
        return math.inf
    magnitudes = [abs(float.fromhex(text)) for text in a]
    scale = max((value for value in magnitudes if math.isfinite(value)), default=0.0)
    worst = max((hex_difference(x, y) for x, y in zip(a, b)), default=0.0)
    if scale == 0.0:
        return worst
    return worst / scale


def objective_difference(a: str, b: str) -> float:
    """|a - b| / |a| of two objective values given as hex-float texts (absolute when a is 0 or not finite)."""
    difference = hex_difference(a, b)
    scale = abs(float.fromhex(a))
    if scale == 0.0 or not math.isfinite(scale):
        return difference
    return difference / scale


def is_identical(a: dict[str, Any], b: dict[str, Any]) -> bool:
    return all(a[key] == b[key] for key in IDENTITY_KEYS)


def worse(current: tuple[float, str], value: float, name: str) -> tuple[float, str]:
    if value > current[0]:
        return (value, name)
    return current


def iteration_ratio(before: dict[str, dict[str, Any]], after: dict[str, dict[str, Any]], names: list[str]) -> float:
    iterations_a = sum(int(before[name]["iterations"]) for name in names)
    iterations_b = sum(int(after[name]["iterations"]) for name in names)
    if iterations_a == 0:
        return 1.0
    return iterations_b / iterations_a


def compare_solutions(a: dict[str, Any], b: dict[str, Any]) -> dict[str, Any]:
    """Compares two dumps of solutions by case name; the order of the cases does not matter.

    identical counts the cases whose status, objective, weights and multipliers are equal as hex-float text,
    that is, bit for bit. The worst differences are relative."""
    before = solutions_by_name(a)
    after = solutions_by_name(b)
    common = sorted(set(before) & set(after))
    worst_weight = (0.0, "")
    worst_objective = (0.0, "")
    not_identical = []
    for name in common:
        left, right = before[name], after[name]
        if not is_identical(left, right):
            not_identical.append(name)
        worst_weight = worse(worst_weight, weight_difference(left["weights"], right["weights"]), name)
        worst_objective = worse(worst_objective, objective_difference(left["objective"], right["objective"]), name)
    return {
        "compared": len(common),
        "identical": len(common) - len(not_identical),
        "not_identical": not_identical,
        "status_changed": [name for name in common if before[name]["status"] != after[name]["status"]],
        "iterations_changed": [name for name in common if before[name]["iterations"] != after[name]["iterations"]],
        "iteration_ratio": iteration_ratio(before, after, common),
        "worst_weight_difference": worst_weight[0],
        "worst_weight_difference_case": worst_weight[1],
        "worst_objective_difference": worst_objective[0],
        "worst_objective_difference_case": worst_objective[1],
        "only_in_a": sorted(set(before) - set(after)),
        "only_in_b": sorted(set(after) - set(before)),
    }


def load(path: str) -> Any:
    with open(path, encoding="utf-8") as file:
        return json.load(file)


def run(out: str, exe: str) -> None:
    executable = Path(exe).resolve()
    if not executable.exists():
        raise SystemExit(f"{executable} does not exist; build it with: cmake --build --preset gcc-bench")
    command = [str(executable), f"--benchmark_filter={QUICK_FILTER}", "--benchmark_min_time=0.3s", "--benchmark_repetitions=3",
               "--benchmark_report_aggregates_only=true", f"--benchmark_out={out}", "--benchmark_out_format=json"]
    started = time.perf_counter()
    subprocess.run(command, check=True)
    print(f"quick suite: {time.perf_counter() - started:.1f} s")


def print_compare(a_path: str, b_path: str) -> None:
    print(f"{'benchmark':<60} {'A':>12} {'B':>12} {'B/A':>7}")
    for name, median_a, median_b, ratio, noise in compare_rows(load(a_path), load(b_path)):
        print(f"{name:<60} {median_a:>12.4g} {median_b:>12.4g} {ratio:>7.3f} {'noise' if noise else ''}")


def print_names(label: str, names: list[str]) -> None:
    if names:
        print(f"{label}: {', '.join(names)}")


def print_solutions(summary: dict[str, Any]) -> None:
    print(f"cases compared: {summary['compared']}")
    print(f"bit-identical: {summary['identical']} of {summary['compared']}")
    print(f"worst relative weight difference: {summary['worst_weight_difference']:.3g}"
          f" ({summary['worst_weight_difference_case'] or '-'})")
    print(f"worst relative objective difference: {summary['worst_objective_difference']:.3g}"
          f" ({summary['worst_objective_difference_case'] or '-'})")
    print(f"iteration ratio B/A: {summary['iteration_ratio']:.4f}")
    print_names("not bit-identical", summary["not_identical"])
    print_names("status changed", summary["status_changed"])
    print_names("iterations changed", summary["iterations_changed"])
    print_names("only in A", summary["only_in_a"])
    print_names("only in B", summary["only_in_b"])


def run_compare_solutions(a_path: str, b_path: str, require_identical: bool) -> int:
    summary = compare_solutions(load(a_path), load(b_path))
    print_solutions(summary)
    complete = not summary["only_in_a"] and not summary["only_in_b"]
    if require_identical and not (complete and summary["identical"] == summary["compared"]):
        return 1
    return 0


def aggregate(name: str, kind: str, value: float) -> dict[str, Any]:
    return {"run_name": name, "run_type": "aggregate", "aggregate_name": kind, "real_time": value}


def selftest_compare() -> None:
    a = {"benchmarks": [aggregate("x", "median", 10.0), aggregate("x", "stddev", 1.0), aggregate("y", "median", 4.0),
                        aggregate("y", "stddev", 0.1), aggregate("only_a", "median", 1.0)]}
    b = {"benchmarks": [aggregate("x", "median", 10.5), aggregate("x", "stddev", 0.2), aggregate("y", "median", 2.0),
                        aggregate("y", "stddev", 0.1)]}
    assert medians(a)["x"] == (10.0, 1.0)
    rows = {row[0]: row for row in compare_rows(a, b)}
    assert set(rows) == {"x", "y"}
    assert rows["x"][4] is True
    assert rows["y"][3] == 0.5
    assert rows["y"][4] is False


def solution(name: str, weights: list[float], iterations: int = 5, status: int = 0) -> dict[str, Any]:
    return {"name": name, "status": status, "iterations": iterations, "objective": (0.5).hex(),
            "weights": [weight.hex() for weight in weights], "multipliers": {"budget": (0.25).hex()},
            "bound_multipliers": [(0.0).hex()], "residual": (0.0).hex()}


def dump(*solutions: dict[str, Any]) -> dict[str, Any]:
    return {"format": SOLUTIONS_FORMAT, "solutions": list(solutions)}


def selftest_solutions_identical() -> None:
    a = dump(solution("x/10", [0.25, 0.75]), solution("y/10", [1.0, 0.0]))
    summary = compare_solutions(a, a)
    assert summary["compared"] == 2
    assert summary["identical"] == 2
    assert summary["worst_weight_difference"] == 0.0


def selftest_solutions_one_ulp() -> None:
    weight = 0.3
    neighbour = math.nextafter(weight, 1.0)
    a = dump(solution("x/10", [weight, 0.7]), solution("y/10", [1.0, 0.0]))
    b = dump(solution("x/10", [neighbour, 0.7]), solution("y/10", [1.0, 0.0]))
    summary = compare_solutions(a, b)
    assert summary["identical"] == 1
    assert summary["not_identical"] == ["x/10"]
    assert summary["worst_weight_difference_case"] == "x/10"
    assert 0.0 < summary["worst_weight_difference"] <= 2.0**-52


def selftest_solutions_swapped_order() -> None:
    a = dump(solution("x/10", [0.25, 0.75]), solution("y/10", [1.0, 0.0]), solution("z/10", [0.5, 0.5]))
    b = dump(*reversed(a["solutions"]))
    summary = compare_solutions(a, b)
    assert summary["compared"] == 3
    assert summary["identical"] == 3


def selftest_solutions_differences() -> None:
    a = dump(solution("x/10", [0.5, 0.5], iterations=4), solution("only_a/10", [1.0]))
    b = dump(solution("x/10", [0.5, 0.5], iterations=6, status=3), solution("only_b/10", [1.0]))
    summary = compare_solutions(a, b)
    assert summary["compared"] == 1
    assert summary["identical"] == 0
    assert summary["status_changed"] == ["x/10"]
    assert summary["iterations_changed"] == ["x/10"]
    assert summary["iteration_ratio"] == 1.5
    assert summary["only_in_a"] == ["only_a/10"]
    assert summary["only_in_b"] == ["only_b/10"]


def selftest_solutions_not_a_dump() -> None:
    try:
        compare_solutions({"benchmarks": []}, dump())
    except ValueError as error:
        assert "not a solutions dump" in str(error)
        return
    raise AssertionError("a Google Benchmark report was accepted as a solutions dump")


def selftest() -> None:
    selftest_compare()
    selftest_solutions_identical()
    selftest_solutions_one_ulp()
    selftest_solutions_swapped_order()
    selftest_solutions_differences()
    selftest_solutions_not_a_dump()
    print("selftest ok")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("selftest")
    run_parser = commands.add_parser("run")
    run_parser.add_argument("--out", required=True)
    run_parser.add_argument("--exe", default=str(DEFAULT_EXE))
    compare_parser = commands.add_parser("compare")
    solutions_parser = commands.add_parser("compare-solutions")
    solutions_parser.add_argument("--require-identical", action="store_true")
    for pair in (compare_parser, solutions_parser):
        pair.add_argument("a")
        pair.add_argument("b")
    return parser


def main() -> int:
    arguments = build_parser().parse_args()
    if arguments.command == "selftest":
        selftest()
    if arguments.command == "run":
        run(arguments.out, arguments.exe)
    if arguments.command == "compare":
        print_compare(arguments.a, arguments.b)
    if arguments.command == "compare-solutions":
        return run_compare_solutions(arguments.a, arguments.b, arguments.require_identical)
    return 0


if __name__ == "__main__":
    sys.exit(main())
