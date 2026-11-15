import pytest

collect_ignore_glob: list[str] = []
UNAVAILABLE = ""

try:
    import portopt  # noqa: F401
except ImportError as error:
    UNAVAILABLE = f"PortOpt native library unavailable, all tests skipped: {error}"
    collect_ignore_glob = ["test_*.py"]


def pytest_terminal_summary(terminalreporter: pytest.TerminalReporter) -> None:
    if UNAVAILABLE:
        terminalreporter.write_line(UNAVAILABLE, yellow=True)
