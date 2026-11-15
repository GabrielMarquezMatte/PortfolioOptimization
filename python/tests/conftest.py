import pytest

try:
    import portopt  # noqa: F401
except ImportError as error:
    pytest.skip(f"PortOpt native library unavailable ({error}); set PORTOPT_LIB", allow_module_level=True)
