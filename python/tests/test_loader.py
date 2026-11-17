import ctypes
from pathlib import Path
from types import SimpleNamespace
from typing import cast

import pytest

from portopt import _native


def test_abi_version_matches() -> None:
    assert _native.lib.po_abi_version() == _native.PO_ABI_VERSION


def test_missing_library_message_names_env_var(tmp_path: Path) -> None:
    with pytest.raises(ImportError) as error:
        _native.locate_library(None, tmp_path)
    assert "PORTOPT_LIB" in str(error.value)
    assert str(tmp_path) in str(error.value)


def test_missing_symbol_raises_import_error_naming_it() -> None:
    missing = "po_portfolio_volatility"
    exported = {name: SimpleNamespace() for name in _native.PROTOTYPES if name != missing}
    stale = cast(ctypes.CDLL, SimpleNamespace(**exported))
    with pytest.raises(ImportError, match=missing) as error:
        _native.declare_prototypes(stale, Path("stale.dll"))
    assert "stale.dll" in str(error.value)


def test_out_of_memory_is_memory_error() -> None:
    with pytest.raises(MemoryError):
        _native.check(_native.PO_OUT_OF_MEMORY)


def test_version_is_the_installed_distribution_version() -> None:
    import importlib.metadata

    import portopt

    assert portopt.__version__ == importlib.metadata.version("pyportopt")
