from pathlib import Path

import pytest

from portopt import _native


def test_abi_version_matches() -> None:
    assert _native.lib.po_abi_version() == _native.PO_ABI_VERSION


def test_missing_library_message_names_env_var(tmp_path: Path) -> None:
    with pytest.raises(ImportError) as error:
        _native.locate_library(None, tmp_path)
    assert "PORTOPT_LIB" in str(error.value)
    assert str(tmp_path) in str(error.value)
