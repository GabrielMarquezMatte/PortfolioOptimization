"""ctypes bindings that mirror include/portopt/portopt.h."""

from __future__ import annotations

import ctypes
import os
from collections.abc import Callable, Mapping
from pathlib import Path
from types import TracebackType
from typing import Any, cast

PO_ABI_VERSION = 1
PO_OK = 0
PO_INVALID_ARGUMENT = 1
PO_UNSUPPORTED = 3
PO_OUT_OF_MEMORY = 4
PO_INTERNAL_ERROR = 5
LIBRARY_ENV = "PORTOPT_LIB"
LIBRARY_DIRECTORY = Path(__file__).resolve().parent / "_lib"
LIBRARY_NAMES = ("libportopt.dll", "portopt.dll", "libportopt.so", "libportopt.dylib")

_HANDLE = ctypes.c_void_p
_OUT_HANDLE = ctypes.POINTER(ctypes.c_void_p)
_STATUS = ctypes.c_int32
_SIZE = ctypes.c_size_t
_I32 = ctypes.c_int32
_DOUBLE = ctypes.c_double
_IN_DOUBLES = ctypes.POINTER(ctypes.c_double)
_OUT_DOUBLE = ctypes.POINTER(ctypes.c_double)
_OUT_SIZE = ctypes.POINTER(ctypes.c_size_t)
_OUT_I32 = ctypes.POINTER(ctypes.c_int32)

PROTOTYPES: dict[str, tuple[Any, list[Any]]] = {
    "po_abi_version": (_I32, []),
    "po_last_error": (ctypes.c_char_p, []),
    "po_spec_create": (_STATUS, [_I32, _OUT_HANDLE]),
    "po_spec_set_bounds": (_STATUS, [_HANDLE, _DOUBLE, _DOUBLE]),
    "po_spec_set_weight_bounds": (_STATUS, [_HANDLE, _IN_DOUBLES, _IN_DOUBLES, _SIZE]),
    "po_spec_set_target_return": (_STATUS, [_HANDLE, _DOUBLE]),
    "po_spec_set_risk_aversion": (_STATUS, [_HANDLE, _DOUBLE]),
    "po_spec_set_importance": (_STATUS, [_HANDLE, _DOUBLE]),
    "po_spec_set_max_volatility": (_STATUS, [_HANDLE, _DOUBLE]),
    "po_spec_set_leverage": (_STATUS, [_HANDLE, _DOUBLE]),
    "po_spec_set_target_volatility": (_STATUS, [_HANDLE, _DOUBLE]),
    "po_spec_set_risk_budget": (_STATUS, [_HANDLE, _IN_DOUBLES, _SIZE]),
    "po_spec_free": (None, [_HANDLE]),
    "po_data_create": (_STATUS, [_SIZE, _IN_DOUBLES, _IN_DOUBLES, _OUT_HANDLE]),
    "po_data_set_risk_free": (_STATUS, [_HANDLE, _DOUBLE]),
    "po_data_set_returns": (_STATUS, [_HANDLE, _IN_DOUBLES, _IN_DOUBLES, _SIZE]),
    "po_data_asset_count": (_STATUS, [_HANDLE, _OUT_SIZE]),
    "po_data_weight_count": (_STATUS, [_HANDLE, _OUT_SIZE]),
    "po_data_free": (None, [_HANDLE]),
    "po_solve": (_STATUS, [_HANDLE, _HANDLE, _OUT_HANDLE]),
    "po_result_status": (_STATUS, [_HANDLE, _OUT_I32]),
    "po_result_weight_count": (_STATUS, [_HANDLE, _OUT_SIZE]),
    "po_result_weights": (_STATUS, [_HANDLE, _OUT_DOUBLE, _SIZE]),
    "po_result_iterations": (_STATUS, [_HANDLE, _OUT_I32]),
    "po_result_multiplier": (_STATUS, [_HANDLE, _I32, _OUT_DOUBLE]),
    "po_result_bound_multipliers": (_STATUS, [_HANDLE, _OUT_DOUBLE, _SIZE]),
    "po_result_residual": (_STATUS, [_HANDLE, _OUT_DOUBLE]),
    "po_result_free": (None, [_HANDLE]),
    "po_portfolio_return": (_STATUS, [_HANDLE, _IN_DOUBLES, _SIZE, _OUT_DOUBLE]),
    "po_portfolio_volatility": (_STATUS, [_HANDLE, _IN_DOUBLES, _SIZE, _OUT_DOUBLE]),
    "po_utility": (_STATUS, [_HANDLE, _IN_DOUBLES, _SIZE, _DOUBLE, _OUT_DOUBLE]),
    "po_tracking_error": (_STATUS, [_HANDLE, _IN_DOUBLES, _SIZE, _OUT_DOUBLE]),
    "po_sharpe_ratio": (_STATUS, [_HANDLE, _IN_DOUBLES, _SIZE, _OUT_DOUBLE]),
    "po_diversification_ratio": (_STATUS, [_HANDLE, _IN_DOUBLES, _SIZE, _OUT_DOUBLE]),
    "po_risk_contributions": (_STATUS, [_HANDLE, _IN_DOUBLES, _SIZE, _OUT_DOUBLE]),
    "po_sample_moments": (_STATUS, [_SIZE, _SIZE, _IN_DOUBLES, _OUT_DOUBLE, _OUT_DOUBLE]),
    "po_bayes_stein": (_STATUS, [_HANDLE, _SIZE, _OUT_DOUBLE, _OUT_DOUBLE, _OUT_DOUBLE]),
    "po_bl_equilibrium": (_STATUS, [_HANDLE, _IN_DOUBLES, _SIZE, _DOUBLE, _OUT_DOUBLE]),
    "po_views_create": (_STATUS, [_SIZE, _SIZE, _IN_DOUBLES, _IN_DOUBLES, _OUT_HANDLE]),
    "po_views_set_omega": (_STATUS, [_HANDLE, _IN_DOUBLES, _SIZE]),
    "po_views_free": (None, [_HANDLE]),
    "po_bl_omega_he_litterman": (_STATUS, [_HANDLE, _HANDLE, _DOUBLE, _OUT_DOUBLE]),
    "po_bl_omega_idzorek": (_STATUS, [_HANDLE, _HANDLE, _DOUBLE, _IN_DOUBLES, _OUT_DOUBLE]),
    "po_bl_posterior": (_STATUS, [_HANDLE, _HANDLE, _DOUBLE, _OUT_DOUBLE, _OUT_DOUBLE]),
}


class PortOptError(Exception):
    """A native call failed; `status` is the PO_* code and the message comes from po_last_error()."""

    def __init__(self, status: int, message: str) -> None:
        super().__init__(message)
        self.status = status


def locate_library(override: str | None, directory: Path) -> Path:
    if override:
        path = Path(override)
        if not path.is_file():
            raise ImportError(f"{LIBRARY_ENV} points to {override}, which is not a file")
        # Absolute: the loader looks a bare file name up on its search path, not in the working directory.
        return path.resolve()
    for name in LIBRARY_NAMES:
        if (directory / name).is_file():
            return directory / name
    raise ImportError(
        f"PortOpt native library not found: set {LIBRARY_ENV} to the library file, "
        f"or place one of {', '.join(LIBRARY_NAMES)} in {directory}"
    )


def load_library(environ: Mapping[str, str], directory: Path) -> ctypes.CDLL:
    path = locate_library(environ.get(LIBRARY_ENV), directory)
    try:
        library = ctypes.CDLL(str(path))
    except OSError as error:
        raise ImportError(f"cannot load the PortOpt native library {path}: {error}") from error
    check_abi(library, path)
    declare_prototypes(library, path)
    return library


def check_abi(library: ctypes.CDLL, path: Path) -> None:
    try:
        declare_function(library, "po_abi_version")
    except AttributeError as error:
        raise ImportError(f"{path} does not export po_abi_version, it is not a PortOpt library") from error
    abi = library.po_abi_version()
    if abi != PO_ABI_VERSION:
        raise ImportError(f"{path} has ABI version {abi}, this package needs {PO_ABI_VERSION}")


def declare_function(library: ctypes.CDLL, name: str) -> None:
    restype, argtypes = PROTOTYPES[name]
    function = getattr(library, name)
    function.restype = restype
    function.argtypes = argtypes


def declare_prototypes(library: ctypes.CDLL, path: Path) -> None:
    for name in PROTOTYPES:
        try:
            declare_function(library, name)
        except AttributeError as error:
            raise ImportError(f"{path} does not export {name}: it is probably an older PortOpt library") from error


lib = load_library(os.environ, LIBRARY_DIRECTORY)


def last_error() -> str:
    message = cast(bytes | None, lib.po_last_error())
    if message is None:
        return ""
    return message.decode("utf-8", errors="replace")


def check(status: int) -> None:
    if status == PO_OK:
        return
    if status == PO_INVALID_ARGUMENT:
        raise ValueError(last_error())
    if status == PO_UNSUPPORTED:
        raise NotImplementedError(last_error())
    if status == PO_OUT_OF_MEMORY:
        raise MemoryError(last_error())
    raise PortOptError(status, last_error())


class Handle:
    """Owns one native handle and frees it exactly once."""

    def __init__(self, pointer: int | None, free: Callable[[int], None]) -> None:
        self._pointer: int | None = pointer
        self._free = free

    @property
    def pointer(self) -> int:
        if self._pointer is None:
            raise ValueError("the native handle is already closed")
        return self._pointer

    def close(self) -> None:
        pointer, self._pointer = self._pointer, None
        if pointer is not None:
            self._free(pointer)

    def __enter__(self) -> Handle:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()

    def __del__(self) -> None:
        self.close()


def create_handle(create: Callable[..., int], free: Callable[[int], None], *args: Any) -> Handle:
    out = ctypes.c_void_p()
    check(create(*args, ctypes.byref(out)))
    if out.value is None:
        raise PortOptError(PO_INTERNAL_ERROR, "the native create call returned a null handle")
    return Handle(out.value, free)


def read_scalar(ctype: Any, function: Callable[..., int], *args: Any) -> Any:
    out = ctype()
    check(function(*args, ctypes.byref(out)))
    return out.value
