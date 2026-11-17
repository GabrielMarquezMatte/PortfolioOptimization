"""Portfolio optimization over the PortOpt C API."""

from importlib.metadata import version

from ._data import Data
from ._native import PortOptError
from ._result import Result, SolveStatus, solve
from ._spec import Objective, Spec

# The version of the installed distribution, the one pip shows.
__version__ = version("pyportopt")

__all__ = ["Data", "Objective", "PortOptError", "Result", "SolveStatus", "Spec", "solve"]
