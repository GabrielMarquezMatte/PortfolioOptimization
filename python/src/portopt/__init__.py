"""Portfolio optimization over the PortOpt C API."""

from importlib.metadata import version

from ._data import Data
from ._estimation import (
    BayesStein,
    BlackLittermanPosterior,
    Views,
    bayes_stein,
    bl_equilibrium,
    bl_omega_he_litterman,
    bl_omega_idzorek,
    bl_posterior,
    sample_moments,
)
from ._native import PortOptError
from ._result import Result, SolveStatus, solve
from ._spec import Objective, Spec

# The version of the installed distribution, the one pip shows.
__version__ = version("pyportopt")

__all__ = [
    "BayesStein",
    "BlackLittermanPosterior",
    "Data",
    "Objective",
    "PortOptError",
    "Result",
    "SolveStatus",
    "Spec",
    "Views",
    "bayes_stein",
    "bl_equilibrium",
    "bl_omega_he_litterman",
    "bl_omega_idzorek",
    "bl_posterior",
    "sample_moments",
    "solve",
]
