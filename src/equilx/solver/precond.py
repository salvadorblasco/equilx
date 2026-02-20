from typing import Final

import numpy as np

from equilx.mytypes import FArray
from ..data_structure import EquilibriumData


def precondition(equil: EquilibriumData, initial_guess: FArray | None = None) -> FArray:
    EXIT_THRESHOLD: Final[float] = 1e-9
    def _exit_condition(sumr: float, sump: float) -> bool:
        return abs(sumr-sump)/(sumr+sump) < EXIT_THRESHOLD


def _sum_react(concentration: FArray, equil: EquilibriumData) -> float:
    ...


def _pick_exponent(equil: EquilibriumData) -> float:
    stoich = equil.stoichiometryx
    return np.argmax(stoich, axis=1)
