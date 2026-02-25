from functools import partial
from typing import Final

import numpy as np
from numpy.typing import NDArray

from equilx.mytypes import FArray, IArray
from ..data_structure import EquilibriumData
from .utils import species_concentration, mass_action


def precondition(equil: EquilibriumData, initial_guess: FArray | None = None) -> FArray:
    EXIT_THRESHOLD: Final[float] = 0.1
    MAX_ITERATIONS: Final[int] = 50

    # TODO remove
    if initial_guess is None:
        initial_guess = 0.01*equil.analytical_concentration
        initial_guess[initial_guess<0.0] = 1e-10

    stoich = _calc_reactive_stoich(equil)
    exponent = _pick_exponent(equil)
    neg_analc = equil.analytical_concentration < 0.0
    theta = np.zeros((equil.n_points, equil.n_components), dtype=float)
    positive_stoich = stoich.copy()
    positive_stoich[positive_stoich<0] = 0
    negative_stoich = -stoich     # negatives -> positives
    negative_stoich[negative_stoich<0] = 0      # cancel former negatives
    assert np.all(negative_stoich>=0)
    abs_analc = np.abs(equil.analytical_concentration)
    calc_ratio = partial(_calc_ratio, positive_stoich= positive_stoich, negative_stoich=negative_stoich,
                         abs_analc=abs_analc, neg_analc=neg_analc)
    best: float = 1e99

    def _exit_condition(sumr: FArray, sump: FArray, cext: FArray) -> bool:
        nonlocal best
        sumsq = float(np.sum(mass_action(equil, cext)**2))
        # print('\t', np.max(abs(sumr-sump)/(sumr+sump)), sumsq)
        if sumsq > best:
            return True
        else:
            best = sumsq

        thres_cond = bool(np.all(abs(sumr-sump)/(sumr+sump) < EXIT_THRESHOLD))
        return thres_cond

    x = initial_guess.copy()
    for iteration in range(MAX_ITERATIONS):
        extconc = species_concentration(x, equil, full=True)
        sump, sumr = calc_ratio(extconc)
        ratio = sump/sumr
        # print(iteration, np.min(ratio), np.max(ratio))
        _cap_ratio(ratio)
        x = _newx(x, ratio, theta, exponent)
        theta = _update_theta(ratio, theta)
        if _exit_condition(sumr, sump, extconc):
            break

    return x


def _cap_ratio(ratio: FArray) -> None:
    CAPPING = 1.0e4
    ratio[ratio>CAPPING] = CAPPING
    ratio[ratio<(1/CAPPING)] = 1/CAPPING


def _calc_ratio(conc: FArray, positive_stoich: IArray, negative_stoich: IArray, abs_analc: FArray,
                neg_analc: NDArray[bool]) -> tuple[FArray, FArray]:
    sum_pr = conc @ positive_stoich
    sum_nr = abs_analc + sum_pr
    sum_np = conc @ negative_stoich
    sum_pp = abs_analc + sum_np
    sum_p = np.where(neg_analc, sum_np, sum_pp)
    sum_r = np.where(neg_analc, sum_nr, sum_pr)
    sum_p[sum_p<=0.0] = 1e-10
    sum_r[sum_r<=0.0] = 1e-10
    return sum_p, sum_r


def _newx(x0, ratio, theta, expn):
    return x0*(theta*ratio**expn[None,:] + (1-theta))


def _update_theta(ratio: FArray, prev_theta: FArray) -> FArray:
    ALPHA: Final[float] = 0.9
    BETA: Final[float] = 0.8
    aux1 = ALPHA - BETA*ratio
    aux2 = ALPHA - BETA/ratio
    new_theta = np.where(ratio>1.0, aux2, aux1)
    return np.where(new_theta>prev_theta, new_theta, prev_theta)


def _pick_exponent(equil: EquilibriumData) -> FArray:
    stoich = equil.stoichiometryx
    aux = np.max(stoich, axis=0)
    return 1/aux


def _calc_reactive_stoich(equil: EquilibriumData) -> IArray:
    return np.concatenate((-np.eye(equil.n_components, dtype=int),
                           equil.stoichiometry),
                          axis=0)
