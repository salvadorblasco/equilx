from functools import partial
from typing import Final, Any, Optional

import numpy as np
from numpy.typing import NDArray

from equilx.mytypes import FArray, IArray, BArray
from equilx.data_structure import EquilibriumData
from equilx.solver.utils import calc_saturation, make_titration, mass_action, species_concentration


def precondition(equil: EquilibriumData, initial_guess: FArray, active_solid_mask: Optional[BArray]=None) -> tuple[FArray, Optional[FArray]]:
    EXIT_THRESHOLD: Final[float] = 1e-9
    MAX_ITERATIONS: Final[int] = 500

    do_solids = (active_solid_mask is not None)
    if active_solid_mask is None:
        active_solid_mask = np.ones((equil.n_points, equil.n_solids), dtype=bool)

    solid_controler = _choose_precipitator(equil)
    exponent = _pick_exponent(equil)

    neg_analc = equil.analytical_concentration < 0.0
    abs_analc = np.abs(equil.analytical_concentration)
    all_stoich = (np.vstack((equil.stoichiometryx, equil.solid_stoichiometry))
                  if do_solids else equil.stoichiometryx) 
    positive_stoich = all_stoich.clip(min=0)
    negative_stoich = (-all_stoich).clip(min=0)     # negatives -> positives
    calc_ratio = partial(_calc_ratio, positive_stoich= positive_stoich, negative_stoich=negative_stoich,
                         abs_analc=abs_analc, neg_analc=neg_analc)

    theta = 0.01+np.zeros((equil.n_points, equil.n_components), dtype=float)

    def _exit_condition(sumr: FArray, sump: FArray) -> bool:
        thres_cond = bool(np.all(abs(sumr-sump)/(sumr+sump) < EXIT_THRESHOLD))
        return thres_cond

    x: FArray = initial_guess
    xs: FArray = np.zeros((equil.n_points, equil.n_solids), dtype=float)

    max_j = _estimate_cai(equil)
    
    for _ in range(MAX_ITERATIONS):
        aqueous = species_concentration(x, equil, full=True)

        if do_solids:
            xs = _update_fixed_solids(aqueous, xs, equil, active_solid_mask, solid_controler)

        all_conc = np.hstack((aqueous, xs)) if do_solids else aqueous
        sump, sumr = calc_ratio(all_conc)

        ratio = sump / sumr

        # Dynamic exponent — always 2D
        expn = np.tile(exponent, (equil.n_points, 1))

        if do_solids:
            si = calc_saturation(x, equil)                    # compute once

            for nsolid, j_ctrl in enumerate(solid_controler):
                ap_j = equil.solid_stoichiometry[nsolid, j_ctrl]
                if ap_j <= 0:
                    continue

                active = active_solid_mask[:, nsolid] & (si[:, nsolid] > 1.0)

                if np.any(active):
                    sum_reac = si[active, nsolid] ** (1.0 / ap_j)
                    ratio[active, j_ctrl] = 1.0 / np.maximum(sum_reac, 1e-12)
                    expn[active, j_ctrl] = 1.0 / ap_j

        x = _newx(x, ratio, theta, expn)

        x = x.clip(min=0.0)
        x = np.where(x > max_j, np.sqrt(max_j * x), x)
        if do_solids:
            xs = xs.clip(min=0.0)

        theta = _update_theta(ratio)
        if _exit_condition(sumr, sump):
            break

    return x, (xs if do_solids else None)


def _estimate_cai(equil: EquilibriumData, hindex: Optional[int]=None):
    SCALE: Final[float] = 100.0
    FLOOR: Final[float] = 1e-12
    maxj = SCALE*np.maximum(np.abs(equil.analytical_concentration), FLOOR)
    if hindex is not None:
        maxj[:, hindex] = np.minimum(maxj[:, hindex], 1.0)
    return maxj


def _update_fixed_solids(aqueous_conc: FArray, solid_conc: FArray, equil: EquilibriumData, active_mask, index_control: list[int]) -> FArray:
    new_xs = solid_conc.copy()
    aux1 = aqueous_conc @ equil.stoichiometryx
    aux2 = solid_conc @ equil.solid_stoichiometry

    for nsolid, jcontrol in enumerate(index_control):
        expn = int(equil.solid_stoichiometry[nsolid, jcontrol])
        sum_other_solids = aux2[:, jcontrol] - solid_conc[:, nsolid]*expn
        anc_j = equil.analytical_concentration[:,jcontrol]
        new_xs[:,nsolid] = (anc_j-aux1[:,jcontrol]-sum_other_solids)/expn

    new_xs[np.logical_not(active_mask)] = 0.0
    return new_xs.clip(min=0.0)


def _calc_ratio(conc: FArray, positive_stoich: IArray, negative_stoich: IArray, abs_analc: FArray,
                neg_analc: NDArray[np.bool_]) -> tuple[FArray, FArray]:
    FLOOR: Final[float] = 1e-10
    sum_pr = conc @ positive_stoich
    sum_nr = abs_analc + sum_pr
    sum_np = conc @ negative_stoich
    sum_pp = abs_analc + sum_np
    sum_p = np.where(neg_analc, sum_np, sum_pp)
    sum_r = np.where(neg_analc, sum_nr, sum_pr)
    return sum_p.clip(min=FLOOR), sum_r.clip(min=FLOOR)


def _newx(x0, ratio, theta, expn):
    return x0*(theta*ratio**expn + (1-theta))


def _update_theta(ratio: FArray) -> FArray:
    ALPHA: Final[float] = 0.9
    BETA: Final[float] = 0.8
    aux1 = ALPHA - BETA*ratio
    aux2 = ALPHA - BETA/ratio
    new_theta = np.where(ratio>1.0, aux2, aux1)
    return new_theta


def _pick_exponent(equil: EquilibriumData) -> FArray:
    stoich = equil.stoichiometryx
    aux = np.max(stoich, axis=0)
    return 1/aux


def _choose_precipitator(eq: EquilibriumData) -> list[int]:
    ss = eq.solid_stoichiometry.copy()
    ss[ss < 1] = 0      # remove negatives
    retval: list[Any] = eq.n_solids * [None]
    for i, j in zip(*np.nonzero(ss)):
        if retval[i] is not None:
            continue
        retval[i] = int(j)
    return retval


if __name__ == "__main__":
    analc = make_titration(initial_conc=np.array([0.05, 0.1, 0.2]),
                           buret_conc=np.array([0.0, 0.0, -1.0]),
                           starting_vol=10.0,
                           titre=np.linspace(0.0, 3.8, 20))
    eq = EquilibriumData(analytical_concentration=analc,
                         stoichiometry=np.array([[1,0,1],
                                            [1,0,2],
                                            [1,1,0],
                                            [0,0,-1]], dtype=int),
                         log_beta=np.array([16.33, 22.68, -12.06, -14.0]),
                         solid_stoichiometry=np.array([[0,1,-2],
                                                       [1,1,0]], dtype=int),
                         log_ks=np.array([23.0, -8.4]))
    xL = np.full((eq.n_points, eq.n_components), 1e-6)
    xS = np.zeros((eq.n_points, eq.n_solids), dtype=float)
    controller = _choose_precipitator(eq)

    c = species_concentration(xL, eq, full=True)
    y0 = mass_action(eq, c, xS)

    xL, _ = precondition(eq, xL)
    c = species_concentration(xL, eq, full=True)
    y1 = mass_action(eq, c, xS)

    si = calc_saturation(xL, eq)
    solid_mask = si > 1.0

    xL, xS = precondition(eq, xL, active_solid_mask=solid_mask)
    c = species_concentration(xL, eq, full=True)
    y2 = mass_action(eq, c, xS)

    print(np.sum(y0**2))
    print(np.sum(y1**2))
    print(np.sum(y2**2))
    print(xL)
    print(xS)

