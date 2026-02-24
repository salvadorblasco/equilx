from typing import Final

from equilx.excepts import TooManyIterations
from equilx.data_structure import EquilibriumData, DebyeHuckelParameters
from equilx.mytypes import FArray
from equilx.solver.ionic import calc_ionic_strength, debye_huckel
from equilx.solver.nr import newton_raphson


def outer_fixed_point(initial_guess, equil: EquilibriumData, dbparms: DebyeHuckelParameters, background_ionic: FArray):
    MAX_ITERATIONS: Final[int] = 10
    activity_changed: bool = True
    previous_concentration = initial_guess
    ref_dh = debye_huckel(equil.reference_ionic, equil, dbparms)
    reference_logbeta = equil.log_beta[None,:] - ref_dh[:,None]

    for iteration in range(MAX_ITERATIONS):
        concentration = newton_raphson(previous_concentration, equil)
        equil_ionic = calc_ionic_strength(concentration, equil)
        total_ionic = background_ionic + equil_ionic
        new_dh = debye_huckel(total_ionic, equil, dbparms)
        new_logbeta = reference_logbeta + new_dh

        if np.max(np.abs(concentration - previous_concentration)) < THRESHOLD:
            break
    else:
        raise TooManyIterations("Outer fixed point reached iteration limit", concentration)
