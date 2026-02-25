from typing import Final

import numpy as np

from equilx.excepts import TooManyIterations
from equilx.data_structure import EquilibriumData, DebyeHuckelParameters
from equilx.mytypes import FArray
from equilx.solver.ionic import calc_ionic_strength, debye_huckel
from equilx.solver.nr import newton_raphson
from equilx.solver.utils import species_concentration


def outer_fixed_point(initial_guess, equil: EquilibriumData, dbparms: DebyeHuckelParameters) -> FArray:
    MAX_ITERATIONS: Final[int] = 10
    THRESHOLD: Final[float] = 1e-9
    previous_concentration = initial_guess
    ref_dh = debye_huckel(np.full(equil.n_points, equil.reference_ionic), equil, dbparms)
    reference_logbeta = equil.log_beta[None,:] - ref_dh
    concentration = np.empty_like(initial_guess)

    for iteration in range(MAX_ITERATIONS):
        concentration[...] = newton_raphson(previous_concentration, equil)
        full_conc = species_concentration(concentration, equil, full=True)
        equil_ionic = calc_ionic_strength(full_conc, equil)
        total_ionic = equil.background_ions + equil_ionic
        new_dh = debye_huckel(total_ionic, equil, dbparms)
        new_logbeta = reference_logbeta + new_dh
        equil.log_beta = new_logbeta

        # test =np.max(np.abs(concentration - previous_concentration))
        # print(iteration, test, np.sum(total_ionic), np.sum(equil_ionic))
        if np.max(np.abs(concentration - previous_concentration)) < THRESHOLD:
            break
        previous_concentration[...] = concentration
    else:
        raise TooManyIterations("Outer fixed point reached iteration limit", concentration)

    return concentration
