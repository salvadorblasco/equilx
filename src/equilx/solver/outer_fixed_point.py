from equilx.data_structure import EquilibriumData, DebyeHuckelParameters
from equilx.mytypes import FArray
from equilx.solver.ionic import calc_ionic_strength
from equilx.solver.nr import newton_raphson


def outer_fixed_point(initial_guess, equil: EquilibriumData, dbparms: DebyeHuckelParameters, background_ionic: FArray):
    activity_changed: bool = True
    previous_concentration = initial_guess
    reference_logbeta = equil.log_beta.copy()
    while not activity_changed:
        concentration = newton_raphson(previous_concentration, equil)
        equil_ionic = calc_ionic_strength(concentration, equil)
        total_ionic = background_ionic + equil_ionic

