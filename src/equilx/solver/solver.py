"Module for solving equilibrium problems."

from typing import Optional

import numpy as np

from equilx.mytypes import FArray
from equilx.data_structure import EquilibriumData, DebyeHuckelParameters

from equilx.solver import fsolver


def Eq(equil: EquilibriumData,
       initial_guess: Optional[FArray],
       debye_params: Optional[DebyeHuckelParameters],
       freeze: Optional[int], **kwargs) -> dict:
    "Main function for equilibrium solving."

    retdict = {}
    control = 0

    # Extract parameters
    if not initial_guess:
        strategy = kwargs.get('initial_guess_strategy', 4)
        initial_guess = guess_starting_point(equil, strategy)

    free_conc = np.asfortranarray(initial_guess)
    anal_conc = np.asfortranarray(equil.analytical_concentration)
    beta = np.asfortranarray(10**equil.log_beta)
    stoichiometry = np.asfortranarray(equil.stoichiometry)

    with_solids = equil.any_solid()
    if with_solids:
        solid_stoichiometry = np.asfortranarray(equil.solid_stoichiometry)
        ks = np.asfortranarray(10**equil.log_ks)
        solid_conc = np.zeros((equil.n_points, equil.n_solids), order='F')

    variable_ionic = (debye_params is not None)
    if variable_ionic:
        ...

    # option 1: simple titration, no solids, no variable ionic strength
    if not with_solids and not variable_ionic:
        # call solve(C, T, B, P, CTRL, NPT, NCP, NEQ)
        fsolver.solve(free_conc, anal_conc, beta, stoichiometry, control)
    elif with_solids and not variable_ionic:
        # call solve_solids
        fsolver.solve_solids(free_conc, solid_conc, anal_conc, beta, stoichiometry,
                             ks, solid_stoichiometry, control)
    elif not with_solids and variable_ionic:
        # call solve_vis
        fsolver.solve_vis(free_conc, anal_conc, beta, stoichiometry, control)
    elif with_solids and variable_ionic:
        # call solve_solids_vis
        fsolver.solve_solids_vis(free_conc, anal_conc, beta, stoichiometry, control)

    return retdict


def guess_starting_point(equil: EquilibriumData, strategy: int=4) -> FArray:
    """Compute initial guess for concentration based on strategy.

    equil: EquilibriumData
        data for equilibrium
    strategy: int
        number of stategy to use (0 to 6)
        defaults to 4 (analytical_fraction_low)
    """
    def __strategy_flat(filln):
        return np.full((equil.n_points, equil.n_components), filln)

    def _strategy_flat_very_low() -> FArray:
        return __strategy_flat(1e-12)

    def _strategy_flat_low() -> FArray:
        return __strategy_flat(1e-9)

    def _strategy_flat_mid() -> FArray:
        return __strategy_flat(1e-6)

    def _strategy_flat_high() -> FArray:
        return __strategy_flat(1e-3)

    def __analytical_fraction(fraction: float, neg_fill: float) -> FArray:
        initial_guess = fraction*equil.analytical_concentration
        initial_guess[initial_guess<=0.0] = neg_fill
        return initial_guess

    def _analytical_fraction_low():
        return __analytical_fraction(1e-6, 1e-12)

    def _analytical_fraction_mid():
        return __analytical_fraction(1e-4, 1e-9)

    def _analytical_fraction_high():
        return __analytical_fraction(1e-2, 1e-6)

    strategies = (_strategy_flat_very_low,
                  _strategy_flat_low,
                  _strategy_flat_mid,
                  _strategy_flat_high,
                  _analytical_fraction_low,
                  _analytical_fraction_mid,
                  _analytical_fraction_high)

    return strategies[strategy]()
