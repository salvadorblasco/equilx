"Module for solving equilibrium problems."

from typing import Optional

import numpy as np

from equilx.mytypes import FArray
from equilx.data_structure import EquilibriumData, DebyeHuckelParameters

from equilx.solver import fsolver


def Eq(equil: EquilibriumData,
       initial_guess: Optional[FArray] = None,
       debye_params: Optional[DebyeHuckelParameters] = None,
       freeze: Optional[int] = -1, **kwargs) -> dict:
    "Main function for equilibrium solving."

    retdict = {}
    control = 0

    funcs = {
         (False, False): fsolver.solve,
         (True,  False): fsolver.solve_solids,
         (False, True):  fsolver.solve_vis,
         (True,  True):  fsolver.solve_solids_vis
    }
    args = []

    with_solids = equil.any_solid()
    variable_ionic = (debye_params is not None)

    # Extract parameters
    if not initial_guess:
        strategy = kwargs.get('initial_guess_strategy', 4)
        initial_guess = guess_starting_point(equil, strategy)

    free_conc = np.zeros((equil.n_points, equil.n_soluble), order='F')
    free_conc[:,:equil.n_components] = initial_guess

    # prepare arguments
    args.append(free_conc)
    args.append(np.asfortranarray(equil.analytical_concentration))
    args.append(np.asfortranarray(10**equil.log_beta))
    args.append(np.asfortranarray(equil.stoichiometry))

    if with_solids:
        args.append(np.zeros((equil.n_points, equil.n_solids), order='F'))
        args.append(np.asfortranarray(equil.solid_stoichiometry))
        args.append(np.asfortranarray(10**equil.log_ks))

    if variable_ionic:
        args.append(equil.reference_ionic)
        args.append(np.asfortranarray(equil.background_ions))
        args.append(np.asfortranarray(equil.chargesx))
        args.append(debye_params.a)
        args.append(debye_params.b)
        args.append(np.zeros((equil.n_equils, 3), order='F'))

    args.append(control)
    #args.append(equil.n_points)
    #args.append(equil.n_components)
    #args.append(equil.n_equils)
    print("NPT=", equil.n_points)     
    print("NCP=", equil.n_components) 
    print("NEQ=", equil.n_equils)     

    # solve           (C, T, B,    P,                                                CTRL, NPT, NCP, NEQ)
    # solve_vis       (C, T, Bref, P,            Iref, Iback, charge, DHA, DHB, DHP, CTRL, NPT, NCP, NEQ)
    # solve_solids    (C, T, B,    P, CS, KS, Q,                                     CTRL, NPT, NCP, NEQ, NSL)
    # solve_solids_vis(C, T, Bref, P, CS, KS, Q, Iref, Iback, charge, DHA, DHB, DHP, CTRL, NPT, NCP, NEQ, NSL)
    
    caller = funcs[(with_solids, variable_ionic)]
    caller(*args)

    retdict['free concentration'] = args[0]
    print(args[0])
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
