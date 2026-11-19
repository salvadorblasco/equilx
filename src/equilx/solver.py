"Module for solving equilibrium problems."

from typing import Optional

import numpy as np

from equilx.mytypes import FArray
from equilx.data_structure import EquilibriumData, DebyeHuckelParameters

from equilx import fsolver
from equilx import utils


def Eq(equil: EquilibriumData,
       initial_guess: FArray | None = None,
       debye_params: DebyeHuckelParameters|None = None,
       freeze: int | None = None, **kwargs) -> dict:
    "Main function for equilibrium solving."

    retdict = {}
    control = 0
    #breakpoint()
    funcs = {
         (False, False, False): fsolver.eq.eq_000,
         (False, True,  False): fsolver.eq.eq_0s0,
         (False, False, True):  fsolver.eq.eq_00i,
         (False, True,  True):  fsolver.eq.eq_0si,
         (True,  False, False): fsolver.eq.eq_f00,
         (True,  True,  False): fsolver.eq.eq_fs0,
         (True,  False, True):  fsolver.eq.eq_f0i,
         (True,  True,  True):  fsolver.eq.eq_fsi
    }
    args = []

    with_solids = equil.any_solid()
    variable_ionic = debye_params is not None
    with_frozen = freeze is not None
    if with_frozen and (freeze < 1 or freeze > equil.n_components):
        raise ValueError("component index out of range")

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

    if freeze is not None:
        args.append(freeze)

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

    caller = funcs[(with_frozen, with_solids, variable_ionic)]
    caller(*args)

    retdict['free concentration'] = args[0]

    if with_solids:
        retdict['solid concentration'] = args[4]
        retdict['saturation'] = utils.calc_saturation(args[0], equil)

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
