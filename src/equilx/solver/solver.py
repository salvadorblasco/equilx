import numpy as np

from ..mytypes import FArray
from ..data_structure import EquilibriumData, DebyeHuckelParameters
from .precond import precondition
from .outer_fixed_point import outer_fixed_point
from .nr import newton_raphson


def EqSolver(equil: EquilibriumData, **kwargs) -> FArray:
    # Extract parameters
    if 'initial_guess' in kwargs:
        initial_guess = kwargs['initial_guess']
    else:
        strategy = kwargs.get('initial_guess_strategy', 4)
        initial_guess = guess_starting_point(equil, strategy)

    # precondition
    if kwargs.get('precondition', True):
        initial_guess[...], _ = precondition(equil, initial_guess)

    # add solids
    initial_guess = np.concatenate((initial_guess, np.full((equil.n_points, equil.n_solids), 1e-10, dtype=float)), axis=1)

    # outer point loop if necessary
    if 'dbh_pars' in kwargs:
        dbh_pars = kwargs['dbh_pars']
        assert isinstance(dbh_pars, DebyeHuckelParameters)
        res = outer_fixed_point(initial_guess, equil, dbh_pars)
    else:
        res = newton_raphson(initial_guess, equil, **kwargs)
    return res


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
