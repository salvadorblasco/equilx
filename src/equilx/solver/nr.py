from typing import Final

import numpy as np
import numpy.typing as npt

from equilx.excepts import FailedCalculateConcentrations, TooManyIterations
from equilx.data_structure import EquilibriumData
from equilx.mytypes import FArray
from .utils import species_concentration, mass_action, amatrix


def newton_raphson(x0: npt.NDArray, equil: EquilibriumData, **kwargs) -> FArray:
    r"""
    Solve the set of equations $J \cdot \delta c = -F$ using Newton-Raphson's method.

    Given an initial guess **x0** and the parameters of the titration, this
    function uses the Newton-Raphson method to solve the free
    concentrations.  It includes several
    convergence tecniques such as scaling, damping and step limiting. For the
    original source of damping and scaling see De Robertis, *et al.*,
    [^1] and also Del Piero, *et al.* [^2].


    [^1]: De Robertis, *et al.*: *Analytica Chimica Acta* 1986, **191**, 385-398
    [^2]: Del Piero, *et al.*: *Annali di Chimica* 2006, 96

    Parameters
    ----------
    x : numpy.ndarray
        The concentration array of shape (n, c+p), where n is the number of points c is the
        number of components and p is the number of precipitates.
    equil: EquilibriumData
        The equilibrium information

    Returns
    -------
    x : numpy.ndarray
        Array of shape (n, c+p) of final concentrations of all components in the system.

    Raises
    ------
    ValueError
        If the input is incorrect.
    FailedCalculateConcentrations
        If the Jacobian or residuals cannot be calculated.
    TooManyIterations
        If too many iterations are performed without convergence.
    """
    OVERSHOOT_CAPPING_LOG: Final[float] = kwargs.get('overshoot_capping_log', 10.0)
    OVERSHOOT_CAPPING_LIN: Final[float] = kwargs.get('overshoot_capping_lin', 100.0)
    MAX_ITERATIONS: Final[int] = kwargs.get('max_iterations', 100)
    ZERO_OFFDIAG: Final[bool] = kwargs.get('zero_offdiag', False)
    THRESHOLD_RESIDUAL: Final[float] = kwargs.get('threshold_residual', 1e-20)
    THRESHOLD_STEP: Final[float] = kwargs.get('threshold_step', 1e-20)
    DEBUG: Final[bool] = kwargs.get('debug', False)
    x: FArray = np.log(x0)
    scomps = slice(0, equil.n_components)
    ssolid = slice(equil.n_components, equil.n_components+equil.n_solids)
    lgsi: FArray = np.array([], dtype=float)
    x[:,ssolid] = x0[:,ssolid]
    breakpoint()

    # ------ main loop ----------
    for iteration in range(MAX_ITERATIONS):
        c_components = np.exp(x[:, scomps])
        c_solids = x[:, ssolid]
        c_soluble = species_concentration(c_components, equil, full=True)

        if equil.n_solids:
            lgsi: FArray = x[:, scomps] @ equil.solid_stoichiometry.T
            assert lgsi.shape == (equil.n_points, equil.n_solids)

        F = objective_function(equil, c_soluble, c_solids, lgsi)
        J = jacobian(c_soluble, c_solids, equil)

        if np.any(np.isnan(J)):
            msg2 = f"could not calculate jacobian (iteration {iteration + 1})"
            raise FailedCalculateConcentrations(msg2, x)
        if np.any(np.isnan(F)):
            msg2 = f"could not calculate residuals (iteration {iteration + 1})"
            raise FailedCalculateConcentrations(msg2, x)

        if ZERO_OFFDIAG:
            J *= np.eye(equil.n_components+equil.n_solids)[None,...]  # zerom

        dx = -np.squeeze(np.linalg.inv(J) @ F[...,None])

        if DEBUG:
            culprit = tuple(map(int, np.unravel_index(np.argmax(np.abs(F)), F.shape)))
            pf = float(np.std(F[:,scomps]))
            ps = float(np.std(F[:,ssolid]))
            print(f"{iteration=:<3} |max F|={np.max(np.abs(F)):.6e}, |max dx|={np.max(np.abs(dx)):.6e} {pf=:.4e} {ps=:.4e} {culprit=}")

        check_residuals = np.all(np.abs(F) < THRESHOLD_RESIDUAL)
        check_step = np.all(np.abs(dx) < THRESHOLD_STEP)

        if check_step or check_residuals:
            return np.clip(_to_linearbase(x, ssolid), a_min=0.0, a_max=None)

        dx[:,scomps] = np.clip(dx[:,scomps], a_min=-OVERSHOOT_CAPPING_LOG, a_max=OVERSHOOT_CAPPING_LOG)
        dx[:,ssolid] = np.clip(dx[:,ssolid], a_min=-OVERSHOOT_CAPPING_LIN, a_max=OVERSHOOT_CAPPING_LIN)

        x += dx

        # if equil.n_solids:
        #     remove_solid = np.logical_and(F[:,ssolid] < 0.0, dx[:, ssolid] < 0.0)
        #     x[:,ssolid][remove_solid] = -49.0

        # correct for solids absent
        # if equil.n_solids > 0:
        #     negative_s = x[:,ssolid] < 0.0
        #     x[:,ssolid][negative_s] = 0.0

    raise TooManyIterations(f"The limit of {MAX_ITERATIONS} has been reached", np.exp(x))


def _to_linearbase(array, solid_slice):
    retval = np.exp(array)
    retval[:,solid_slice] = array[:,solid_slice]
    return retval


def objective_function(equil: EquilibriumData, free_concentration: FArray,
                       solid_concentration: FArray, lgsi: FArray) -> FArray:
    """
    Calculate the objective function for a given set of parameters.

    Parameters
    ----------
    concentration : numpy.ndarray
        The concentration array of shape (n, c+p), where n is the number of points c is the number
        of components and p is the number of solid species.
    total_concentration : numpy.ndarray
        The total concentration vector with shape (n,c), where n is the number of points and c is
        the number of components.

    Returns
    -------
    delta : numpy.ndarray
        The objective function values with shape (n, nt), where nt is the total number of
        components (n + f).
    """
    ss_positive = _get_oversaturated(solid_concentration)
    ffunc = mass_action(equil, free_concentration, ss_positive)
    if equil.n_solids:
        ss_negative = _get_undersaturated(solid_concentration)
        gfunc = lgsi - equil.log_ks - ss_negative
        return np.concatenate((ffunc, gfunc), axis=1)
    return ffunc


def _get_oversaturated(array):
    return np.where(array>0.0, array, np.zeros_like(array))


def _get_undersaturated(array):
    return np.where(array<0.0, array, np.zeros_like(array))


def jacobian(free_concentration: FArray, solid_concentration: FArray, equil: EquilibriumData) -> FArray:
    """
    Compute the Jacobian matrix for the given system of equations.

    Parameters
    ----------
    concentration : numpy.ndarray
        The concentration array of shape (n, c+p), where n is the number of points c is the number
        of components and p is the number of solid species.
    stoichiometry : numpy.ndarray
        The stoichiometry matrix of shape (c, s), representing the stoichiometric coefficients of
        the soluble components.
    solid_stoichiometry : numpy.ndarray
        The solid stoichiometry matrix of shape (n, p), representing the stoichiometric coefficients
        of the solid components.

    Returns
    -------
    numpy.ndarray
        The Jacobian matrix
    """
    # dim 1 -> number of titration points
    # dim 2 -> function f1, f2, ... g1, g2 ...
    # dim 3 -> variable c1, c2, ... s1, s2 ...
    dfdc = amatrix(free_concentration, equil.stoichiometryx)
    assert dfdc.shape == (equil.n_points, equil.n_components, equil.n_components)
    if equil.n_solids:
        mask_positive = _get_oversaturated(solid_concentration)
        mask_positive[mask_positive>0.0] = 1.0
        dfds = equil.solid_stoichiometry.T[None,...] *  mask_positive[:,None,:]
        assert dfds.shape == (equil.n_points, equil.n_components, equil.n_solids)

        dgdc = equil.solid_stoichiometry[None, ...]
        assert dgdc.shape == (1,  equil.n_solids, equil.n_components)

        mask_negative = _get_undersaturated(solid_concentration)
        mask_negative[mask_negative<0.0] = 1.0
        dgds = np.repeat(-np.eye(equil.n_solids)[None,...], equil.n_points, axis=0) * mask_negative[:, None, :]
        assert dgds.shape == (equil.n_points,  equil.n_solids, equil.n_solids)

        return np.block([[dfdc, dfds],
                         [np.repeat(dgdc, equil.n_points, axis=0), dgds]])
    return dfdc
