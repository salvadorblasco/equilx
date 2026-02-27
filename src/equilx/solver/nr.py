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
        The concentration array of shape (n, c+p), where n is the number of points c is the number of components and p is the number of precipitates.
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
    OVERSHOOT_CAPPING: Final[float] = 1.0
    MAX_ITERATIONS: Final[int] = kwargs.get('max_iterations', 100)
    ZERO_OFFDIAG: Final[bool] = kwargs.get('zero_offdiag', False)
    THRESHOLD_RESIDUAL: Final[float] = kwargs.get('threshold_residual', 1e-10)
    THRESHOLD_STEP: Final[float] = kwargs.get('threshold_step', 1e-10)
    DEBUG = False
    x = np.log(x0)
    scomps = slice(0, equil.n_components)
    ssolid = slice(equil.n_components, equil.n_components+equil.n_solids)
    x[:,ssolid] = x0[:,ssolid]

    # ------ main loop ----------
    for iteration in range(MAX_ITERATIONS):
        c_components = np.exp(x[:, scomps])
        c_solids = x[:, ssolid]
        c_soluble = species_concentration(c_components, equil, full=True)

        if equil.n_solids:
            si = np.prod(c_components[...,None]**equil.solid_stoichiometry.T[None,...],axis=1) / 10**equil.log_ks
        else:
            si = np.array([])

        F = objective_function(equil, c_soluble, c_solids, si)
        J = jacobian(c_soluble, c_solids, equil, si)

        if np.any(np.isnan(J)):
            msg2 = f"could not calculate jacobian (iteration {iteration + 1})"
            raise FailedCalculateConcentrations(msg2, x)
        if np.any(np.isnan(F)):
            msg2 = f"could not calculate residuals (iteration {iteration + 1})"
            raise FailedCalculateConcentrations(msg2, x)

        if ZERO_OFFDIAG:
            J *= np.eye(equil.n_components+equil.n_solids)[None,...]  # zerom

        dx = -np.squeeze(np.linalg.pinv(J) @ F[...,None])

        if DEBUG:
            culprit = tuple(map(int, np.unravel_index(np.argmax(np.abs(F)), F.shape)))
            pf = float(np.std(F[:,scomps]))
            ps = float(np.std(F[:,ssolid]))
            print(f"{iteration=} |max F|={np.max(np.abs(F))}, |max dx|={np.max(np.abs(dx))} {pf=} {ps=} {culprit=}")

        check_residuals = np.all(np.abs(F) < THRESHOLD_RESIDUAL)
        check_step = np.all(np.abs(dx) < THRESHOLD_STEP)

        if check_step or check_residuals:
            return _to_linearbase(x, ssolid)

        # cap step for soluble species
        overshoots = np.abs(dx[:,scomps]) > OVERSHOOT_CAPPING
        if np.any(overshoots):
            dx[:,scomps][overshoots] = np.sign(dx[:,scomps][overshoots]) * OVERSHOOT_CAPPING

        x += dx
        # correct for solids absent
        if equil.n_solids > 0:
            negative_s = x[:,ssolid] < 0.0
            x[:,ssolid][negative_s] = 0.0

    raise TooManyIterations(f"The limit of {MAX_ITERATIONS} has been reached",
                            _to_linearbase(x, ssolid))


def _to_linearbase(array, solid_slice):
    retval = np.exp(array)
    retval[:,solid_slice] = array[:,solid_slice]
    return retval


def objective_function(equil: EquilibriumData, free_concentration: FArray, solid_concentration: FArray, si: FArray) -> FArray:
    """
    Calculate the objective function for a given set of parameters.

    Parameters
    ----------
    concentration : numpy.ndarray
        The concentration array of shape (n, c+p), where n is the number of points c is the number of components and p is the number of solid species.
    log_ks : numpy.ndarray
        The logarithm of the solubility products with shape (n, p), where p is the number of solid species.
    stoichiometry : numpy.ndarray
        The stoichiometry matrix with shape (n, s), where s is the number of soluble species.
    solid_stoichiometry : numpy.ndarray
        The stoichiometry matrix for solid species with shape (n, p), where p is the number of solid species.
    total_concentration : numpy.ndarray
        The total concentration vector with shape (n,c), where n is the number of points and c is the number of components.

    Returns
    -------
    delta : numpy.ndarray
        The objective function values with shape (n, nt), where nt is the total number of components (n + f).
    """
    ffunc = mass_action(equil, free_concentration, solid_concentration)
    if equil.n_solids:
        free = free_concentration[:,:equil.n_components]
        si = np.prod(free[...,None]**equil.solid_stoichiometry.T[None,...],axis=1) / 10**equil.log_ks
        #print(si[25:39].T)
        gfunc = 1.0 - si
        # gfunc[si<=1.0] = 0.0
        # gfunc[solid_concentration==0.0] = 0.0
        gfunc[np.logical_and(solid_concentration==0.0, si <= 1.0)] = 0.0
        return np.concatenate((ffunc, gfunc), axis=1)
    else:
        return ffunc


def jacobian(free_concentration: FArray, solid_concentration: FArray, equil: EquilibriumData, si: FArray):
    """
    Compute the Jacobian matrix for the given system of equations.

    Parameters
    ----------
    concentration : numpy.ndarray
        The concentration array of shape (n, c+p), where n is the number of points c is the number of components and p is the number of solid species.
    stoichiometry : numpy.ndarray
        The stoichiometry matrix of shape (c, s), representing the stoichiometric coefficients of the soluble components.
    solid_stoichiometry : numpy.ndarray
        The solid stoichiometry matrix of shape (n, p), representing the stoichiometric coefficients of the solid components.

    Returns
    -------
    numpy.ndarray
        The Jacobian matrix of shape (n, nt, nt), where nt is the total number of components (n + f).
    """
    # dim 1 -> number of titration points
    # dim 2 -> function f1, f2, ... g1, g2 ...
    # dim 3 -> variable c1, c2, ... s1, s2 ...
    dfdc = amatrix(free_concentration, equil.stoichiometryx)
    if equil.n_solids:
        dfds = np.where(solid_concentration[:,None,:]==0.0,
                        np.zeros((equil.n_points, equil.n_components, equil.n_solids), dtype=float), 
                        np.ones((equil.n_points,1,1), dtype=float)*equil.solid_stoichiometry.T)
        dgdc = -equil.solid_stoichiometry[None,...] * si[...,None]
        dgdc[np.logical_and(solid_concentration==0.0, si <= 1.0)] = 0.0
        dgds = np.zeros((equil.n_points, equil.n_solids, equil.n_solids))
        return np.block([[dfdc, dfds],[dgdc, dgds]])
    return dfdc
