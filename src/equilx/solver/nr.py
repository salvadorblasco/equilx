from typing import Any, Dict, List

import numpy as np
import numpy.typing as npt

from libeq.excepts import FailedCalculateConcentrations, TooManyIterations
from libeq.utils import species_concentration, amatrix, _extend_stoichiometry

FArray = npt.NDArray[float]
IArray = npt.NDArray[int]


def newton_raphson(
    x0: npt.NDArray,
    log_beta: npt.NDArray,
    log_ks: npt.NDArray,
    stoichiometry: npt.NDArray,
    solid_stoichiometry: npt.NDArray,
    total_concentration: npt.NDArray,
    solids_idx: List[int] = [],
    max_iterations: int = 200,
    threshold_residual: float = 1e-10,
    threshold_step: float = 1e-8,
    damping: bool = False,
    forcer: bool = False,
    scaling: bool = False,
    step_limiter: bool = False,
    zero_offdiag: bool = False,
    debug: bool = False,
    panic: bool = True,
    **kwargs: Dict[str, Any],
) -> npt.NDArray:
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
    log_beta : numpy.ndarray
        The logarithm of the equilibrium constants with shape (n, s), where s is the number of soluble species.
    log_ks : numpy.ndarray
        The logarithm of the solubility products with shape (n, p), where p is the number of solid species.
    stoichiometry : numpy.ndarray
        The stoichiometry matrix.
    solid_stoichiometry : numpy.ndarray
        The stoichiometry matrix with shape (n, p), where s is the number of precipitable species.
    total_concentration : numpy.ndarray
        The total concentration vector with shape (n, c), where n is the number of points c is the number of components..
    solids_idx : List[int]
        Indices of solid species to compute.
    max_iterations : int
        Maximum number of iterations.
    threshold : float
        Convergence threshold.
    damping : bool
        Whether to apply damping.
    forcer : bool
        Whether to use a line search algorithm.
    scaling : bool
        Whether to apply scaling to the Jacobian matrix.
    step_limiter : bool
        Whether to limit the step size.
    zero_offdiag : bool
        Whether to zero out off-diagonal elements of the Jacobian matrix.
    debug : bool
        Whether to print debug information.
    panic : bool
        Whether to save intermediate results in case of failure.
    **kwargs : Dict[str, Any]
        Additional keyword arguments.

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
    if zero_offdiag and scaling:
        raise ValueError(
            "Options scaling and zero_offdiag are not compatible with each other"
        )

    OVERSHOOT_CAPPING: float = 1.0
    n_variables = x0.shape[1]
    n_solids = solid_stoichiometry.shape[0]
    n_components = stoichiometry.shape[1]
    stoichiometryx = _extend_stoichiometry(stoichiometry)
    ks = np.atleast_2d(10**log_ks)
    x = np.log(x0)
    scomps = slice(0, n_components)
    ssolid = slice(n_components, n_components+n_solids)
    x[:,ssolid] = x0[:,ssolid]

    # ------ main loop ----------
    for iteration in range(max_iterations):
        c_components = np.exp(x[:, scomps])
        c_solids = x[:, ssolid]
        c_soluble = species_concentration(c_components, log_beta, stoichiometry, full=True)

        F = objective_function(total_concentration,
                               c_soluble,
                               stoichiometryx,
                               c_solids,
                               solid_stoichiometry,
                               ks)
        culprit = tuple(map(int, np.unravel_index(np.argmax(np.abs(F)), F.shape)))
        pf = float(np.std(F[:,scomps]))
        ps = float(np.std(F[:,ssolid]))
        gfunc = F[:,-1, None] if n_solids == 1 else F[:, ssolid]
        J = jacobian(c_soluble, stoichiometryx, c_solids, solid_stoichiometry, gfunc)

        if np.any(np.isnan(J)):
            msg2 = f"could not calculate jacobian (iteration {iteration + 1})"
            raise FailedCalculateConcentrations(msg2, x)
        if np.any(np.isnan(F)):
            msg2 = f"could not calculate residuals (iteration {iteration + 1})"
            raise FailedCalculateConcentrations(msg2, x)

        if zero_offdiag:
            J *= np.eye(n_variables)[None,...]  # zerom

        dx = -np.squeeze(np.linalg.pinv(J) @ F[...,None])

        print(f"{iteration=} |max F|={np.max(np.abs(F))}, |max dx|={np.max(np.abs(dx))} {pf=} {ps=} {culprit=}")

        check_residuals = np.all(np.abs(F) < threshold_residual)
        check_step = np.all(np.abs(dx) < threshold_step)

        if check_step or check_residuals:
            # retval = np.exp(x)
            # retval[:,ssolid] = x[:,ssolid]
            # return retval
            return _to_linearbase(x, ssolid)

        # cap step for soluble species
        overshoots = dx[:,scomps] > OVERSHOOT_CAPPING
        if np.any(overshoots):
            dx[:,scomps][overshoots] = OVERSHOOT_CAPPING

        x += dx
        # correct for solids absent
        if n_solids > 0:
            negative_s = x[:,ssolid] < 0.0
            x[:,ssolid][negative_s] = 0.0

    raise TooManyIterations(f"The limit of {max_iterations} has been reached",
                            _to_linearbase(x, ssolid))


def _to_linearbase(array, solid_slice):
    retval = np.exp(array)
    retval[:,solid_slice] = array[:,solid_slice]
    return retval


def limit_step(x, dx):
    r"""
    Limit step.

    Given a state (x) and a step (dx), the next state is expected to be
    (x+dx). However, in some cases negative values are forbidden. This
    may happen for small values of the state. In the case of x being small
    we can approximate ($x+1 \approx e^{-x}$)

    Parameters
    ----------
    x : numpy.ndarray
        The state. It must be 1D.
    dx : numpy.ndarray
        The step. It must have the same length as x.

    Returns
    -------
    new_dx : numpy.ndarray
        The limited step.
    """
    if len(x) != len(dx):
        raise ValueError("both arguments must have the same size")
    # new = x + dx
    # who = new > 0.0
    # return np.where(who, new, x*np.exp(dx))
    one_over_del = -dx*x / (0.5 * x)
    rev_del = 1 / np.where(one_over_del > 1, one_over_del, 1)
    return dx * rev_del


def DRScaling(J, F):
    r"""Apply scaling to both jacobian and objective function.

    Applies scaling technique according to De Robertis, *et al.* [^1]
    The scaling technique overcomes the problem of divergence
    when the jacobian (here called G) matrix is near singular.
    "... scaling was applied to matrix *G* and to vector *e*
    (residuals, $e = C_{k, calcd} - C_k$
    according to the equations $g_{kj}^* = g_{kj}(g_{kk}g_{jj})^{-1/2}$
    and $(e_k^* = e_kg_{kk}^{-1/2}$ where $g_{kj}^*)$
    and $(e_k^*)$ are the elements of the scaled matrix and vector
    respectively."

    $$
    J^*_{kj} = J_{kj}(J_{kk}J_{jj})^{-1/2}\\
    F^*_{k} = F_{k}J_{kk}^{-1/2}
    $$

    [^1]: De Robertis, *et al.*: *Analytica Chimica Acta* 1986, **191**, 385-398


    Parameters
    ----------
    J : numpy.ndarray
        jacobian array, which will be modified.
        It can be of any dimensionality provided that the last two are
        of the same size.
    F : numpy.ndarray
        residuals array, which will be modified.
        If must have one dimmension less than J and the rest of the
        axes be of the same size than J.

    Returns
    -------
    d : numpy.ndarray
        The diagonal of the jacobian, $(J_{ii})$
        to scale back the result.
    """
    d = J[np.arange(J.shape[0])[:, None], np.eye(J.shape[1], dtype=bool)]
    J /= np.sqrt(d[..., np.newaxis] * d[..., np.newaxis, :])
    F /= np.sqrt(d)
    return d


def objective_function(analyticalc: FArray,
                       free_concentration: FArray,
                       stoichiometryx: IArray,
                       solid_concentration: FArray,
                       solid_stoichiometry: IArray,
                       Ks: FArray):
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
    n_components = analyticalc.shape[1]
    ffunc = free_concentration @ stoichiometryx - analyticalc
    if all(solid_stoichiometry.shape):
        # solids present
        ffunc += solid_concentration @ solid_stoichiometry
        free = free_concentration[:,:n_components]
        gfunc = np.prod(free[...,None]**solid_stoichiometry[...,None],axis=1) / Ks - 1.0
        gfunc[solid_concentration==0.0] = 0.0
        return np.concatenate((ffunc, gfunc), axis=1)
    else:
        return ffunc


def jacobian(free_concentration, stoichiometryx, solid_concentration, solid_stoichiometry, gfunc):
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
    n_solids = solid_stoichiometry.shape[0]
    n_titpoints = free_concentration.shape[0]
    dfdc = amatrix(free_concentration, stoichiometryx)
    if all(solid_stoichiometry.shape):
        # dfds = solid_concentration[:,None,:]*solid_stoichiometry[...,None]
        dfds = np.where(solid_concentration==0.0,
                        np.zeros_like(solid_concentration), 
                        np.ones_like(solid_concentration))[:,None,:]*solid_stoichiometry[...,None]
        #dfds = np.ones_like(solid_concentration)[:,None,:]*solid_stoichiometry[...,None]
        dgdc = solid_stoichiometry[None,...] * (gfunc[...,None]+1)
        dgdc[solid_concentration==0.0] = 0.0
        dgds = np.zeros((n_titpoints, n_solids, n_solids))
        return np.block([[dfdc, dfds],[dgdc, dgds]])
    else:
        return dfdc
