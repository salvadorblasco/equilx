import numpy as np

from ..mytypes import FArray, IArray
from ..data_structure import EquilibriumData


def freeze_concentration(equil: EquilibriumData, reference_component: int) -> EquilibriumData:
    ...


def species_concentration(concentration, equil: EquilibriumData, full=False) -> FArray:
    r"""Calculate the species concentrations through the mass action law.

    $$
    S_{i} = \beta_i \prod_{j=1}^{N_c} C_j^{p_{ij}}
    $$

    With $S_i$ being the concentration of the species $i$, $\beta_i$ the equilibrium constant of the species $i$,
    $C_j$ the concentration of the component $j$, and $p_{ij}$ the stoichiometric coefficient of the component $j$ in the species $i$.

    Parameters
    ----------
    concentration : numpy.ndarray
        The concentration array of shape (n, c+p), where n is the number of points c is the number of components and p is the number of solid species.
    equil : EquilibriumData
        Data for equilibrium
    full : bool, optional
        If True, return the concentrations of all species including the original concentrations.
        If False, return only the concentrations of the new species.

    Returns
    -------
    numpy.ndarray
        The calculated species concentrations.
    """
    cext = 10**equil.log_beta*np.prod(np.power(concentration[...,None], equil.stoichiometry.T[None,...]), axis=1) 
    if full:
        return np.concatenate((concentration, cext), axis=1)
    return cext


def percent_distribution(concentration: FArray, equil: EquilibriumData, reference: int) -> FArray:
    """Transform free concentrations to relative concentrations.

    This function converts absolute concentrations to relative percent of
    concentrations with respect to a given species.

    Parameters
    ----------
    concentration: numpy.ndarray
        The extended free concentrations array.
    equil : EquilibriumData
        Data containing equilibrium information
    reference: int
        the index of the reference species

    Returns
    -------
    numpy.ndarray: The concentrations expressed as percent with
        respect to the reference species. The first dimmension is
        unchanged with respect to that of **C** but the second one is
        shorter as all components whose stoichiometric coefficient with
        respect to the reference one is zero have been removed.
    """
    stoich_ref = equil.stoichiometry[:, reference].copy()
    stoich_ref[stoich_ref == 0] = np.inf

    analytc_ref = equil.analytical_concentration[:, reference]
    new_c = concentration / analytc_ref[None, :]
    return 100*new_c*stoich_ref


def amatrix(concentration, stoichiometryx):
    r"""Calculate the matrix **A**.

    **A** is a matrix that apperars commonly in many equations for the elaboration
    of the jacobian. It is defined  as follows:

    .. math::

        A_{nij} = c_{nk}\delta_{ki} + \sum_{j=1}^E p_{ki}p_{jk}c_{n,j+S}
    """
    return np.einsum('ji,jk,...j->...ik', stoichiometryx, stoichiometryx, concentration)


def calc_saturation(concentration: FArray, solid_stoichiometry: IArray, ks: FArray) -> FArray:
    return np.prod(concentration[...,None]**solid_stoichiometry[...,None],axis=1) / ks 


def make_titration(initial_conc: FArray, buret_conc: FArray, starting_vol: float, titre: FArray) -> FArray:
    ...


def make_distribution(initial_conc: FArray, final_conc: FArray, px: list[bool], n_points: int) -> FArray:
    assert len(initial_conc) == len(final_conc)
    retval = np.empty((n_points, len(initial_conc)), dtype=float)
    log_which = np.array(px, dtype=bool)
    lin_which = np.logical_not(log_which)

    retval[:, lin_which] = np.linspace(initial_conc[lin_which], final_conc[lin_which], n_points, endpoint=True)
    retval[:, log_which] = np.logspace(initial_conc[log_which], final_conc[log_which], n_points, endpoint=True)
    return retval


def mass_action(equil: EquilibriumData, free_concentration: FArray, solid_concentration: FArray|None=None) -> FArray:
    ffunc = free_concentration @ equil.stoichiometryx - equil.analytical_concentration
    if solid_concentration is not None: # solids present
        ffunc += solid_concentration @ equil.solid_stoichiometry
    return ffunc
