import numpy as np

from ..mytypes import FArray
from ..data_structure import EquilibriumData


def unpack_eqsolve(concentration: FArray, equil: EquilibriumData) -> dict:
    c_components = concentration[:, :equil.n_components]
    c_extended = species_concentration(c_components, equil)
    if equil.n_solids:
        c_solids = concentration[:, equil.n_components:]
        saturation = calc_saturation(c_components, equil)
    else:
        c_solids = None
        saturation = None
    return {
        'concentration components': c_components,
        'concentration species': c_extended,
        'concentration solids': c_solids,
        'concentration full': np.hstack((c_components, c_extended)),
        'saturation index': saturation
    }


def unfreeze_concentration(original_equil: EquilibriumData, frozen_equil: EquilibriumData,
                           free_concentration: FArray, reference_component: int) -> EquilibriumData:
    ...


def freeze_concentration(equil: EquilibriumData, reference_component: int) -> EquilibriumData:
    r"""Convert one component to independent variable.

    When solving the equilibrium, sometimes those concentrations are plotted
    as function of the concentration of one of them, typically the pH. That
    component, therefore, can be converted into an independent variable
    and removed from the unknowns.

    .. math::

    c_{i+S} = \beta_ic_{\omega}^{p_{i\omega}}\prod_{j\ne\omega}^Sc_j^{p_{ij}}
            = \beta_i'\prod_{j\ne\omega}^Sc_j^{p_{ij}}

    Parameters
    ----------
    equil: EquilibriumData
        The data containing equilibrium information
    reference_component: int
        The index of the reference component

    Returns
    -------
    EquilibriumData
        The new frozen equilibrium data.
    """
    new_x = equil.analytical_concentration[:, reference_component]
    beta = 10**equil.log_beta
    beta_prime = beta[np.newaxis, ...] * \
        new_x[..., np.newaxis]**equil.stoichiometry[:, reference_component]
    stoich_new = np.delete(equil.stoichiometry, reference_component, axis=1)
    analc_new = np.delete(equil.analytical_concentration, reference_component, axis=1)
    # TODO add solids as well.
    frozen_equil = EquilibriumData(
        log_beta=np.log10(beta_prime),
        stoichiometry=stoich_new,
        analytical_concentration=analc_new)
    return frozen_equil


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


def calc_saturation(concentration: FArray, equil: EquilibriumData) -> FArray:
    return np.prod(concentration[...,None]**equil.solid_stoichiometry.T[None,...],axis=1) / (10**equil.log_ks) 


def make_titration(initial_conc: FArray, buret_conc: FArray, starting_vol: float, titre: FArray) -> FArray:
    moles = initial_conc[None,:]*starting_vol +  buret_conc[None,:]*titre[:,None]
    conc = moles / (titre[:, None] + starting_vol)
    return conc


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
    if equil.n_solids and solid_concentration is not None: # solids present
        ffunc += solid_concentration @ equil.solid_stoichiometry
    return ffunc
