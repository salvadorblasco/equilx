"Fixtures for ionic strength fixing"

import numpy as np

from equilx.data_structure import EquilibriumData
from equilx.solver.utils import make_titration

def case_phosphoric():
    labels = ('H', 'FOS')
    log_beta = np.array([12.38, 19.578, 21.726, -13.73])
    stoich = np.array([[1,1],[1,2],[1,3],[0,-1]], dtype=int)
    charge = np.array([1,-3], dtype=int)
    background_ions = 0.1
    titre_ions = 0.1
    return stoich, log_beta


def large_fixture1():
    labels = "L Cu Im H".split()
    n_comps = len(labels)
    logbeta = np.fromiter(map(float, "10.66 20.07 28.31 35.67 42.57 48.44 7.09 4.2 7.73 10.6 13   44.36 40.54 34.41 26.97 17.21 25.18 17.97 50.43 45.82 38.97 36.3 29.79 22.28 -13.77".split()), dtype=float)
    n_equil = len(logbeta)
    stoich = np.fromiter(map(float, "1 0 0 1 1 0 0 2 1 0 0 3 1 0 0 4 1 0 0 5 1 0 0 6 0 0 1 1 0 1 1 0 0 1 2 0 0 1 3 0 0 1 4 0 1 1 0 4 1 1 0 3 1 1 0 2 1 1 0 1 1 1 0 0 1 2 0 0 1 2 0 -1 1 1 1 4 1 1 1 3 1 1 1 2 1 1 1 1 1 2 1 0 1 2 1 -1 0 0 0 -1".split()), dtype=int).reshape(n_equil, n_comps)
    # analc = make_distribution(np.array([0.001, 0.001, 0.001, 2.0]),
    #                           np.array([0.001, 0.001, 0.001, 11.0]),
    #                           [False, False, False, True], 10)
    analc = np.empty((20,4), dtype=float)
    analc[:,0:3] = 0.001
    analc[:,3] = np.linspace(0.005, -0.005, 20)
    charge = np.array([0,2,-1,1], dtype=int)
    return EquilibriumData(analc, logbeta, stoich, charges=charge)

def case_aluminium():
    analc = make_titration(initial_conc=np.array([1e-8, 5e-5]),
                           buret_conc=np.array([0.0, -0.0001]),
                           starting_vol=1.0,
                           titre=np.linspace(0.0, 0.98, 50))
    stoich = np.array([[1, -1],
                       [1, -2],
                       [1, -3],
                       [1, -4],
                       [0, -1]], dtype=int)
    log_beta = np.array([-5.0, -9.3, -15.0, -23.0, -14.0])
    sstoich = np.array([[1,-3]], dtype=int)
    log_ks = np.array([5.5])
    charges = np.array([3, 1], dtype=int)
    eq = EquilibriumData(
        analytical_concentration=analc,
        log_beta=log_beta,
        stoichiometry=stoich,
        log_ks=log_ks,
        solid_stoichiometry=sstoich,
        charges=charges)
    return eq
